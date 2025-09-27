package vrpd.algorithm.model;

import vrpd.algorithm.model.Customer;
import vrpd.algorithm.model.Solution;

import java.util.*;

/**
 * Static evaluator compatible with MOEADSolver.evaluate(child, customers).
 *
 * Usage:
 *  - Optionally call Evaluator.setParameters(...) at startup to override defaults.
 *  - Then call Evaluator.evaluate(solution, customers).
 *
 * Notes:
 *  - Expects depot customer with id == 0 in customers map (used for coordinates).
 *  - Assumes each route in solution.truckRoutes contains 0 at start and 0 at end (MOEADSolver.split does that).
 */
public class Evaluator {

    // --- default parameters (you can override with setParameters) ---
    public static int TOTAL_EVAL = 50000;
    public static double truckCapacity = 1300.0; // default, override
    public static double droneCapacity = 10.0; // default, override
    public static double truckSpeed = 35.0; // distance units per minute
    public static double droneSpeed = 50.0; // distance units per minute
    public static double droneEndurance = 30.0; // minutes (flight+service+flight)
    public static double sLaunch = 1.0;   // launch setup (minutes)
    public static double sRecover = 1.0;  // receive setup (minutes)

    // emission params
    public static double WAER = 1.2603;   // kg CO2 per distance unit for truck
    public static double PGFER = 3.773e-4; // kg CO2 per Wh (grid)
    public static double AER = 3.3333;    // Wh per distance unit (drone energy per distance)

    // service-level target (global)
    public static double serviceLevelTarget = 0.8;

    // small penalty for infeasible (when we want to sort but keep numeric)
    public static double INFEASIBLE_PENALTY = 1e6;

    /**
     * Optional: set parameters before running MOEA/D.
     */
    public static void setParameters(int truckCapacity_,
                                     double truckSpeed_, double droneSpeed_, double droneEndurance_,
                                     double sLaunch_, double sRecover_,
                                     double WAER_, double PGFER_, double AER_,
                                     double serviceLevelTarget_) {
        truckCapacity = truckCapacity_;
        truckSpeed = truckSpeed_;
        droneSpeed = droneSpeed_;
        droneEndurance = droneEndurance_;
        sLaunch = sLaunch_;
        sRecover = sRecover_;
        WAER = WAER_;
        PGFER = PGFER_;
        AER = AER_;
        serviceLevelTarget = serviceLevelTarget_;
    }

    /**
     * Main evaluator used by MOEADSolver.
     * It mutates Solution's objective fields in-place.
     */
    public static void evaluate(Solution sol, Map<Integer, Customer> customers) {
        if(TOTAL_EVAL == 0) {
            throw new OutOfMemoryError();
        }
        // sanity: customers must contain depot id 0
        if (!customers.containsKey(0)) {
            throw new IllegalArgumentException("Customers map must contain depot with id 0");
        }

        // accumulate emission and finish times
        double totalEmission = 0.0;
        int numRoutes = Math.max(1, sol.truckRoutes.size());
        double[] truckFinishTimes = new double[numRoutes];

        // track service-level sums and counting customers once
        double sumServiceLevel = 0.0;
        Set<Integer> accounted = new HashSet<>();

        boolean capacityOK = true;
        boolean enduranceOK = true;

        // For each route
        for (int t = 0; t < sol.truckRoutes.size(); t++) {
            List<Integer> routeOriginal = sol.truckRoutes.get(t);
            if (routeOriginal == null || routeOriginal.size() == 0) {
                truckFinishTimes[t] = 0.0;
                continue;
            }

            // Build truckStops: points where truck actually stops (exclude drone-served customers)
            List<Integer> truckStops = new ArrayList<>();
            truckStops.add(0); // start depot
            for (int idx = 1; idx < routeOriginal.size() - 1; idx++) {
                int cid = routeOriginal.get(idx);
                if (!sol.droneCustomers.contains(cid)) {
                    truckStops.add(cid);
                }
            }
            truckStops.add(0); // end depot

            // Capacity check (sum demands of truck-served customers on this route)
            double load = 0;
            for (int cid : truckStops) {
                if (cid == 0) continue;
                Customer c = customers.get(cid);
                if(c != null) {
                    if (!sol.droneCustomers.contains(cid)) {
                        load += c.demand;
                    } else {
                        if(c.demand > droneCapacity) {
                            capacityOK = false;
                            break;
                        }
                    }
                }
            }
            if (load > truckCapacity) capacityOK = false;

            // Map drone jobs grouped by arc (launch#return -> list of jobs)
            Map<String, List<Integer>> droneJobsByArc = new HashMap<>();
            for (int pos = 1; pos < routeOriginal.size() - 1; pos++) {
                int j = routeOriginal.get(pos);
                if (!sol.droneCustomers.contains(j)) continue;
                Customer cj = customers.get(j);
                if (cj == null) continue;
                // find previous truck stop (or depot 0)
                int launch = 0;
                for (int p = pos - 1; p >= 0; p--) {
                    int cand = routeOriginal.get(p);
                    if (cand == 0) { launch = 0; break; }
                    if (!sol.droneCustomers.contains(cand)) { launch = cand; break; }
                }
                // find next truck stop (or depot 0)
                int ret = 0;
                for (int p = pos + 1; p < routeOriginal.size(); p++) {
                    int cand = routeOriginal.get(p);
                    if (cand == 0) { ret = 0; break; }
                    if (!sol.droneCustomers.contains(cand)) { ret = cand; break; }
                }
                String key = launch + "#" + ret;
                droneJobsByArc.computeIfAbsent(key, k -> new ArrayList<>()).add(j);
                // we'll account service-level later when scheduling
            }

            // Simulate truck movement and schedule drone missions serially across the route
            double time = 0.0; // truck time starts at 0 at depot
            // route-level latestRecovery to ensure only 1 drone active at any time on this truck
            double latestRecoveryRoute = 0.0;

            // iterate over truckStops arcs
            for (int i = 0; i < truckStops.size() - 1; i++) {
                int prev = truckStops.get(i);
                int next = truckStops.get(i + 1);
                double legDist = dist(customers.get(prev), customers.get(next));
                double truckTravel = legDist / truckSpeed;

                // key for drone jobs that launch at prev and return at next
                String arcKey = prev + "#" + next;
                List<Integer> arcDroneJobs = droneJobsByArc.getOrDefault(arcKey, Collections.emptyList());

                if (!arcDroneJobs.isEmpty()) {
                    // SERIAL scheduling: only 1 drone per truck, schedule missions one by one
                    double launchStart = time;              // earliest time we can start first setup at this stop
                    double lastLaunchStart = launchStart;   // track last actual launch start time
                    double cumulativeSetup = 0.0;          // cumulative setup time used as baseline

                    for (int jobIdx = 0; jobIdx < arcDroneJobs.size(); jobIdx++) {
                        int j = arcDroneJobs.get(jobIdx);
                        Customer cj = customers.get(j);
                        if (cj == null) continue;

                        // earliest we can start setup: either after previous setups at this stop
                        double earliestSetupStart = launchStart + cumulativeSetup;
                        // but we also must ensure drone has returned from previous missions on this route
                        double setupStart = Math.max(earliestSetupStart, latestRecoveryRoute);
                        double t_launch = setupStart + sLaunch; // actual launch time (after setup)

                        double fly1 = dist(customers.get(prev), cj) / droneSpeed;
                        double serve = cj.service;
                        double fly2 = dist(cj, customers.get(next)) / droneSpeed;

                        // endurance check: flight + serve must be <= droneEndurance
                        if (fly1 + serve + fly2 > droneEndurance + 1e-9) {
                            enduranceOK = false;
                        }

                        double droneRecovery = t_launch + fly1 + serve + fly2 + sRecover;
                        // update latest recovery time across route (serial)
                        latestRecoveryRoute = Math.max(latestRecoveryRoute, droneRecovery);

                        // service level measured at drone arrival
                        double arrivalAtJ = t_launch + fly1;
                        if (!accounted.contains(j)) {
                            double Sj = serviceLevelAt(arrivalAtJ, cj);
                            sumServiceLevel += Sj;
                            accounted.add(j);
                        }

                        // emission for drone mission (energy->CO2): PGFER * AER * distance flown
                        double uavDist = dist(customers.get(prev), cj) + dist(cj, customers.get(next));
                        totalEmission += PGFER * AER * uavDist;

                        // update cumulative setup so that next job's earliestSetupStart accounts for previous setups at same stop
                        cumulativeSetup = (t_launch - launchStart); // essentially last launch start offset
                        lastLaunchStart = t_launch;
                    }

                    // truck departs after the last launch setup has started (we model truck leaving at last launch time)
                    double truckDepart = lastLaunchStart;
                    double truckArrivalNext = truckDepart + truckTravel;
                    // truck may need to wait for the latest drone recovery
                    time = Math.max(truckArrivalNext, latestRecoveryRoute);

                    // truck emission for this leg
                    totalEmission += WAER * legDist;

                } else {
                    // No drone jobs on this leg: truck simply goes, then possibly services next
                    time += truckTravel;
                    if (next != 0) {
                        Customer cn = customers.get(next);
                        if (!accounted.contains(next)) {
                            double Sn = serviceLevelAt(time, cn);
                            sumServiceLevel += Sn;
                            accounted.add(next);
                        }
                        time += cn.service;
                    }
                    // truck emission
                    totalEmission += WAER * legDist;
                }
            } // end iterate truckStops

            truckFinishTimes[t] = time;
        } // end for each route

        // compute makespan and average service level
        double makespan = 0.0;
        for (double ft : truckFinishTimes) makespan = Math.max(makespan, ft);

        // number of customers served (expect all non-depot customers appear in routes)
        int totalCustomers = 0;
        // collect unique customers from routes
        Set<Integer> allCust = new HashSet<>();
        for (List<Integer> r : sol.truckRoutes) {
            for (int id : r) if (id != 0) allCust.add(id);
        }
        totalCustomers = allCust.size();

        double avgServiceLevel = 0.0;
        if (totalCustomers > 0) avgServiceLevel = sumServiceLevel / ((double) totalCustomers);

        boolean feasible = capacityOK && enduranceOK && (avgServiceLevel + 1e-9 >= serviceLevelTarget);

        // assign objectives back to solution (with penalty for infeasible)
        if (feasible) {
            sol.makespan = makespan;
            sol.carbonEmission = totalEmission;
        } else {
            sol.makespan = makespan + INFEASIBLE_PENALTY;
            sol.carbonEmission = totalEmission + INFEASIBLE_PENALTY;
        }
        sol.serviceLevel = avgServiceLevel;
        sol.feasible = feasible;

        sol.makespan = Math.round(sol.makespan * 100.0) / 100.0;
        sol.carbonEmission = Math.round(sol.carbonEmission * 100.0) / 100.0;
        TOTAL_EVAL--;
    }

    // Euclidean distance helper (uses Customer coords)
    private static double dist(Customer a, Customer b) {
        if (a == null || b == null) return 0.0;
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        return Math.hypot(dx, dy);
    }

    /**
     * Service level S_i(t) as in the paper:
     * - arrival < earliness: 0
     * - earliness <= arrival < twStart: 1 - ((arrival - earliness)/(twStart - earliness))^5
     * - twStart <= arrival < twEnd: 1
     * - twEnd <= arrival < lateness: 1 - ((lateness - arrival)/(lateness - twEnd))^5
     * - arrival >= lateness: 0
     */
    private static double serviceLevelAt(double arrival, Customer c) {
        if (arrival < c.earliness) return 0.0;
        if (arrival < c.twStart) {
            double denom = Math.max(1e-9, (c.twStart - c.earliness));
            double t = (arrival - c.earliness) / denom;
            return 1.0 - Math.pow(t, 5.0);
        }
        if (arrival < c.twEnd) return 1.0;
        if (arrival < c.lateness) {
            double denom = Math.max(1e-9, (c.lateness - c.twEnd));
            double t = (c.lateness - arrival) / denom;
            return 1.0 - Math.pow(t, 5.0);
        }
        return 0.0;
    }
}
