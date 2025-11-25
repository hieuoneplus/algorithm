package vrpd.algorithm.mode;

import vrpd.algorithm.model.Evaluator;
import vrpd.algorithm.model.Solution;
import vrpd.algorithm.model.Customer;
import java.util.*;

import static vrpd.algorithm.model.Evaluator.applyNormalization;

/**
 * MODESolver - Multi-Objective Differential Evolution for VRPD
 *
 * Encoding: double[] individual = [K_1..K_n, D_1..D_n] where
 *  - K: random keys for permutation (sort ascending -> giant tour)
 *  - D: drone assignment values; D_i > threshold => customer i assigned to drone
 *
 * Selection: constrained-dominance (feasible > infeasible), Pareto dominance, then sum-objectives tie-breaker.
 * Uses Evaluator.evaluate(solution, customers) to compute objectives and feasibility.
 */
public class MODESolver {

    private final int popSize;
    private final int maxGens;
    private final int numTrucks;
    private final Map<Integer, Customer> customers;
    private final int nCustomers; // excluding depot (we expect ids 1..n)

    // DE parameters
    private double F = 0.5;
    private double CR = 0.6;

    // encoding threshold for drone mask
    private double droneThreshold = 0.5;

    // local search flag
    private boolean useLocalSearch = false;

    // archive size (if <=0, keep all nondominated)
    private int archiveMax = 200;

    private final Random rnd;

    // population: arrays of real vectors + decoded Solutions
    private List<double[]> populationX;
    private List<Solution> populationS;

    // external archive (nondominated)
    private List<Solution> archive;

    public MODESolver(int popSize, int maxGens, int numTrucks, Map<Integer, Customer> customers, long seed) {
        this.popSize = popSize;
        this.maxGens = maxGens;
        this.numTrucks = numTrucks;
        this.customers = customers;
        this.rnd = new Random();
        // count customers excluding depot 0
        this.nCustomers = (int) customers.keySet().stream().filter(id -> id != 0).count();
        this.archiveMax = popSize;
    }

    // setters for parameters
    public void setF(double F) { this.F = F; }
    public void setCR(double CR) { this.CR = CR; }
    public void setDroneThreshold(double t) { this.droneThreshold = t; }
    public void setUseLocalSearch(boolean v) { this.useLocalSearch = v; }
    public void setArchiveMax(int m) { this.archiveMax = m; }

    /**
     * Run MODE. Returns the external archive (approximate Pareto front).
     */
    public List<Solution> run() {
        archive = new ArrayList<>();
        try {
            initPopulation();
            for (Solution s : populationS) updateArchive(s);

            // main loop
            for (int gen = 0; gen < maxGens; gen++) {
                for (int i = 0; i < popSize; i++) {
                    // mutation (rand/1)
                    int r1, r2, r3;
                    do { r1 = rnd.nextInt(popSize); } while (r1 == i);
                    do { r2 = rnd.nextInt(popSize); } while (r2 == i || r2 == r1);
                    do { r3 = rnd.nextInt(popSize); } while (r3 == i || r3 == r1 || r3 == r2);

                    double[] xr1 = populationX.get(r1);
                    double[] xr2 = populationX.get(r2);
                    double[] xr3 = populationX.get(r3);

                    double[] v = mutateRand1(xr1, xr2, xr3);

                    // crossover
                    double[] u = crossoverBinomial(populationX.get(i), v);

                    // decode -> solution, repair, local search, evaluate
                    Solution su = decodeToSolution(u);
                    repairSolution(su); // enforce no-adjacent-drones, etc.
                    if (useLocalSearch) localSearch(su);
                    Evaluator.evaluate(su, customers);

                    // selection between su and current solution xi
                    Solution xi = populationS.get(i);
                    if (constrainedDominates(su, xi)) {
                        // replace vector and solution
                        populationX.set(i, u);
                        populationS.set(i, su);
                    }
//                    else if (!constrainedDominates(xi, su)) {
//                        // nondominated wrt each other -> tie-breaker by sum objectives (smaller wins)
//                        double sumU = su.makespan + su.carbonEmission;
//                        double sumXi = xi.makespan + xi.carbonEmission;
//                        if (sumU < sumXi) {
//                            populationX.set(i, u);
//                            populationS.set(i, su);
//                        }
//                    }

                    // update archive with su
                    updateArchive(su);
                }
            }
        } catch (OutOfMemoryError e) {
            return new ArrayList<>(applyNormalization(archive));
        }
        // return archive copy
        return new ArrayList<>(applyNormalization(archive));
    }

    /* --------------------------
       Initialization & utility
       -------------------------- */

    private void initPopulation() {
        populationX = new ArrayList<>(popSize);
        populationS = new ArrayList<>(popSize);

        List<Integer> custIds = new ArrayList<>();
        for (int id : customers.keySet()) if (id != 0) custIds.add(id);

        for (int i = 0; i < popSize; i++) {
            double[] ind = new double[2 * nCustomers];
            // K part random [0,1)
            for (int k = 0; k < nCustomers; k++) ind[k] = rnd.nextDouble();
            // D part random [0,1)
            for (int k = 0; k < nCustomers; k++) ind[nCustomers + k] = rnd.nextDouble();
            populationX.add(ind);

            // decode & evaluate initial solution
            Solution s = decodeToSolution(ind);
            repairSolution(s);
            if (useLocalSearch) localSearch(s);
            Evaluator.evaluate(s, customers);
            populationS.add(s);
        }
    }

    private double[] mutateRand1(double[] xr1, double[] xr2, double[] xr3) {
        int L = xr1.length;
        double[] v = new double[L];
        for (int j = 0; j < L; j++) {
            double val = xr1[j] + F * (xr2[j] - xr3[j]);
            // clip to [0,1]
            if (val < 0.0) val = 0.0;
            if (val > 1.0) val = 1.0;
            v[j] = val;
        }
        return v;
    }

    private double[] crossoverBinomial(double[] x, double[] v) {
        int L = x.length;
        double[] u = new double[L];
        int jrand = rnd.nextInt(L);
        for (int j = 0; j < L; j++) {
            if (rnd.nextDouble() < CR || j == jrand) u[j] = v[j];
            else u[j] = x[j];
        }
        return u;
    }

    /* --------------------------
       Decode & repair & local search
       -------------------------- */

    /**
     * Decode real vector (K+D) to Solution (routes + drone set).
     * - K part general permutation of customers 1..n
     * - D part threshold -> drone boolean assignment
     * Returned solution routes include depot 0 at start and end per route.
     */
    private Solution decodeToSolution(double[] ind) {
        // build key list
        List<KeyPair> keys = new ArrayList<>(nCustomers);
        int idx = 0;
        // map customer ids sorted by natural order (we expect customers keys maybe not contiguous; create mapping)
        List<Integer> custIds = new ArrayList<>();
        for (int id : customers.keySet()) if (id != 0) custIds.add(id);
        Collections.sort(custIds); // ensure stable order

        for (int k = 0; k < nCustomers; k++) {
            keys.add(new KeyPair(custIds.get(k), ind[k]));
        }
        // sort ascending by key value
        keys.sort(Comparator.comparingDouble(a -> a.key));

        // permutation list
        List<Integer> perm = new ArrayList<>(nCustomers);
        for (KeyPair kp : keys) perm.add(kp.id);

        // drone mask
        Set<Integer> droneSet = new HashSet<>();
        for (int k = 0; k < nCustomers; k++) {
            double dval = ind[nCustomers + k];
            if (dval > droneThreshold) {
                int cid = custIds.get(k);
                droneSet.add(cid);
            }
        }

        // split permutation into numTrucks routes (same split heuristic as MOEAD)
        List<List<Integer>> routes = split(perm, numTrucks);

        // wrap into Solution (routes contain depot 0 at start and end)
        Solution s = new Solution();
        s.truckRoutes = routes;
        s.droneCustomers = new HashSet<>(droneSet);
        return s;
    }

    /**
     * Repair: enforce no adjacent drone-served customers within each route.
     * (If adjacency occurs we remove drone assignment of latter customer.)
     * You can extend repair to address capacity, endurance etc.
     */
    private void repairSolution(Solution sol) {
        if (sol.droneCustomers == null || sol.truckRoutes == null) return;
        for (List<Integer> route : sol.truckRoutes) {
            boolean prevIsDrone = false;
            for (int idx = 1; idx < route.size() - 1; idx++) {
                int cid = route.get(idx);
                if (sol.droneCustomers.contains(cid)) {
                    if (prevIsDrone) {
                        // unassign later to avoid adjacent drone tasks
                        sol.droneCustomers.remove(cid);
                        prevIsDrone = false;
                    } else {
                        prevIsDrone = true;
                    }
                } else {
                    prevIsDrone = false;
                }
            }
        }
    }

    /**
     * Simple local search: apply 2-opt to each truck route (ignoring depot positions).
     * This is lightweight but helps solution quality.
     */
    private void localSearch(Solution s) {
        for (List<Integer> route : s.truckRoutes) {
            // route has depot at first and last positions (0 ... 0)
            int m = route.size();
            if (m <= 4) continue; // no interior edges to improve
            boolean improved = true;
            while (improved) {
                improved = false;
                // try 2-opt between interior nodes 1..m-2
                for (int i = 1; i < m - 2; i++) {
                    for (int k = i + 1; k < m - 1; k++) {
                        double delta = twoOptDelta(route, i, k);
                        if (delta < -1e-6) {
                            // perform 2-opt between positions i..k
                            reverseSublist(route, i, k);
                            improved = true;
                        }
                    }
                    if (improved) break;
                }
            }
        }
    }

    // compute change in truck distance if perform 2-opt between i..k (positions in route list)
    private double twoOptDelta(List<Integer> route, int i, int k) {
        int a = route.get(i - 1);
        int b = route.get(i);
        int c = route.get(k);
        int d = route.get(k + 1);
        Customer A = customers.get(a), B = customers.get(b), C = customers.get(c), D = customers.get(d);
        double before = dist(A, B) + dist(C, D);
        double after = dist(A, C) + dist(B, D);
        return after - before;
    }

    private void reverseSublist(List<Integer> list, int i, int k) {
        while (i < k) {
            int tmp = list.get(i);
            list.set(i, list.get(k));
            list.set(k, tmp);
            i++; k--;
        }
    }

    /* --------------------------
       Selection helpers & archive
       -------------------------- */

    private boolean constrainedDominates(Solution a, Solution b) {
        if (a == null || b == null) return false;
        if (a.feasible && !b.feasible) return true;
        if (!a.feasible && b.feasible) return false;
        // both feasible -> Pareto dominance
        if (a.feasible && b.feasible) {
            boolean le = a.makespan <= b.makespan && a.carbonEmission <= b.carbonEmission;
            boolean lt = a.makespan < b.makespan || a.carbonEmission < b.carbonEmission;
            return le && lt;
        }
        // both infeasible -> compare sum of objectives (lower is better as approximated)
//        double sa = a.makespan + a.carbonEmission;
//        double sb = b.makespan + b.carbonEmission;
//        return sa < sb;
        boolean betterOrEqual = a.makespan <= b.makespan && a.carbonEmission <= b.carbonEmission;
        boolean strictly = a.makespan < b.makespan || a.carbonEmission < b.carbonEmission;
        return betterOrEqual && strictly;
    }

    private void updateArchive(Solution s) {
        // if any in archive dominates s -> discard
        for (Solution e : archive) {
            if (constrainedDominates(e, s)) return;
        }
        // remove archive members dominated by s
        Iterator<Solution> it = archive.iterator();
        while (it.hasNext()) {
            Solution e = it.next();
            if (constrainedDominates(s, e)) it.remove();
        }
        // add clone of s (to avoid aliasing)
        archive.add(s.cloneSolution());

        // truncate archive by crowding or sum-objectives if needed
        if (archiveMax > 0 && archive.size() > archiveMax) {
            Map<Solution, Double> cd = computeCrowdingDistance(archive);
            archive.sort((a, b) -> Double.compare(cd.get(b), cd.get(a))); // sort descending by CD
            archive = new ArrayList<>(archive.subList(0, archiveMax));
        }
    }

    private Map<Solution, Double> computeCrowdingDistance(List<Solution> front) {
        int n = front.size();
        Map<Solution, Double> cd = new HashMap<>();
        if (n == 0) return cd;
        for (Solution s : front) cd.put(s, 0.0);

        // 2 objectives: makespan, carbonEmission
        // 1) sort by makespan
        front.sort(Comparator.comparingDouble(s -> s.makespan));
        cd.put(front.get(0), Double.POSITIVE_INFINITY);
        cd.put(front.get(n - 1), Double.POSITIVE_INFINITY);
        double minF = front.get(0).makespan, maxF = front.get(n - 1).makespan;
        for (int i = 1; i < n - 1; i++) {
            double prev = front.get(i - 1).makespan;
            double next = front.get(i + 1).makespan;
            double dist = (next - prev) / Math.max(1e-9, maxF - minF);
            cd.put(front.get(i), cd.get(front.get(i)) + dist);
        }

        // 2) sort by carbonEmission
        front.sort(Comparator.comparingDouble(s -> s.carbonEmission));
        cd.put(front.get(0), Double.POSITIVE_INFINITY);
        cd.put(front.get(n - 1), Double.POSITIVE_INFINITY);
        minF = front.get(0).carbonEmission;
        maxF = front.get(n - 1).carbonEmission;
        for (int i = 1; i < n - 1; i++) {
            double prev = front.get(i - 1).carbonEmission;
            double next = front.get(i + 1).carbonEmission;
            double dist = (next - prev) / Math.max(1e-9, maxF - minF);
            cd.put(front.get(i), cd.get(front.get(i)) + dist);
        }
        return cd;
    }
    /* --------------------------
       Helpers: splitting, distance, mapping
       -------------------------- */

    private List<List<Integer>> split(List<Integer> seq, int k) {
        int chunk = (int) Math.ceil((double) seq.size() / Math.max(1, k));
        List<List<Integer>> m = new ArrayList<>();
        for (int i = 0; i < k; i++) {
            int s = i * chunk, e = Math.min(s + chunk, seq.size());
            List<Integer> r = new ArrayList<>();
            r.add(0); // depot
            if (s < e) r.addAll(seq.subList(s, e));
            r.add(0);
            m.add(r);
        }
        return m;
    }

    private double dist(Customer a, Customer b) {
        if (a == null || b == null) return 0.0;
        double dx = a.x - b.x;
        double dy = a.y - b.y;
        return Math.hypot(dx, dy);
    }

    /* --------------------------
       Small utility class
       -------------------------- */
    private static class KeyPair { int id; double key; KeyPair(int id, double key) { this.id = id; this.key = key; } }

}
