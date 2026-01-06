package vrpd.algorithm.moead;

import vrpd.algorithm.model.Solution;
import vrpd.algorithm.model.Customer;
import vrpd.algorithm.model.Evaluator;

import java.util.*;

import static vrpd.algorithm.model.Evaluator.applyNormalization;

/**
 * MOEADSolver - classical MOEA/D implementation for VRPD (permutation + drone mask)
 *
 * Key features:
 *  - One subproblem per weight vector (popSize)
 *  - Tchebycheff aggregation with ideal point z
 *  - Neighborhood B(i) (size T) used for mating restriction and update
 *  - Permutation encoding for truck routes + drone mask (Set of customer ids)
 *  - Variation operators:
 *      * Order Crossover (OX) for permutation
 *      * Swap mutation for permutation
 *      * Uniform crossover + flip mutation for drone mask
 *  - Repair: enforce "no two adjacent customers in same route both drone-served"
 *  - External population (EP) maintained (nondominated)
 *
 * Usage:
 *   MOEADSolver solver = new MOEADSolver(popSize, maxGens, numTrucks, customers);
 *   List<Solution> pareto = solver.run();
 */
public class MOEADSolver {

    private final int popSize;
    private final int maxGens;
    private final int numTrucks;
    private final Map<Integer, Customer> customers;

    // MOEA/D params
    private final int T;               // neighborhood size
    private final double delta;        // prob select parents from neighborhood
    private final int nr = 2;          // max replacements in neighborhood (optional, here we replace all better)
    private final Random rnd;

    // Evolution params (sensible defaults for permutation+mask)
    private final double mutationSwapProb = 0.2;
    private final double mutationFlipProb = 0.08;
    private final double oxProb = 1.0; // always do OX when reproducing two parents

    // decomposition
    private List<double[]> weights;    // lambda vectors (size popSize)
    private List<int[]> neighborhood;  // B(i) neighbors indices

    // population: one solution per subproblem
    private List<Solution> population;     // current solution for each subproblem i
    private List<List<Integer>> populationPerms; // stored permutation (flat) corresponding to each solution (for efficiency)
    private List<Set<Integer>> populationDrones; // stored drone sets corresponding to each solution

    private double[] idealPoint;  // z*

    // external population (EP)
    private List<Solution> externalPop;

    public MOEADSolver(int popSize, int maxGens, int numTrucks, Map<Integer, Customer> customers, long seed) {
        this.popSize = popSize;
        this.maxGens = maxGens;
        this.numTrucks = numTrucks;
        this.customers = customers;
        this.rnd = new Random();
        this.T = Math.max(2, (int) Math.round(popSize * 0.1));
        this.delta = 0.9;
        initWeights();
        initNeighborhood();
    }

    /* -------------------------
       Main entry
       ------------------------- */
    public List<Solution> run() {

        externalPop = new ArrayList<>();
        try {
            initPopulation();
            updateIdealPoint();
            for (Solution s : population) updateExternalPopulation(s);

            for (int gen = 0; gen < maxGens; gen++) {
                for (int i = 0; i < popSize; i++) {
                    // 1) choose mating pool (2 parents) either from neighborhood or whole population
                    List<Integer> parentsIdx = chooseMatingPool(i);

                    // 2) reproduce -> child (perm + drone set)
                    List<Integer> parentPerm1 = populationPerms.get(parentsIdx.get(0));
                    List<Integer> parentPerm2 = populationPerms.get(parentsIdx.get(1));
                    Set<Integer> parentDrone1 = populationDrones.get(parentsIdx.get(0));
                    Set<Integer> parentDrone2 = populationDrones.get(parentsIdx.get(1));

                    List<Integer> childPerm = orderCrossover(parentPerm1, parentPerm2);
                    // mutation on permutation
                    if (rnd.nextDouble() < mutationSwapProb) swapMutate(childPerm);

                    Set<Integer> childDrone = uniformCrossoverDrone(parentDrone1, parentDrone2, childPerm);
                    flipMutateDrone(childDrone);

                    // enforce no adjacent drones in same route
                    Solution childSol = new Solution();
                    childSol.truckRoutes = split(childPerm, numTrucks);
                    childSol.droneCustomers = new HashSet<>(childDrone);
                    enforceNoAdjacentDrones(childSol);

                    // evaluate child
                    Evaluator.evaluate(childSol, customers);

                    // update ideal point
                    updateIdeal(childSol);

                    // update neighborhood B(i) using Tchebycheff
                    updateNeighborhood(i, childSol);

                    // update external population
                    updateExternalPopulation(childSol);
                }
            }

            // return EP (Pareto approximation)
            return new ArrayList<>(applyNormalization(externalPop));
        } catch (OutOfMemoryError e) {
            return new ArrayList<>(applyNormalization(externalPop));
        }
    }

    /* -------------------------
       Initialize weights and neighborhood
       ------------------------- */
    private void initWeights() {
        weights = new ArrayList<>(popSize);
        for (int i = 0; i < popSize; i++) {
            double w1 = (double) i / Math.max(1, (popSize - 1));
            double w2 = 1.0 - w1;
            weights.add(new double[]{w1, w2});
        }
    }

    private void initNeighborhood() {
        neighborhood = new ArrayList<>(popSize);
        for (int i = 0; i < popSize; i++) {
            double[] wi = weights.get(i);
            DoubleIndex[] d = new DoubleIndex[popSize];
            for (int j = 0; j < popSize; j++) {
                double[] wj = weights.get(j);
                double dist = Math.hypot(wi[0] - wj[0], wi[1] - wj[1]);
                d[j] = new DoubleIndex(dist, j);
            }
            Arrays.sort(d, Comparator.comparingDouble(a -> a.value));
            int[] nb = new int[Math.min(T, popSize)];
            for (int k = 0; k < nb.length; k++) nb[k] = d[k].index;
            neighborhood.add(nb);
        }
    }

    /* -------------------------
       Initialization of population
       ------------------------- */
    private void initPopulation() {
        population = new ArrayList<>(popSize);
        populationPerms = new ArrayList<>(popSize);
        populationDrones = new ArrayList<>(popSize);

        // list of customer ids excluding depot 0 (stable order)
        List<Integer> custList = new ArrayList<>();
        for (int id : customers.keySet()) if (id != 0) custList.add(id);

        for (int i = 0; i < popSize; i++) {
            // random permutation
            List<Integer> perm = new ArrayList<>(custList);
            Collections.shuffle(perm, new Random(System.nanoTime() + i));

            // create solution and random drone assignment (with repair)
            Solution s = new Solution();
            s.truckRoutes = split(perm, numTrucks);
            Set<Integer> droneSet = new HashSet<>();
            for (List<Integer> route : s.truckRoutes) {
                boolean prevWasDrone = false;
                for (int idx = 1; idx < route.size() - 1; idx++) {
                    int cid = route.get(idx);
                    // random choose drone for this customer with small prob,
                    // but prevent adjacency in same route
                    if (customers.get(cid).droneServe && !prevWasDrone && rnd.nextDouble() < 0.35) {
                        droneSet.add(cid);
                        prevWasDrone = true;
                    } else prevWasDrone = false;
                }
            }
            s.droneCustomers = droneSet;

            // evaluate
            Evaluator.evaluate(s, customers);

            // store both representation and solution
            population.add(s);
            populationPerms.add(new ArrayList<>(perm));
            populationDrones.add(new HashSet<>(droneSet));
        }
    }

    /* -------------------------
       Mating: choose parents indices
       ------------------------- */
    private List<Integer> chooseMatingPool(int idx) {
        int[] neigh = neighborhood.get(idx);
        List<Integer> pool = new ArrayList<>(2);
        if (rnd.nextDouble() < delta) {
            pool.add(neigh[rnd.nextInt(neigh.length)]);
            pool.add(neigh[rnd.nextInt(neigh.length)]);
        } else {
            pool.add(rnd.nextInt(popSize));
            pool.add(rnd.nextInt(popSize));
        }
        return pool;
    }

    /* -------------------------
       Variation operators for permutation + drone mask
       ------------------------- */

    // Order Crossover (OX) for permutation representation
    private List<Integer> orderCrossover(List<Integer> a, List<Integer> b) {
        int n = a.size();
        int i = rnd.nextInt(n);
        int j = rnd.nextInt(n);
        if (i > j) { int t = i; i = j; j = t; }
        Set<Integer> mid = new HashSet<>(a.subList(i, j + 1));
        List<Integer> child = new ArrayList<>(Collections.nCopies(n, -1));
        for (int k = i; k <= j; k++) child.set(k, a.get(k));
        int pos = (j + 1) % n;
        for (int k = 0; k < n; k++) {
            int idx = (j + 1 + k) % n;
            int gene = b.get(idx);
            if (!mid.contains(gene)) {
                child.set(pos, gene);
                pos = (pos + 1) % n;
            }
        }
        // safety: fill any -1 with missing genes
        if (child.contains(-1)) {
            Set<Integer> present = new HashSet<>(child);
            List<Integer> missing = new ArrayList<>();
            for (int g : a) if (!present.contains(g)) missing.add(g);
            for (int k = 0; k < n; k++) if (child.get(k) == -1) child.set(k, missing.remove(0));
        }
        return child;
    }

    private void swapMutate(List<Integer> perm) {
        int n = perm.size();
        int i = rnd.nextInt(n);
        int j = rnd.nextInt(n);
        Collections.swap(perm, i, j);
    }

    // uniform crossover for drone sets: inherit each customer's drone bit randomly
    private Set<Integer> uniformCrossoverDrone(Set<Integer> d1, Set<Integer> d2, List<Integer> perm) {
        Set<Integer> all = new HashSet<>(perm);
        Set<Integer> child = new HashSet<>();
        for (int id : all) {
            boolean in1 = d1 != null && d1.contains(id);
            boolean in2 = d2 != null && d2.contains(id);
            if (in1 && in2) child.add(id);
            else if (in1 || in2) if (rnd.nextDouble() < 0.5) child.add(id);
        }
        return child;
    }

    private void flipMutateDrone(Set<Integer> droneSet) {
        // iterate over a random ordering of customers to flip with small prob
        List<Integer> ids = new ArrayList<>(droneSet);
        // also include non-drone customers possibly to flip on
        for (int id : customers.keySet()) if (id != 0 && !droneSet.contains(id)) ids.add(id);
        Collections.shuffle(ids, rnd);
        for (int id : ids) {
            if (rnd.nextDouble() < mutationFlipProb) {
                if (droneSet.contains(id)) droneSet.remove(id);
                else droneSet.add(id);
            }
        }
    }

    /* -------------------------
       Neighborhood update using Tchebycheff
       ------------------------- */
    private void updateNeighborhood(int subIdx, Solution child) {
        int[] neigh = neighborhood.get(subIdx);
        for (int idx : neigh) {
            Solution curr = population.get(idx);
            double fChild = tchebycheff(child, weights.get(idx), idealPoint);
            double fCurr = tchebycheff(curr, weights.get(idx), idealPoint);
            if (fChild < fCurr) {
                // replace population idx with child
                population.set(idx, child.cloneSolution());
                // also update stored perm and drone set by encoding child back to vector pieces
                List<Integer> flat = flatten(child.truckRoutes);
                populationPerms.set(idx, new ArrayList<>(flat));
                populationDrones.set(idx, new HashSet<>(child.droneCustomers));
            }
        }
    }

    private double tchebycheff(Solution s, double[] w, double[] z) {
        double eps = 1e-9;
        double v1 = w[0] * Math.abs(s.makespan - z[0]);
        double v2 = w[1] * Math.abs(s.carbonEmission - z[1]);
        return Math.max(v1, v2) + eps * (s.makespan + s.carbonEmission);
    }

    /* -------------------------
       Ideal point update
       ------------------------- */
    private void updateIdealPoint() {
        idealPoint = new double[]{Double.POSITIVE_INFINITY, Double.POSITIVE_INFINITY};
        for (Solution s : population) {
            idealPoint[0] = Math.min(idealPoint[0], s.makespan);
            idealPoint[1] = Math.min(idealPoint[1], s.carbonEmission);
        }
    }
    private void updateIdeal(Solution s) {
        idealPoint[0] = Math.min(idealPoint[0], s.makespan);
        idealPoint[1] = Math.min(idealPoint[1], s.carbonEmission);
    }

    /* -------------------------
       External population (EP) maintenance (nondominated)
       ------------------------- */
    private void updateExternalPopulation(Solution cand) {
        // discard if any in EP dominates cand
        for (Solution e : externalPop) {
            if (dominates(e, cand)) return;
        }
        // remove any EP members dominated by cand
        externalPop.removeIf(e -> dominates(cand, e));
        // add clone
        externalPop.add(cand.cloneSolution());
    }

    private boolean dominates(Solution a, Solution b) {
        if (a == null || b == null) return false;
        boolean le = a.makespan <= b.makespan && a.carbonEmission <= b.carbonEmission;
        boolean lt = a.makespan < b.makespan || a.carbonEmission < b.carbonEmission;
        return le && lt;
    }

    /* -------------------------
       Repair: no adjacent drone-served customers in same route
       ------------------------- */
    private void enforceNoAdjacentDrones(Solution sol) {
        if (sol == null || sol.truckRoutes == null || sol.droneCustomers == null) return;
        sol.droneCustomers.removeIf(
                cusId -> !customers.get(cusId).droneServe
        );
        for (List<Integer> route : sol.truckRoutes) {
            boolean prevIsDrone = false;
            for (int idx = 1; idx < route.size() - 1; idx++) {
                int cid = route.get(idx);
                if (sol.droneCustomers.contains(cid)) {
                    if (prevIsDrone) {
                        // remove drone assignment to avoid adjacent drone customers
                        sol.droneCustomers.remove(cid);
                        prevIsDrone = false;
                    } else prevIsDrone = true;
                } else prevIsDrone = false;
            }
        }
    }

    /* -------------------------
       Helpers: split, flatten, distance, 2-opt (optional)
       ------------------------- */
    private List<List<Integer>> split(List<Integer> seq, int k) {
        int chunk = (int) Math.ceil((double) seq.size() / Math.max(1, k));
        List<List<Integer>> m = new ArrayList<>();
        for (int i = 0; i < k; i++) {
            int s = i * chunk, e = Math.min(s + chunk, seq.size());
            List<Integer> r = new ArrayList<>();
            r.add(0);
            if (s < e) r.addAll(seq.subList(s, e));
            r.add(0);
            m.add(r);
        }
        return m;
    }

    private List<Integer> flatten(List<List<Integer>> routes) {
        List<Integer> flat = new ArrayList<>();
        for (List<Integer> r : routes) for (int id : r) if (id != 0) flat.add(id);
        return flat;
    }

    private double dist(Customer a, Customer b) {
        if (a == null || b == null) return 0.0;
        return Math.hypot(a.x - b.x, a.y - b.y);
    }

    private double twoOptDelta(List<Integer> route, int i, int k) {
        int a = route.get(i - 1), b = route.get(i), c = route.get(k), d = route.get(k + 1);
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

    private void localSearch(Solution s) {
        if (s == null || s.truckRoutes == null) return;
        for (List<Integer> route : s.truckRoutes) {
            int m = route.size();
            if (m <= 4) continue;
            boolean improved = true;
            while (improved) {
                improved = false;
                for (int i = 1; i < m - 2; i++) {
                    for (int k = i + 1; k < m - 1; k++) {
                        double delta = twoOptDelta(route, i, k);
                        if (delta < -1e-9) {
                            reverseSublist(route, i, k);
                            improved = true;
                        }
                    }
                    if (improved) break;
                }
            }
        }
    }

    /* -------------------------
       Utility classes
       ------------------------- */
    private static class DoubleIndex { double value; int index; DoubleIndex(double v, int i) { value = v; index = i; } }
}
