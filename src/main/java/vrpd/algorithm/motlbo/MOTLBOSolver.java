package vrpd.algorithm.motlbo;

import vrpd.algorithm.model.Customer;
import vrpd.algorithm.model.Evaluator;
import vrpd.algorithm.model.Solution;

import java.util.*;

/**
 * MOTLBO solver for VRPD (Multi-Objective Teaching-Learning-Based Optimization)
 *
 * Encoding: double[] individual = [K_1..K_n, D_1..D_n]
 * - K: random keys for permutation -> sort ascending -> giant tour
 * - D: drone assignment values; D_i > droneThreshold -> customer i assigned to drone
 *
 * Teacher phase:
 *  - teacher chosen from archive (nondominated set)
 *  - update each individual: X_new = X + r * (T - TF * mean)
 *
 * Learner phase:
 *  - each individual pairs with random partner and learns:
 *    if partner dominates -> X_new = X + r*(partner - X)
 *    else X_new = X + r*(X - partner)
 *
 * Selection: constrained dominance (feasible > infeasible), then Pareto dominance,
 * tie-break by sum objectives (lower better).
 *
 * Local search: lightweight 2-opt per truck route.
 *
 * Assumes Evaluator.evaluate(sol, customers) exists and sets sol.makespan, sol.carbonEmission, sol.feasible.
 */
public class MOTLBOSolver {
    private final int popSize;
    private final int maxGens;
    private final int numTrucks;
    private final Map<Integer, Customer> customers;
    private final int nCustomers; // excluding depot 0

    // encoding threshold for drone mask
    private double droneThreshold = 0.5;

    // archive size (<=0 means unlimited)
    private int archiveMax = 200;

    // local search flag
    private boolean useLocalSearch = false;

    // random seed
    private final Random rnd;

    // population: real vectors + decoded Solutions
    private List<double[]> populationX;
    private List<Solution> populationS;

    // external archive (nondominated)
    private List<Solution> archive;

    /**
     * Constructor.
     * @param popSize population size
     * @param maxGens number of generations
     * @param numTrucks number of trucks (used in decoder split)
     * @param customers map of customers (must include depot id=0)
     * @param seed random seed (use System.nanoTime() for non-reproducible runs)
     */
    public MOTLBOSolver(int popSize, int maxGens, int numTrucks, Map<Integer, Customer> customers, long seed) {
        this.popSize = popSize;
        this.maxGens = maxGens;
        this.numTrucks = numTrucks;
        this.customers = customers;
        this.rnd = new Random(seed);
        this.archiveMax = popSize;
        this.nCustomers = (int) customers.keySet().stream().filter(id -> id != 0).count();
    }

    // parameter setters
    public void setDroneThreshold(double t) { this.droneThreshold = t; }
    public void setArchiveMax(int m) { this.archiveMax = m; }
    public void setUseLocalSearch(boolean flag) { this.useLocalSearch = flag; }

    /**
     * Run MOTLBO and return external nondominated archive.
     */
    public List<Solution> run() {

        archive = new ArrayList<>();
        try {
            initPopulation();

            for (Solution s : populationS) updateArchive(s);

            for (int gen = 0; gen < maxGens; gen++) {
                // compute population mean (on real vector encoding) for teacher phase
                double[] mean = computePopulationMean(populationX);

                // choose teacher: sample from archive; if archive empty use best by sum objective from population
                double[] teacherX = selectTeacherVector();

                // Teaching Factor TF randomly chosen 1 or 2 per generation (common TLBO)
                int TF = (rnd.nextDouble() < 0.5) ? 1 : 2;

                // TEACHER PHASE: update each individual
                for (int i = 0; i < popSize; i++) {
                    double[] xi = populationX.get(i);
                    double[] xt = teacherX;
                    double[] xNew = teacherPhaseUpdate(xi, xt, mean, TF);
                    // decode/eval/select
                    Solution sNew = decodeToSolution(xNew);
                    repairSolution(sNew);
                    if (useLocalSearch) localSearch(sNew);
                    Evaluator.evaluate(sNew, customers);

                    Solution sCurr = populationS.get(i);
                    if (constrainedDominates(sNew, sCurr)) {
                        populationX.set(i, xNew);
                        populationS.set(i, sNew);
                    }
//                    else if (!constrainedDominates(sCurr, sNew)) {
//                        double sumNew = sNew.makespan + sNew.carbonEmission;
//                        double sumCurr = sCurr.makespan + sCurr.carbonEmission;
//                        if (sumNew < sumCurr) {
//                            populationX.set(i, xNew);
//                            populationS.set(i, sNew);
//                        }
//                    }
                    updateArchive(sNew);
                }

                // After teacher phase, optionally recompute mean or choose new teacher; we proceed to learner phase
                // LEARNER PHASE: each learner pairs with a random partner
                for (int i = 0; i < popSize; i++) {
                    int j;
                    do { j = rnd.nextInt(popSize); } while (j == i);
                    double[] xi = populationX.get(i);
                    double[] xj = populationX.get(j);

                    // decide update direction based on constrained dominance between decoded solutions
                    Solution si = populationS.get(i);
                    Solution sj = populationS.get(j);
                    double[] xNew;
                    if (constrainedDominates(sj, si)) {
                        // learn from better partner
                        xNew = learnerUpdate(xi, xj);
                    } else {
                        // explore opposite direction
                        xNew = learnerUpdate(xi, xi, xj); // xi + r*(xi - xj)
                    }

                    Solution sNew = decodeToSolution(xNew);
                    repairSolution(sNew);
                    if (useLocalSearch) localSearch(sNew);
                    Evaluator.evaluate(sNew, customers);

                    Solution sCurr = populationS.get(i);
                    if (constrainedDominates(sNew, sCurr)) {
                        populationX.set(i, xNew);
                        populationS.set(i, sNew);
                    }
//                    else if (!constrainedDominates(sCurr, sNew)) {
//                        double sumNew = sNew.makespan + sNew.carbonEmission;
//                        double sumCurr = sCurr.makespan + sCurr.carbonEmission;
//                        if (sumNew < sumCurr) {
//                            populationX.set(i, xNew);
//                            populationS.set(i, sNew);
//                        }
//                    }
                    updateArchive(sNew);
                }
                // generation done
            }
            // return archive copy
            return new ArrayList<>(archive);
        } catch (OutOfMemoryError e) {
            return new ArrayList<>(archive);
        }
    }

    /* ---------------------------
       Initialization & helpers
       --------------------------- */

    private void initPopulation() {
        populationX = new ArrayList<>(popSize);
        populationS = new ArrayList<>(popSize);

        // stable ordered list of customer ids excluding depot
        List<Integer> custIds = new ArrayList<>();
        for (int id : customers.keySet()) if (id != 0) custIds.add(id);
        Collections.sort(custIds);

        for (int i = 0; i < popSize; i++) {
            double[] ind = new double[2 * nCustomers];
            // K part random
            for (int k = 0; k < nCustomers; k++) ind[k] = rnd.nextDouble();
            // D part random
            for (int k = 0; k < nCustomers; k++) ind[nCustomers + k] = rnd.nextDouble();
            populationX.add(ind);

            Solution s = decodeToSolution(ind);
            repairSolution(s);
            if (useLocalSearch) localSearch(s);
            Evaluator.evaluate(s, customers);
            populationS.add(s);
        }
    }

    // compute mean vector across populationX
    private double[] computePopulationMean(List<double[]> popX) {
        int L = 2 * nCustomers;
        double[] mean = new double[L];
        Arrays.fill(mean, 0.0);
        for (double[] x : popX) {
            for (int j = 0; j < L; j++) mean[j] += x[j];
        }
        for (int j = 0; j < L; j++) mean[j] /= popX.size();
        return mean;
    }

    // select teacher vector: prefer random from archive; if archive empty choose best by sum objectives in population
    private double[] selectTeacherVector() {
        if (archive != null && !archive.isEmpty()) {
            // randomly choose one nondominated solution from archive
            Solution t = archive.get(rnd.nextInt(archive.size()));
            // need to find its vector representation in populationX or create by encoding from its routes/drone mask
            // easiest: try to find a matching solution in populationS by comparing routes+droneCustomers; else create an approximate vector by random
            for (int i = 0; i < populationS.size(); i++) {
                Solution s = populationS.get(i);
                if (solutionsEquivalent(s, t)) {
                    return populationX.get(i).clone();
                }
            }
            // fallback: build a pseudo-vector by mapping permutation back to keys (uniformly spaced) and drone mask from t
            return encodeSolutionToVector(t);
        } else {
            // choose best by sum-objective in population
            double bestSum = Double.POSITIVE_INFINITY;
            int bestIdx = 0;
            for (int i = 0; i < populationS.size(); i++) {
                Solution s = populationS.get(i);
                double sum = s.makespan + s.carbonEmission;
                if (sum < bestSum) { bestSum = sum; bestIdx = i; }
            }
            return populationX.get(bestIdx).clone();
        }
    }

    // teacher phase update: X_new = X + r*(T - TF * M)
    private double[] teacherPhaseUpdate(double[] x, double[] teacher, double[] mean, int TF) {
        int L = x.length;
        double[] out = new double[L];
        for (int j = 0; j < L; j++) {
            double r = rnd.nextDouble();
            double val = x[j] + r * (teacher[j] - TF * mean[j]);
            // clip to [0,1]
            if (val < 0.0) val = 0.0;
            if (val > 1.0) val = 1.0;
            out[j] = val;
        }
        return out;
    }

    // learner update: X_new = X + r*(partner - X)
    private double[] learnerUpdate(double[] x, double[] partner) {
        int L = x.length;
        double[] out = new double[L];
        for (int j = 0; j < L; j++) {
            double r = rnd.nextDouble();
            double val = x[j] + r * (partner[j] - x[j]);
            if (val < 0.0) val = 0.0;
            if (val > 1.0) val = 1.0;
            out[j] = val;
        }
        return out;
    }

    // alternate learner update for the case xi + r*(xi - xj)
    private double[] learnerUpdate(double[] xi, double[] unused, double[] xj) {
        int L = xi.length;
        double[] out = new double[L];
        for (int j = 0; j < L; j++) {
            double r = rnd.nextDouble();
            double val = xi[j] + r * (xi[j] - xj[j]);
            if (val < 0.0) val = 0.0;
            if (val > 1.0) val = 1.0;
            out[j] = val;
        }
        return out;
    }

    /* ---------------------------
       Decode/encode, repair, local search
       --------------------------- */

    // decode vector -> Solution (perm + drone mask), similar to MODE/MOEAD earlier
    private Solution decodeToSolution(double[] ind) {
        // build list of customer ids (stable order)
        List<Integer> custIds = new ArrayList<>();
        for (int id : customers.keySet()) if (id != 0) custIds.add(id);
        Collections.sort(custIds);

        List<KeyPair> keys = new ArrayList<>();
        for (int k = 0; k < nCustomers; k++) keys.add(new KeyPair(custIds.get(k), ind[k]));
        keys.sort(Comparator.comparingDouble(a -> a.key));
        List<Integer> perm = new ArrayList<>();
        for (KeyPair kp : keys) perm.add(kp.id);

        Set<Integer> droneSet = new HashSet<>();
        for (int k = 0; k < nCustomers; k++) {
            double dval = ind[nCustomers + k];
            if (dval > droneThreshold) {
                droneSet.add(custIds.get(k));
            }
        }
        Solution s = new Solution();
        s.truckRoutes = split(perm, numTrucks);
        s.droneCustomers = new HashSet<>(droneSet);
        return s;
    }

    // encode a Solution approximately back to vector: K keys from route order, D from drone set (0/1 mapping)
    private double[] encodeSolutionToVector(Solution s) {
        double[] ind = new double[2 * nCustomers];
        List<Integer> custIds = new ArrayList<>();
        for (int id : customers.keySet()) if (id != 0) custIds.add(id);
        Collections.sort(custIds);
        // map customer -> position in giant tour (0..n-1)
        Map<Integer, Integer> pos = new HashMap<>();
        List<Integer> flat = new ArrayList<>();
        for (List<Integer> r : s.truckRoutes) {
            for (int id : r) if (id != 0) flat.add(id);
        }
        for (int i = 0; i < flat.size(); i++) pos.put(flat.get(i), i);
        // assign K as normalized position
        for (int k = 0; k < nCustomers; k++) {
            int cid = custIds.get(k);
            int p = pos.getOrDefault(cid, k);
            ind[k] = (double)p / Math.max(1, nCustomers - 1);
        }
        // D part
        for (int k = 0; k < nCustomers; k++) {
            int cid = custIds.get(k);
            ind[nCustomers + k] = s.droneCustomers.contains(cid) ? 0.9 : 0.1;
        }
        return ind;
    }

    /**
     * Repair: enforce no adjacent drone-served customers within each route.
     */
    private void repairSolution(Solution sol) {
        if (sol.droneCustomers == null || sol.truckRoutes == null) return;
        for (List<Integer> route : sol.truckRoutes) {
            boolean prevIsDrone = false;
            for (int idx = 1; idx < route.size() - 1; idx++) {
                int cid = route.get(idx);
                if (sol.droneCustomers.contains(cid)) {
                    if (prevIsDrone) {
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
     * Simple local search: 2-opt per truck route
     */
    private void localSearch(Solution s) {
        if (s.truckRoutes == null) return;
        for (List<Integer> route : s.truckRoutes) {
            int m = route.size();
            if (m <= 4) continue;
            boolean improved = true;
            while (improved) {
                improved = false;
                for (int i = 1; i < m - 2; i++) {
                    for (int k = i + 1; k < m - 1; k++) {
                        double delta = twoOptDelta(route, i, k);
                        if (delta < -1e-6) {
                            reverseSublist(route, i, k);
                            improved = true;
                        }
                    }
                    if (improved) break;
                }
            }
        }
    }

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

    /* ---------------------------
       Archive & dominance helpers
       --------------------------- */

    private boolean constrainedDominates(Solution a, Solution b) {
        if (a == null || b == null) return false;
        if (a.feasible && !b.feasible) return true;
        if (!a.feasible && b.feasible) return false;
        if (a.feasible && b.feasible) {
            boolean le = a.makespan <= b.makespan && a.carbonEmission <= b.carbonEmission;
            boolean lt = a.makespan < b.makespan || a.carbonEmission < b.carbonEmission;
            return le && lt;
        }
//        double sa = a.makespan + a.carbonEmission;
//        double sb = b.makespan + b.carbonEmission;
//        return sa < sb;

        boolean betterOrEqual = a.makespan <= b.makespan && a.carbonEmission <= b.carbonEmission;
        boolean strictly = a.makespan < b.makespan || a.carbonEmission < b.carbonEmission;
        return betterOrEqual && strictly;
    }

    private void updateArchive(Solution s) {
        // discard if dominated by any
        for (Solution e : archive) {
            if (constrainedDominates(e, s)) return;
        }
        // remove dominated
        Iterator<Solution> it = archive.iterator();
        while (it.hasNext()) {
            Solution e = it.next();
            if (constrainedDominates(s, e)) it.remove();
        }
        archive.add(s.cloneSolution());
        // truncate if needed (keep best by sum objectives)
        if (archiveMax > 0 && archive.size() > archiveMax) {
            Map<Solution, Double> cd = computeCrowdingDistance(archive);
            archive.sort((a, b) -> Double.compare(cd.get(b), cd.get(a))); // sort descending by CD
            archive = new ArrayList<>(archive.subList(0, archiveMax));
        }
    }

    // helper to judge if two solutions are equivalent by simple route+drone sets
    private boolean solutionsEquivalent(Solution a, Solution b) {
        if (a == null || b == null) return false;
        if (a.truckRoutes == null || b.truckRoutes == null) return false;
        if (a.droneCustomers == null || b.droneCustomers == null) return false;
        // compare flattened sequences
        List<Integer> fa = flatten(a.truckRoutes);
        List<Integer> fb = flatten(b.truckRoutes);
        return fa.equals(fb) && a.droneCustomers.equals(b.droneCustomers);
    }

    private List<Integer> flatten(List<List<Integer>> routes) {
        List<Integer> flat = new ArrayList<>();
        for (List<Integer> r : routes) for (int id : r) if (id != 0) flat.add(id);
        return flat;
    }

    /* ---------------------------
       Utilities: split, dist, encode helper
       --------------------------- */

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

    /* ---------------------------
       Encoding helper class
       --------------------------- */
    private static class KeyPair { int id; double key; KeyPair(int id, double key) { this.id = id; this.key = key; } }
}
