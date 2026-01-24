package vrpd.algorithm.util;

import vrpd.algorithm.model.Solution;
import vrpd.algorithm.model.Evaluator;

import java.util.*;

/**
 * IGD utility for multi-algorithm comparison.
 *
 * Usage:
 *   Map<String,List<Solution>> results = new LinkedHashMap<>();
 *   results.put("MOEAD", moeadPareto);
 *   results.put("MODE", modePareto);
 *   results.put("MOTLBO", motlboPareto);
 *   results.put("NSGAII", nsgaPareto);
 *
 *   Map<String, Double> igds = IGDUtil.computeIGDs(results, true);
 */
public class IGDUtil {

    /**
     * Compute IGD for each algorithm in the map.
     * @param algResults map algorithmName -> list of Solutions (final nondominated set from that algorithm)
     * @param onlyFeasible if true consider only solutions with sol.feasible==true; if false, filter out obviously penalized infeasible
     * @return map algorithmName -> IGD value (Double.NaN if cannot compute)
     */
    public static Map<String, Double> computeIGDs(Map<String, List<Solution>> algResults, boolean onlyFeasible) {
        // 1) Build union of all solutions (filtered)
        List<Solution> union = new ArrayList<>();
        for (List<Solution> list : algResults.values()) {
            if (list == null) continue;
            for (Solution s : list) {
                if (s == null) continue;
                if (onlyFeasible) {
                    if (s.feasible) union.add(s);
                } else {
                    // exclude clearly-penalized solutions
                    double cutoff = Evaluator.INFEASIBLE_PENALTY;
                    if (s.makespan < cutoff && s.carbonEmission < cutoff) union.add(s);
                }
            }
        }

        // if union empty -> cannot build reference set
        if (union.isEmpty()) {
            Map<String, Double> out = new LinkedHashMap<>();
            for (String name : algResults.keySet()) out.put(name, Double.NaN);
            return out;
        }

        // 2) Extract reference Pareto P* from union (non-dominated)
        List<Solution> pref = extractNonDominated(union, onlyFeasible);

        if (pref.isEmpty()) {
            Map<String, Double> out = new LinkedHashMap<>();
            for (String name : algResults.keySet()) out.put(name, Double.NaN);
            return out;
        }

        // 3) determine normalization bounds (min/max) on union (not only P*)
        double minF1 = Double.POSITIVE_INFINITY, maxF1 = Double.NEGATIVE_INFINITY;
        double minF2 = Double.POSITIVE_INFINITY, maxF2 = Double.NEGATIVE_INFINITY;
        for (Solution s : union) {
            minF1 = Math.min(minF1, s.makespan);
            maxF1 = Math.max(maxF1, s.makespan);
            minF2 = Math.min(minF2, s.carbonEmission);
            maxF2 = Math.max(maxF2, s.carbonEmission);
        }
        // avoid zero-range
        double rangeF1 = Math.max(1e-12, maxF1 - minF1);
        double rangeF2 = Math.max(1e-12, maxF2 - minF2);

        // 4) normalize P* points (list of double[2])
        List<double[]> prefNorm = new ArrayList<>(pref.size());
        for (Solution v : pref) {
            prefNorm.add(new double[]{
                    (v.makespan - minF1) / rangeF1,
                    (v.carbonEmission - minF2) / rangeF2
            });
        }

        // 5) For each algorithm compute IGD: average over v in P* of min_{u in P_alg} dist(v,u)
        Map<String, Double> igdMap = new LinkedHashMap<>();
        for (Map.Entry<String, List<Solution>> e : algResults.entrySet()) {
            String name = e.getKey();
            List<Solution> algList = e.getValue();
            // filter algList similarly to union step
            List<Solution> algFiltered = new ArrayList<>();
            if (algList != null) {
                for (Solution s : algList) {
                    if (s == null) continue;
                    if (onlyFeasible) {
                        if (s.feasible) algFiltered.add(s);
                    } else {
                        double cutoff = Evaluator.INFEASIBLE_PENALTY;
                        if (s.makespan < cutoff && s.carbonEmission < cutoff) algFiltered.add(s);
                    }
                }
            }
            // if algFiltered empty -> IGD undefined (set to +inf or NaN). We'll use Double.NaN.
            if (algFiltered.isEmpty()) {
                igdMap.put(name, Double.NaN);
                continue;
            }
            // normalize alg points
            List<double[]> algNorm = new ArrayList<>(algFiltered.size());
            for (Solution u : algFiltered) {
                algNorm.add(new double[]{
                        (u.makespan - minF1) / rangeF1,
                        (u.carbonEmission - minF2) / rangeF2
                });
            }

            // compute IGD
            double sum = 0.0;
            for (double[] v : prefNorm) {
                double minDist = Double.POSITIVE_INFINITY;
                for (double[] u : algNorm) {
                    double dx = v[0] - u[0];
                    double dy = v[1] - u[1];
                    double d = Math.hypot(dx, dy);
                    if (d < minDist) minDist = d;
                }
                sum += minDist;
            }
            double igd = sum / prefNorm.size();
            igdMap.put(name, igd);
        }

        return igdMap;
    }

    /* -------------------------
       Helpers
       ------------------------- */

    /**
     * Extract non-dominated solutions from list (using constrained dominance if onlyFeasible true)
     */
    private static List<Solution> extractNonDominated(List<Solution> list, boolean onlyFeasible) {
        List<Solution> out = new ArrayList<>();
        for (Solution s : list) {
            boolean dominated = false;
            Iterator<Solution> it = out.iterator();
            while (it.hasNext()) {
                Solution cur = it.next();
                if (constrainedDominates(cur, s)) {
                    dominated = true;
                    break;
                }
                if (constrainedDominates(s, cur)) {
                    it.remove();
                }
            }
            if (!dominated) out.add(s);
        }
        return out;
    }

    /**
     * Constrained dominance as used elsewhere:
     * - feasible > infeasible
     * - if both feasible -> Pareto dominance
     * - if both infeasible -> compare sum of objectives
     */
    private static boolean constrainedDominates(Solution a, Solution b) {
        if (a == null || b == null) return false;
        if (a.feasible && !b.feasible) return true;
        if (!a.feasible && b.feasible) return false;
        if (a.feasible && b.feasible) {
            boolean le = a.makespan <= b.makespan && a.carbonEmission <= b.carbonEmission;
            boolean lt = a.makespan < b.makespan || a.carbonEmission < b.carbonEmission;
            return le && lt;
        }
        // both infeasible -> smaller sum better
        double sa = a.makespan + a.carbonEmission;
        double sb = b.makespan + b.carbonEmission;
        return sa < sb;
    }

    /* -------------------------
       Printing helpers
       ------------------------- */

    public static void printIGDMap(Map<String, Double> igd) {
        System.out.println("IGD results:");
        for (Map.Entry<String, Double> e : igd.entrySet()) {
            Double v = e.getValue();
            if (v == null || Double.isNaN(v)) System.out.printf("%15s : %s%n", e.getKey(), "NaN");
            else System.out.printf("%15s : %.6f%n", e.getKey(), v);
        }
    }
}
