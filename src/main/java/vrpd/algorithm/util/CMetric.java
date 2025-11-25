package vrpd.algorithm.util;

import vrpd.algorithm.model.Solution;
import vrpd.algorithm.model.Evaluator;

import java.util.*;


public class CMetric {
    public static double computeC(List<Solution> A, List<Solution> B, boolean onlyFeasible) {
        if (B == null || B.isEmpty()) return 0.0;
        List<Solution> Af = filterSolutions(A, onlyFeasible);
        List<Solution> Bf = filterSolutions(B, onlyFeasible);

        if (Bf.isEmpty()) return 0.0;

        int dominatedCount = 0;
        for (Solution b : Bf) {
            boolean covered = false;
            for (Solution a : Af) {
                if (constrainedDominates(a, b)) { covered = true; break; }
            }
            if (covered) dominatedCount++;
        }
        return (double) dominatedCount / (double) Bf.size();
    }

    public static Map<String, Map<String, Double>> computePairwiseCoverage(Map<String, List<Solution>> results, boolean onlyFeasible) {
        Map<String, Map<String, Double>> out = new LinkedHashMap<>();
        List<String> names = new ArrayList<>(results.keySet());
        for (String a : names) {
            out.put(a, new LinkedHashMap<>());
        }
        for (String a : names) {
            for (String b : names) {
                double c = computeC(results.get(a), results.get(b), onlyFeasible);
                out.get(a).put(b, c);
            }
        }
        return out;
    }

    /**
     * Helper: filter solutions according to feasibility or penalty cutoff.
     */
    private static List<Solution> filterSolutions(List<Solution> list, boolean onlyFeasible) {
        List<Solution> out = new ArrayList<>();
        if (list == null) return out;
        if (onlyFeasible) {
            for (Solution s : list) if (s != null && s.feasible) out.add(s);
        } else {
            double cutoff = Evaluator.INFEASIBLE_PENALTY;
            for (Solution s : list) {
                if (s == null) continue;
                // remove clearly penalized solutions (objectives >= cutoff)
                if (s.makespan < cutoff && s.carbonEmission < cutoff) out.add(s);
            }
        }
        return out;
    }

    /**
     * Constrained dominance:
     * - feasible > infeasible
     * - if both feasible: Pareto dominance (minimize)
     * - if both infeasible: compare sum of objectives (smaller better)
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
        // both infeasible -> compare sum objectives as heuristic
        double sa = a.makespan + a.carbonEmission;
        double sb = b.makespan + b.carbonEmission;
        return sa < sb;
    }

    /* --------------------
       Example printing helper
       -------------------- */
    public static void printPairwiseMatrix(Map<String, Map<String, Double>> matrix) {
        List<String> names = new ArrayList<>(matrix.keySet());
        // header
        System.out.print(String.format("%15s", ""));
        for (String b : names) System.out.print(String.format("%12s", b));
        System.out.println();
        for (String a : names) {
            System.out.print(String.format("%15s", a));
            for (String b : names) {
                Double v = matrix.get(a).get(b);
                System.out.print(String.format("%12.3f", v));
            }
            System.out.println();
        }
    }
}
