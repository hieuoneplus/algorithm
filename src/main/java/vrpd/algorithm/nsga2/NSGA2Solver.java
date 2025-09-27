package vrpd.algorithm.nsga2;

import vrpd.algorithm.model.Customer;
import vrpd.algorithm.model.Evaluator;
import vrpd.algorithm.model.Solution;

import java.util.*;

public class NSGA2Solver {
    private int popSize, gens;
    private Map<Integer, Customer> customers;
    private int numTrucks;

    private final Random rnd = new Random();

    public NSGA2Solver(int popSize, int gens, int numTrucks, Map<Integer, Customer> customers) {
        this.popSize = popSize;
        this.gens = gens;
        this.numTrucks = numTrucks;
        this.customers = customers;
    }

    public List<Solution> run() {
        List<Solution> pop = initialize();
        try {
            pop.forEach(s -> Evaluator.evaluate(s, customers));
            for (int g = 0; g < gens; g++) {
                List<Solution> off = new ArrayList<>();
                while (off.size() < popSize) {
                    Solution p1 = tournament(pop);
                    Solution p2 = tournament(pop);
                    List<Solution> kids = crossover(p1, p2);
                    for (Solution k : kids) {
                        mutate(k);
                        Evaluator.evaluate(k, customers);
                        off.add(k);
                        if (off.size() >= popSize) break;
                    }
                }
                pop.addAll(off);
                pop = selectNext(pop);

                // NSGA-II cho drone ở đây
                pop = optimizeDroneAssignment(pop, 1);
            }
            return getParetoFront(pop);
        } catch (OutOfMemoryError e) {
            return getParetoFront(pop);
        }


    }

    private List<Solution> initialize() {
        List<Solution> pop = new ArrayList<>();
        List<Integer> custs = new ArrayList<>(customers.keySet());
        custs.remove(Integer.valueOf(0));
        for (int i = 0; i < popSize; i++) {
            Collections.shuffle(custs);
            Solution sol = new Solution();
            int chunk = (int)Math.ceil((double)custs.size()/numTrucks);
            for (int t = 0; t < numTrucks; t++) {
                int start = t*chunk;
                int end = Math.min(start+chunk, custs.size());
                if (start>=end) break;
                List<Integer> route = new ArrayList<>();
                route.add(0);
                route.addAll(custs.subList(start,end));
                route.add(0);
                sol.truckRoutes.add(route);
            }
            // smart drone assignment: pick customers with largest savings
//            for (List<Integer> route : sol.truckRoutes) {
//                for (int idx = 1; idx < route.size()-1; idx++) {
//                    int id = route.get(idx);
//                    Customer prev = customers.get(route.get(idx-1));
//                    Customer next = customers.get(route.get(idx+1));
//                    Customer cur = customers.get(id);
//                    double saving = prev.distanceTo(next) - (prev.distanceTo(cur) + cur.distanceTo(next));
//                    if (saving > 0) sol.droneCustomers.add(id);
//                }
//            }
            pop.add(sol);
        }
        return pop;
    }

    private Solution tournament(List<Solution> pop) {
        Solution best = null;
        Random r = new Random();
        for (int i = 0; i < 2; i++) {
            Solution s = pop.get(r.nextInt(pop.size()));
            if (best==null || dominates(s,best)) best = s;
        }
        return best;
    }

    private List<Solution> crossover(Solution p1, Solution p2) {
        // two-point crossover on flattened routes
        List<Integer> seq1 = flatten(p1.truckRoutes);
        List<Integer> seq2 = flatten(p2.truckRoutes);
        int n = seq1.size();
        Random r = new Random();
        int i = r.nextInt(n), j = r.nextInt(n);
        if (i>j) {int tmp=i;i=j;j=tmp;}
        Set<Integer> seg = new HashSet<>(seq1.subList(i,j));
        List<Integer> child = new ArrayList<>();
        child.addAll(seq1.subList(i,j));
        for (int x : seq2) if (!seg.contains(x)) child.add(x);
        child = child.subList(0,n);
        Solution off = new Solution();
        off.truckRoutes = split(child, p1.truckRoutes.size());
        off.droneCustomers = new HashSet<>(p1.droneCustomers);
        return Arrays.asList(off);
    }

    private void mutate(Solution s) {
        Random r = new Random();
        for (List<Integer> route : s.truckRoutes) {
            if (r.nextDouble()<0.2) { // swap mutation
                int i = 1 + r.nextInt(route.size()-2);
                int j = 1 + r.nextInt(route.size()-2);
                Collections.swap(route,i,j);
            }
        }
        // flip drone assignment
//        for (Integer id : new ArrayList<>(s.droneCustomers)) {
//            if (r.nextDouble()<0.1) s.droneCustomers.remove(id);
//        }
    }

    private List<Solution> selectNext(List<Solution> population) {
        // 1. Fast non-dominated sort
        List<List<Solution>> fronts = fastNonDominatedSort(population);
        List<Solution> nextPop = new ArrayList<>();
        int i = 0;
        while (nextPop.size() + fronts.get(i).size() <= popSize) {
            assignCrowdingDistance(fronts.get(i));
            nextPop.addAll(fronts.get(i));
            if(fronts.size() - 1 != i) {
                i++;
            }
        }
        // fill remaining slots from front i
        List<Solution> lastFront = fronts.get(i);
        assignCrowdingDistance(lastFront);
        lastFront.sort((s1, s2) -> Double.compare(s2.crowdingDistance, s1.crowdingDistance));
        int remaining = popSize - nextPop.size();
        nextPop.addAll(lastFront.subList(0, remaining));
        return nextPop;
    }

    /**
     * Fast non-dominated sorting
     */
    private List<List<Solution>> fastNonDominatedSort(List<Solution> pop) {
        Map<Solution, List<Solution>> S = new HashMap<>();
        Map<Solution, Integer> n = new HashMap<>();
        List<List<Solution>> fronts = new ArrayList<>();
        fronts.add(new ArrayList<>());

        for (Solution p : pop) {
            S.put(p, new ArrayList<>());
            n.put(p, 0);
            for (Solution q : pop) {
                if (dominates(p, q)) {
                    S.get(p).add(q);
                } else if (dominates(q, p)) {
                    n.put(p, n.get(p) + 1);
                }
            }
            if (n.get(p) == 0) fronts.get(0).add(p);
        }
        int idx = 0;
        while (idx < fronts.size()) {
            List<Solution> next = new ArrayList<>();
            for (Solution p : fronts.get(idx)) {
                for (Solution q : S.get(p)) {
                    n.put(q, n.get(q) - 1);
                    if (n.get(q) == 0) next.add(q);
                }
            }
            if (!next.isEmpty()) fronts.add(next);
            idx++;
        }
        return fronts;
    }

    /**
     * Assign crowding distance to each solution in a front
     */
    private void assignCrowdingDistance(List<Solution> front) {
        int l = front.size();
        if (l == 0) return;
        for (Solution s : front) s.crowdingDistance = 0;
        // objectives: makespan, carbonEmission
        front.sort(Comparator.comparingDouble(s -> s.makespan));
        front.get(0).crowdingDistance = front.get(l-1).crowdingDistance = Double.POSITIVE_INFINITY;
        double minM = front.get(0).makespan, maxM = front.get(l-1).makespan;
        for (int i = 1; i < l-1; i++) {
            front.get(i).crowdingDistance += (front.get(i+1).makespan - front.get(i-1).makespan) / (maxM - minM);
        }
        front.sort(Comparator.comparingDouble(s -> s.carbonEmission));
        front.get(0).crowdingDistance = front.get(l-1).crowdingDistance = Double.POSITIVE_INFINITY;
        double minC = front.get(0).carbonEmission, maxC = front.get(l-1).carbonEmission;
        for (int i = 1; i < l-1; i++) {
            front.get(i).crowdingDistance += (front.get(i+1).carbonEmission - front.get(i-1).carbonEmission) / (maxC - minC);
        }
    }

    private List<Solution> getParetoFront(List<Solution> pop) {
        List<Solution> front = new ArrayList<>();
        for (Solution s: pop) {
            boolean dom=false;
            for (Solution o: pop) if (dominates(o,s) && !o.equals(s)) dom=true;
            if (!dom) front.add(s);
        }
        return front;
    }

    private boolean dominates(Solution a, Solution b) {
        return (a.makespan<=b.makespan&&a.carbonEmission<b.carbonEmission)||
                (a.makespan<b.makespan&&a.carbonEmission<=b.carbonEmission);
    }

    private List<Integer> flatten(List<List<Integer>> routes) {
        List<Integer> flat = new ArrayList<>();
        for (List<Integer> r: routes) for (int id: r) if (id!=0) flat.add(id);
        return flat;
    }
    private List<List<Integer>> split(List<Integer> seq, int k) {
        int chunk = (int)Math.ceil((double)seq.size()/k);
        List<List<Integer>> m = new ArrayList<>();
        for (int i=0;i<k;i++){
            int s=i*chunk,e=Math.min(s+chunk,seq.size());
            List<Integer> r=new ArrayList<>();r.add(0);
            if(s<e) r.addAll(seq.subList(s,e));
            r.add(0);
            m.add(r);
        }
        return m;
    }

    public List<Solution> optimizeDroneAssignment(List<Solution> pop, int droneGens) {
        var popAbsolute = new ArrayList<Solution>();
        for (Solution pt : pop) {
            var truckSol = pt.cloneSolution();
            truckSol.droneCustomers.clear();
            List<Solution> dronePop = initDronePop(truckSol);
            dronePop.forEach(d -> evaluateCombined(d, truckSol));

            for (int g = 0; g < droneGens; g++) {
                List<Solution> off = new ArrayList<>();
                while (off.size() < popSize) {
                    List<Solution> sol = null;
                    do {
                        Solution p1 = tournament(dronePop);
                        Solution p2 = tournament(dronePop);
                        sol = crossoverDrone(p1, p2);
                    } while (sol == null);
                    for (Solution child : sol) {
                        mutateDrone(child);
                        evaluateCombined(child, truckSol);  //x2
                        off.add(child);
                        if (off.size() >= popSize) break;
                    }

                }
                dronePop.addAll(off);
                dronePop = selectNext(dronePop);
            }

            // gán drone tốt nhất vào truckSol
            Solution best = getParetoFront(dronePop).get(0).cloneSolution(); // lấy front tốt nhất
            popAbsolute.add(best);
//            truckSol.droneCustomers = new HashSet<>(best.droneCustomers);
//            truckSol.makespan = best.makespan;
//            truckSol.carbonEmission = best.carbonEmission;
//            truckSol.serviceLevel = best.serviceLevel;
        }
        return popAbsolute;
    }

    private List<Solution> initDronePop(Solution truckSol) {
        List<Solution> pop = new ArrayList<>();

        int n = 0;
        while (pop.size() < popSize) {
            Solution s = truckSol.cloneSolution();

            for (List<Integer> route : truckSol.truckRoutes) {
                var c = new ArrayList<>(route.subList(1, route.size() - 1));
                List<Integer> selected = getRandomNonAdjacentSelection(c);
                s.droneCustomers.addAll(selected);
            }

            if (!isDuplicate(s.droneCustomers, pop)) {
                pop.add(s);
            }
            n++;
        }
//        while (pop.size() < popSize) {
//            Solution s = truckSol.cloneSolution();
//            for (List<Integer> route : s.truckRoutes) {
//                boolean prevWasDrone = false;
//                for (int idx = 1; idx < route.size() - 1; idx++) {
//                    int custId = route.get(idx);
//                    // check eligibility if Customer model provides it
//                    Customer c = customers.get(custId);
//                    boolean eligible = (c == null) || c.droneEligible;
//                    if (eligible && !prevWasDrone && rnd.nextDouble() < 0.9) {
//                        s.droneCustomers.add(custId);
//                        prevWasDrone = true;
//                    } else {
//                        prevWasDrone = false;
//                    }
//                }
//            }
//            pop.add(s);
//        }

        return pop;
    }

    private List<Solution> crossoverDrone(Solution p1, Solution p2) {
        Solution c = p1.cloneSolution();
        c.droneCustomers.clear();
        var a = p1.cloneSolution();
        var b = p2.cloneSolution();
        List<Integer> a1 = getDroneGen(a);
        List<Integer> a2 = getDroneGen(b);

        for(int i=0;i<a1.size();i++) {
            if(Math.random() < 0.3 && !a1.get(i).equals(a2.get(i))) {
                var temp = a1.get(i);
                a1.set(i, a2.get(i));
                a2.set(i, temp);
            }
        }
        var valid1 = isValidDroneGen(a, a1);
        var valid2 = isValidDroneGen(b, a2);
        if( valid1 != null && valid2 != null) {
            a.droneCustomers = valid1;
            b.droneCustomers = valid2;
            return List.of(a,b);
        }
        return null;
    }

    private void mutateDrone(Solution s) {
        Random r = new Random();
        List<Integer> all = new ArrayList<>();
        for (List<Integer> rList : s.truckRoutes) for (int id : rList) if (id != 0) all.add(id);
        for (int i = 0;i<all.size();i++) {
            if (r.nextDouble() < 0.1) {
                if (s.droneCustomers.contains(all.get(i))) s.droneCustomers.remove(all.get(i));
                else if((i > 0 && !s.droneCustomers.contains(all.get(i-1))) && (i < all.size() - 1 && !s.droneCustomers.contains(all.get(i+1))))
                    s.droneCustomers.add(all.get(i));
            }
        }
    }

    private void evaluateCombined(Solution droneSol, Solution truckSol) {
        Solution merged = truckSol.cloneSolution();
        merged.droneCustomers = new HashSet<>(droneSol.droneCustomers);
        Evaluator.evaluate(merged, customers);
        droneSol.makespan = merged.makespan;
        droneSol.carbonEmission = merged.carbonEmission;
        droneSol.serviceLevel = merged.serviceLevel;
    }

    public static <T> List<T> getRandomNonAdjacentSelection(List<T> list) {
        int n = list.size();
        if (n == 0) return Collections.emptyList();

        Random random = new Random();
        int maxAttempts = 1000;

        for (int attempt = 0; attempt < maxAttempts; attempt++) {
            int mask = random.nextInt((1 << n)); // random từ 0 đến 2^n - 1

            if (mask == 0 || hasAdjacentOnes(mask)) continue;

            List<T> combination = new ArrayList<>();
            for (int i = 0; i < n; i++) {
                if ((mask & (1 << i)) != 0) {
                    combination.add(list.get(i));
                }
            }
            return combination;
        }

        // Không tìm được tổ hợp hợp lệ sau nhiều lần thử
        return Collections.emptyList();
    }

    private static boolean hasAdjacentOnes(int mask) {
        return (mask & (mask >> 1)) != 0;
    }

    private static boolean isDuplicate(Set<Integer> droneCustomers, List<Solution> pop) {
        for (Solution existing : pop) {
            if (existing.droneCustomers.equals(droneCustomers)) {
                return true;
            }
        }
        return false;
    }

    public static List<Integer> getDroneGen(Solution solution) {
        List<Integer> a1 = new ArrayList<>();
        for(var i : solution.truckRoutes) {
            for(var j : i) {
                if(solution.droneCustomers.contains(j)) {
                    a1.add(1);
                } else {
                    a1.add(0);
                }
            }
        }
        return a1;
    }

    public Set<Integer> isValidDroneGen(Solution solution, List<Integer> list) {
        Set<Integer> rs = new HashSet<>();

        var truckRoute = solution.truckRoutes;
        int i = 0;
        for (List<Integer> integers : truckRoute) {
            for (int j = 0; j < truckRoute.get(0).size(); j++) {
                if (list.get(i) == 1) {
                    if ((i > 0 && list.get(i - 1) == 1) || (i < list.size() - 1 && list.get(i + 1) == 1)) {
                        return null;
                    } else {
                        rs.add(integers.get(j));
                        i++;
                    }
                }
            }
        }
        return rs;
    }
}