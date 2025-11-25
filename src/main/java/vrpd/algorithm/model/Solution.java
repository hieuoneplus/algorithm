package vrpd.algorithm.model;

import lombok.Data;

import java.util.*;
@Data
public class Solution {
    public List<List<Integer>> truckRoutes = new ArrayList<>();
    public Set<Integer> droneCustomers = new HashSet<>();
    public double makespan;
    public double carbonEmission;
    public double serviceLevel;
    public boolean feasible;

    public double crowdingDistance;
    public Solution() {
        this.truckRoutes = new ArrayList<>();
        this.droneCustomers = new HashSet<>();
        this.makespan = Double.POSITIVE_INFINITY;
        this.carbonEmission = Double.POSITIVE_INFINITY;
        this.serviceLevel = 0.0;
        this.feasible = false;
    }

//    public Solution cloneSolution() {
//        Solution s = new Solution();
//        for (List<Integer> r : this.truckRoutes) s.truckRoutes.add(new ArrayList<>(r));
//        s.droneCustomers = new HashSet<>(this.droneCustomers);
//        return s;
//    }


    public Solution cloneSolution() {
        Solution c = new Solution();
        // deep copy routes
        c.truckRoutes = new ArrayList<>(this.truckRoutes.size());
        for (List<Integer> r : this.truckRoutes) {
            c.truckRoutes.add(new ArrayList<>(r));
        }
        c.droneCustomers = new HashSet<>(this.droneCustomers);
        c.makespan = this.makespan;
        c.carbonEmission = this.carbonEmission;
        c.serviceLevel = this.serviceLevel;
        c.feasible = this.feasible;
        return c;
    }

    @Override
    public String toString() {
        StringBuilder sb = new StringBuilder();
        sb.append("Solution{makespan=").append(makespan)
                .append(", carbon=").append(carbonEmission)
                .append(", S=").append(serviceLevel)
                .append(", feasible=").append(feasible)
                .append(", routes=");
        for (List<Integer> r : truckRoutes) {
            sb.append(r.toString()).append(";");
        }
        sb.append(", drone=").append(droneCustomers);
        sb.append("}");
        return sb.toString();
    }

    @Override
    public boolean equals(Object o) {
        if (this == o) return true;
        if (!(o instanceof Solution)) return false;
        Solution that = (Solution) o;

        // So sánh truckRoutes (deep compare)
        if (this.truckRoutes.size() != that.truckRoutes.size()) return false;
        for (int i = 0; i < this.truckRoutes.size(); i++) {
            List<Integer> route1 = this.truckRoutes.get(i);
            List<Integer> route2 = that.truckRoutes.get(i);
            if (!Objects.equals(route1, route2)) {
                return false;
            }
        }

        // So sánh droneCustomers (Set có sẵn equals deep)
        return Objects.equals(this.droneCustomers, that.droneCustomers);
    }

    @Override
    public int hashCode() {
        return Objects.hash(truckRoutes, droneCustomers);
    }
}
