package vrpd.algorithm.model;

public class Customer {
    public int id;
    public double x, y;
    public double demand;
    public double twStart, twEnd;
    public double earliness;  // EET_i
    public double lateness;   // ELT_i
    public boolean droneServe;
    public double service;    // service time s_i (minutes)

    public Customer(int id, double x, double y, double demand, double twStart, double twEnd) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.demand = demand;
        this.twStart = twStart;
        this.twEnd = twEnd;
    }

    public Customer(int id, double x, double y, double demand, double twStart, double twEnd, boolean droneServe) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.demand = demand;
        this.twStart = twStart;
        this.twEnd = twEnd;
        this.droneServe = droneServe;
    }
    public Customer(int id, double x, double y,
                    int demand, double service,
                    double twStart, double twEnd,
                    double earliness, double lateness,
                    boolean droneServe) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.demand = demand;
        this.service = service;
        this.twStart = twStart;
        this.twEnd = twEnd;
        this.earliness = earliness;
        this.lateness = lateness;
        this.droneServe = droneServe;
    }

    public double distanceTo(Customer other) {
        return Math.abs(this.x - other.x) + Math.abs(this.y - other.y);
    }
    @Override
    public String toString() {
        return "C[" + id + "](x=" + x + ",y=" + y + ")";
    }
}
