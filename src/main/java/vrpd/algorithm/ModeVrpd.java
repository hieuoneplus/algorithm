package vrpd.algorithm;

import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import vrpd.algorithm.mode.MODESolver;
import vrpd.algorithm.model.Customer;
import vrpd.algorithm.model.Evaluator;
import vrpd.algorithm.model.Solution;
import vrpd.algorithm.nsga2.NSGA2Solver;
import vrpd.algorithm.util.CommonService;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

import static vrpd.algorithm.util.CommonService.loadCustomersV2;

@SpringBootApplication
public class ModeVrpd {
    public static void main(String[] args) throws Exception {
//        Map<Integer, Customer> customers = loadCustomers("src/data/customers.csv");

        Map<Integer, Customer> customers = loadCustomersV2("src/data/h100c101.csv");
        MODESolver solver = new MODESolver(50, 100000000, 4, customers, 42L);
        List<Solution> pareto = solver.run();
        pareto.forEach(s -> System.out.println(
                "Makespan=" + s.makespan + ", CO2=" + s.carbonEmission + ", truck route: " + s.truckRoutes.toString()
                        + "\n" + ",drone route:" + NSGA2Solver.getDroneGen(s).toString()
        ));
        System.out.println("HV: " + CommonService.calculateHypervolume(pareto));
        CommonService.drawImg(pareto, "MODE");
    }
    public static List<Solution> run(String cusPath) throws IOException {
        if(Evaluator.TOTAL_EVAL == 0) {
            Evaluator.TOTAL_EVAL = 50000;
        }
        Map<Integer, Customer> customers = loadCustomersV2(cusPath);
        MODESolver solver = new MODESolver(50, 100000000, 4, customers, 42L);
        return solver.run();
    }

    private static Map<Integer, Customer> loadCustomers(String path) throws IOException {
        Map<Integer, Customer> map = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            while ((line = br.readLine())!=null) {
                if (line.startsWith("id")) continue;
                String[] t = line.split(",");
                int id = Integer.parseInt(t[0]);
                double x = Double.parseDouble(t[1]), y = Double.parseDouble(t[2]);
                int d = Integer.parseInt(t[3]);
                double tws = Double.parseDouble(t[4]), twe = Double.parseDouble(t[5]);
                map.put(id, new Customer(id,x,y,d,tws,twe));
            }
        }
        return map;
    }

}
