package vrpd.algorithm;

import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import org.springframework.boot.autoconfigure.SpringBootApplication;
import vrpd.algorithm.model.Customer;
import vrpd.algorithm.model.Solution;
import vrpd.algorithm.moead.MOEADSolver;
import vrpd.algorithm.nsga2.NSGA2Solver;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.*;

@SpringBootApplication
public class MoeadVrpd {

    public static void main(String[] args) throws Exception {
//        Map<Integer, Customer> customers = loadCustomers("src/data/customers.csv");

        Map<Integer, Customer> customers = loadCustomersV2("src/data/cus_1.csv");
        MOEADSolver solver = new MOEADSolver(50, 100000000, 4, customers, 42L);
        List<Solution> pareto = solver.run();
        pareto.forEach(s -> System.out.println(
                "Makespan=" + s.makespan + ", CO2=" + s.carbonEmission + ", truck route: " + s.truckRoutes.toString()
                        + "\n" + ",drone route:" + NSGA2Solver.getDroneGen(s).toString()
        ));
        System.out.println("HV: " + calculateHypervolume(pareto));
        drawImg(pareto);
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
    private static Map<Integer, Customer> loadCustomersV2(String path) throws IOException {
        Map<Integer, Customer> map = new HashMap<>();
        try (BufferedReader br = new BufferedReader(new FileReader(path))) {
            String line;
            int id = 0;
            while ((line = br.readLine())!=null) {
                if (line.startsWith("x")) continue;
                String[] t = line.split(",");

                double x = Double.parseDouble(t[0]), y = Double.parseDouble(t[1]);
                double d = Double.parseDouble(t[2]);
                double tws = Double.parseDouble(t[3]), twe = Double.parseDouble(t[4]);
                map.put(id, new Customer(id,x,y,d,tws,twe));
                id++;
            }
        }
        return map;
    }
    public static double calculateHypervolume(List<Solution> pareto) {
        if (pareto == null || pareto.isEmpty()) return 0.0;

        // 1. tìm điểm tham chiếu r = (max f1, max f2)
        double maxMakespan = Double.NEGATIVE_INFINITY;
        double maxCarbon = Double.NEGATIVE_INFINITY;
        for (Solution s : pareto) {
            maxMakespan = Math.max(maxMakespan, s.makespan);
            maxCarbon = Math.max(maxCarbon, s.carbonEmission);
        }
        double ref1 = 0.0;
        double ref2 = 0.0;
        double[] referencePoint = {ref1, ref2};
        // 2. sắp xếp tập nghiệm theo f1 (makespan) tăng dần
        List<Solution> sorted = new ArrayList<>(pareto);
        sorted.sort(Comparator.comparingDouble(s -> s.makespan));

        double hypervolume = 0.0;
        double prevFb = referencePoint[1];
        for (int i=0;i<sorted.size();i++) {
            var ind = sorted.get(i);
            if(i != sorted.size()-1) {
                var indLbPre = sorted.get(i+1);
                double leng = Math.abs(indLbPre.makespan - ind.makespan);
                double width = Math.abs(prevFb - ind.carbonEmission);
                hypervolume += leng*width;
            } else {
                double leng = Math.abs(prevFb - ind.makespan);
                double width = Math.abs(prevFb - ind.carbonEmission);
                hypervolume += leng*width;
            }
        }
        return Math.round(hypervolume*100.0)/100.0;
    }

    public static void drawImg(List<Solution> pareto) {
        XYChart chart = new XYChartBuilder().width(800).height(600).title("Pareto MOEAD").xAxisTitle("CO2").yAxisTitle("Makespan").build();
        chart.getStyler().setDefaultSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Scatter);

        int rank = 0;

            List<Double> Lb = new ArrayList<>();
            List<Double> ratio = new ArrayList<>();
            draw(pareto, Lb, ratio);

            // Thêm dữ liệu vào biểu đồ

            // Nếu không phải lần lặp đầu tiên, nối điểm hiện tại với điểm trước đó
            chart.addSeries("Rank " + rank, ratio, Lb);//.setXYSeriesRenderStyle(XYSeries.XYSeriesRenderStyle.Line);


        new SwingWrapper<>(chart).displayChart();

    }
    public static void draw(List<Solution> list, List<Double> Lb, List<Double> ratio) {
        for(int i=0; i<list.size();i++) {

                Lb.add(Math.round(list.get(i).makespan * 1000.0) / 1000.0);
                ratio.add(Math.round(list.get(i).carbonEmission * 1000.0) / 1000.0);

        }
    }
}
