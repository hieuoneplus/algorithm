package vrpd.algorithm.util;

import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import vrpd.algorithm.model.Customer;
import vrpd.algorithm.model.Solution;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.*;

public class CommonService {
    public static String calculateHypervolume(List<Solution> pareto) {
        if (pareto == null || pareto.isEmpty()) return "0.0";

        // Điểm tham chiếu cố định
        final double REF1 = 2.0;
        final double REF2 = 2.0;

        // (Tùy chọn) lọc điểm nằm ngoài hộp tham chiếu
        List<Solution> pts = new ArrayList<>();
        for (Solution s : pareto) {
            if (s != null && s.makespan <= REF1 && s.carbonEmission <= REF2) {
                pts.add(s);
            }
        }
        if (pts.isEmpty()) return "0.0";

        // Sắp xếp theo f1 (makespan) tăng dần
        pts.sort(Comparator.comparingDouble(s -> s.makespan));

        // Quét ngược để tránh chồng lấp (HV 2D, minimization)
        double hv = 0.0;
        double prevX = REF1; // biên phải ban đầu là ref.x

        for (int i = pts.size() - 1; i >= 0; i--) {
            Solution s = pts.get(i);

            // Cắt trong hộp tham chiếu (an toàn)
            double x = Math.min(s.makespan, REF1);
            double y = Math.min(s.carbonEmission, REF2);

            double width  = Math.max(0.0, prevX - x);
            double height = Math.max(0.0, REF2 - y);

            hv += width * height;
            prevX = x; // cập nhật biên phải cho ô tiếp theo
        }

        // Format kết quả
        DecimalFormat df = new DecimalFormat("#,###.####");
        return df.format(Math.round(hv * 10000.0) / 10000.0);
    }


    public static void drawImg(List<Solution> pareto, String algor) {
        XYChart chart = new XYChartBuilder().width(800).height(600).title("Pareto " + algor).xAxisTitle("CO2").yAxisTitle("Makespan").build();
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

    public static Map<Integer, Customer> loadCustomersV2(String path) throws IOException {
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
                boolean droneServe = Integer.parseInt(t[6]) == 1;
                map.put(id, new Customer(id,x,y,d,tws,twe,droneServe));
                id++;
            }
        }
        return map;
    }

    public static Map<String, List<Solution>> copyCommon(Map<String, List<Solution>> originalMap) {
        Map<String, List<Solution>> deepCopy = new LinkedHashMap<>();

        // Duyệt qua từng entry của Map gốc
        for (Map.Entry<String, List<Solution>> entry : originalMap.entrySet()) {
            List<Solution> copiedList = new ArrayList<>();

            // Duyệt qua từng phần tử trong List của mỗi entry và sao chép nó
            for (Solution item : entry.getValue()) {
                copiedList.add(item.cloneSolution()); // Gọi phương thức cloneSolution() để sao chép item
            }
            deepCopy.put(entry.getKey(), copiedList); // Thêm vào bản sao
        }

        return deepCopy;
    }

    public static String printPA(List<Solution> solutions) {
        final String[] ss = {""};
        solutions.forEach(s -> {
            ss[0] = ss[0] + s.toString() + "\n";
        });
        return ss[0];
    }

    public static List<Solution> getParetoFront(List<Solution> pop) {
        List<Solution> front = new ArrayList<>();
        for (Solution s: pop) {
            boolean dom=false;
            for (Solution o: pop) if (dominates(o,s) && !o.equals(s)) dom=true;
            if (!dom) front.add(s);
        }
        return front;
    }
    private static boolean dominates(Solution a, Solution b) {
        return (a.makespan<=b.makespan&&a.carbonEmission<b.carbonEmission)||
                (a.makespan<b.makespan&&a.carbonEmission<=b.carbonEmission);
    }
}
