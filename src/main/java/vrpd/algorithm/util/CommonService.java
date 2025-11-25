package vrpd.algorithm.util;

import org.knowm.xchart.SwingWrapper;
import org.knowm.xchart.XYChart;
import org.knowm.xchart.XYChartBuilder;
import org.knowm.xchart.XYSeries;
import vrpd.algorithm.model.Solution;

import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;

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
}
