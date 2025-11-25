package vrpd.algorithm;

import org.apache.poi.ss.usermodel.*;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import vrpd.algorithm.model.Solution;

import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static vrpd.algorithm.util.CommonService.*;

public class CommonRun {
    public static void main(String[] args) throws IOException {
        String folderPath = "src/data";
        File folder = new File(folderPath);

        // Tên file Excel sẽ lưu kết quả
        String outExcel = "src/output/results.xlsx";

        // Tạo workbook và sheet
        try (Workbook workbook = new XSSFWorkbook()) {
            Sheet sheet = workbook.createSheet("Results");

            // Tạo header
            Row header = sheet.createRow(0);
            header.createCell(0).setCellValue("Data");
            header.createCell(1).setCellValue("HV Moead");
            header.createCell(2).setCellValue("HV Mode");
            header.createCell(3).setCellValue("HV Motlbo");
            header.createCell(4).setCellValue("HV NSGA-II");

            // Tùy chọn format số (nếu cần)
            DecimalFormat df = new DecimalFormat("#,###.######"); // điều chỉnh pattern nếu muốn

            File[] csvFiles = folder.listFiles((dir, name) -> name.toLowerCase().endsWith(".csv"));

            if (csvFiles != null) {
                int rowIdx = 1;
                for (File file : csvFiles) {
                    var par1 = MoeadVrpd.run(Paths.get(folderPath, file.getName()).toString());
                    var par2 = ModeVrpd.run(Paths.get(folderPath, file.getName()).toString());
                    var par3 = MotlboVrpd.run(Paths.get(folderPath, file.getName()).toString());
                    var par4 = NsgaIIVrpd.run(Paths.get(folderPath, file.getName()).toString());


                    Map<String, List<Solution>> results = new LinkedHashMap<>();
                    results.put("MOEAD", par1);
                    results.put("MODE", par2);
                    results.put("MOTLBO", par3);
                    results.put("NSGAII", par4);

                    Row row = sheet.createRow(rowIdx++);

                    // Cột filename
                    row.createCell(0).setCellValue(file.getName());

                    // Cột HV: nếu HV là số (Double), ghi trực tiếp; nếu object/array, chuyển sang String
                    writeCell(row, 1, calculateHypervolume(par1), df);
                    writeCell(row, 2, calculateHypervolume(par2), df);
                    writeCell(row, 3, calculateHypervolume(par3), df);
                    writeCell(row, 4, calculateHypervolume(par4), df);
                }

                // autosize columns (tối đa 10 cột tránh chiếm quá lâu nếu nhiều cột)
                for (int c = 0; c <= 4; c++) {
                    sheet.autoSizeColumn(c);
                }

                // Ghi workbook ra file
                try (FileOutputStream fos = new FileOutputStream(outExcel)) {
                    workbook.write(fos);
                }

                System.out.println("Đã ghi kết quả vào " + outExcel + " — tổng " + csvFiles.length + " file.");
            } else {
                System.out.println("Không tìm thấy thư mục hoặc không có file nào.");
            }
        }
    }

    /**
     * Hàm helper để ghi cell: nếu là Number -> ghi là numeric cell,
     * nếu là mảng hoặc object -> gọi toString().
     */
    private static void writeCell(Row row, int columnIndex, Object value, DecimalFormat df) {
        Cell cell = row.createCell(columnIndex);

        if (value == null) {
            cell.setCellValue("");
            return;
        }

        // Nếu trả về kiểu Number (Double, Integer, Long, Float, BigDecimal ...)
        if (value instanceof Number) {
            // ghi numeric cell
            double d = ((Number) value).doubleValue();
            cell.setCellValue(d);

            // tùy chọn: set format hiển thị số (thousand separator)
            CellStyle style = row.getSheet().getWorkbook().createCellStyle();
            DataFormat fmt = row.getSheet().getWorkbook().createDataFormat();
            style.setDataFormat(fmt.getFormat("#,##0.######"));
            cell.setCellStyle(style);
            return;
        }

        // Nếu là mảng (ví dụ double[] hay Double[])
        if (value.getClass().isArray()) {
            // convert to string như "[a, b, c]" hoặc join
            StringBuilder sb = new StringBuilder();
            int len = java.lang.reflect.Array.getLength(value);
            for (int i = 0; i < len; i++) {
                Object el = java.lang.reflect.Array.get(value, i);
                if (el instanceof Number) sb.append(df.format(((Number)el).doubleValue()));
                else sb.append(el == null ? "" : el.toString());
                if (i < len - 1) sb.append(", ");
            }
            cell.setCellValue(sb.toString());
            return;
        }

        // Mặc định: ghi toString()
        cell.setCellValue(value.toString());
    }

}
