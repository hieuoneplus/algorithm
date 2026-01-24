package vrpd.algorithm;

import org.apache.poi.ss.usermodel.*;
import org.apache.poi.xssf.usermodel.XSSFWorkbook;
import vrpd.algorithm.model.CountGeneration;
import vrpd.algorithm.model.Solution;
import vrpd.algorithm.util.IGDUtil;

import java.io.*;
import java.nio.file.Files;
import java.nio.file.Path;
import java.nio.file.Paths;
import java.text.DecimalFormat;
import java.util.Arrays;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import static vrpd.algorithm.util.CMetric.computePairwiseCoverage;
import static vrpd.algorithm.util.CommonService.*;

public class CommonRun {
    private static DecimalFormat df = new DecimalFormat("#,###.######");
    public static void main(String[] args) throws IOException {
        String folderPath = "src/data/100";
        File folder = new File(folderPath);

        // Tên file Excel sẽ lưu kết quả
        String outExcel = "src/output/results_100_1.xlsx";

        // Tạo workbook và sheet
        try (Workbook workbook = new XSSFWorkbook()) {
            Sheet sheetHV = workbook.createSheet("HyperVolume");
            Sheet sheetC = workbook.createSheet("C_metric");
            Sheet sheetPA = workbook.createSheet("PA");
            Sheet sheetRT = workbook.createSheet("Running Time");
            Sheet sheetGenCount = workbook.createSheet("Generation Count");
            Sheet sheetIGD = workbook.createSheet("Inverted Generational Distance");
            // Tạo header
            createHeaderHV(sheetHV);
            createHeaderC(sheetC);
            createHeaderPA(sheetPA);
            createHeaderRT(sheetRT);
            createHeaderGC(sheetGenCount);
            createHeaderIGD(sheetIGD);
            // điều chỉnh pattern nếu muốn

            File[] csvFiles = folder.listFiles((dir, name) -> name.toLowerCase().endsWith(".csv"));
            if (csvFiles != null) {
                int rowIdx = 1;

                for (File file : csvFiles) {


                    long start1 = System.nanoTime();
                    var par1 = MoeadVrpd.run(Paths.get(folderPath, file.getName()).toString());
                    long end1 = System.nanoTime();
                    double timeMoeadMs = (end1 - start1) / 1_000_000.0;

                    long start2 = System.nanoTime();
                    var par2 = ModeVrpd.run(Paths.get(folderPath, file.getName()).toString());
                    long end2 = System.nanoTime();
                    double timeModeMs = (end2 - start2) / 1_000_000.0;

                    long start3 = System.nanoTime();
                    var par3 = MotlboVrpd.run(Paths.get(folderPath, file.getName()).toString());
                    long end3 = System.nanoTime();
                    double timeMotlboMs = (end3 - start3) / 1_000_000.0;

                    long start4 = System.nanoTime();
                    var par4 = NsgaIIVrpd.run(Paths.get(folderPath, file.getName()).toString());
                    long end4 = System.nanoTime();
                    double timeNsgaIIMs = (end4 - start4) / 1_000_000.0;

                    Map<String, Double> rt = new LinkedHashMap<>();
                    rt.put("MOEAD", timeMoeadMs);
                    rt.put("MODE", timeModeMs);
                    rt.put("MOTLBO", timeMotlboMs);
                    rt.put("NSGAII", timeNsgaIIMs);

                    Map<String, List<Solution>> results = new LinkedHashMap<>();
                    results.put("MOEAD", par1);
                    results.put("MODE", par2);
                    results.put("MOTLBO", par3);
                    results.put("NSGAII", par4);

                    createDataHV(copyCommon(results), sheetHV, rowIdx, file);
                    createDataC(copyCommon(results), sheetC, rowIdx, file);
                    createDataPA(copyCommon(results), sheetPA, rowIdx, file);
                    createDataRT(rt, sheetRT, rowIdx, file);
                    createDataGC(sheetGenCount, rowIdx, file);
                    createDataIGD(copyCommon(results), sheetIGD, rowIdx,file);
                    rowIdx++;
                }

                // autosize columns (tối đa 10 cột tránh chiếm quá lâu nếu nhiều cột)
                for (int c = 0; c <= 4; c++) {
                    sheetHV.autoSizeColumn(c);
                    sheetC.autoSizeColumn(c);
                    sheetPA.autoSizeColumn(c);
                }

                // Ghi workbook ra file
                try (FileOutputStream fos = new FileOutputStream(outExcel)) {
                    workbook.write(fos);
                }

                System.out.println("Đã ghi kết quả vào " + outExcel + " — tổng " + csvFiles.length + " file.");
            } else {

            }
        }
    }

    public static void createHeaderHV(Sheet sheet) {
        Row header = sheet.createRow(0);
        header.createCell(0).setCellValue("Data");
        header.createCell(1).setCellValue("HV Moead");
        header.createCell(2).setCellValue("HV Mode");
        header.createCell(3).setCellValue("HV Motlbo");
        header.createCell(4).setCellValue("HV NSGA-II");
    }
    public static void createHeaderC(Sheet sheet) {
        Row header = sheet.createRow(0);
        header.createCell(0).setCellValue("Data");
        header.createCell(1).setCellValue("C /Mode");
        header.createCell(2).setCellValue("C /Motlbo");
        header.createCell(3).setCellValue("C /NSGA-II");
    }
    public static void createHeaderPA(Sheet sheet) {
        Row header = sheet.createRow(0);
        header.createCell(0).setCellValue("Data");
        header.createCell(1).setCellValue("PA Moead");
        header.createCell(2).setCellValue("PA Mode");
        header.createCell(3).setCellValue("PA Motlbo");
        header.createCell(4).setCellValue("PA NSGA-II");
    }

    public static void createHeaderRT(Sheet sheet) {
        Row header = sheet.createRow(0);
        header.createCell(0).setCellValue("Data");
        header.createCell(1).setCellValue("Running time Moead ms");
        header.createCell(2).setCellValue("Running time Mode ms");
        header.createCell(3).setCellValue("Running time Motlbo ms");
        header.createCell(4).setCellValue("Running time NSGA-II ms");
    }

    public static void createHeaderIGD(Sheet sheet) {
        Row header = sheet.createRow(0);
        header.createCell(0).setCellValue("Data");
        header.createCell(1).setCellValue("IGD Moead");
        header.createCell(2).setCellValue("IGD Mode");
        header.createCell(3).setCellValue("IGD Motlbo");
        header.createCell(4).setCellValue("IGD NSGA-II");
    }

    public static void createHeaderGC(Sheet sheet) {
        Row header = sheet.createRow(0);
        header.createCell(0).setCellValue("Data");
        header.createCell(1).setCellValue("Generation Count Moead ms");
        header.createCell(2).setCellValue("Generation Count Mode ms");
        header.createCell(3).setCellValue("Generation Count Motlbo ms");
        header.createCell(4).setCellValue("Generation Count NSGA-II ms");
    }

    public static void createDataHV(Map<String, List<Solution>> results, Sheet sheet, int rowIdx, File file) {
        Row row = sheet.createRow(rowIdx);

        // Cột filename
        row.createCell(0).setCellValue(file.getName());

        // Cột HV: nếu HV là số (Double), ghi trực tiếp; nếu object/array, chuyển sang String
        writeCell(row, 1, calculateHypervolume(results.get("MOEAD")), df);
        writeCell(row, 2, calculateHypervolume(results.get("MODE")), df);
        writeCell(row, 3, calculateHypervolume(results.get("MOTLBO")), df);
        writeCell(row, 4, calculateHypervolume(results.get("NSGAII")), df);
    }
    public static void createDataC(Map<String, List<Solution>> results, Sheet sheet, int rowIdx, File file) {
        Row row = sheet.createRow(rowIdx);

        var map = computePairwiseCoverage(results, true);
        // Cột filename
        row.createCell(0).setCellValue(file.getName());

        // Cột C: nếu HV là số (Double), ghi trực tiếp; nếu object/array, chuyển sang String
        writeCell(row, 1, map.get("MOEAD").get("MODE") - map.get("MODE").get("MOEAD"), df);
        writeCell(row, 2, map.get("MOEAD").get("MOTLBO") - map.get("MOTLBO").get("MOEAD"), df);
        writeCell(row, 3, map.get("MOEAD").get("NSGAII") - map.get("NSGAII").get("MOEAD"), df);
    }
    public static void createDataPA(Map<String, List<Solution>> results, Sheet sheet, int rowIdx, File file) {
        Row row = sheet.createRow(rowIdx);

        row.createCell(0).setCellValue(file.getName());

        // Ghi số lượng nghiệm Pareto vào Excel
        writeCell(row, 1, results.get("MOEAD").size(), df);
        writeCell(row, 2, results.get("MODE").size(), df);
        writeCell(row, 3, results.get("MOTLBO").size(), df);
        writeCell(row, 4, results.get("NSGAII").size(), df);

        // Ghi PA chi tiết ra file txt
        writePAtoTxt(results.get("MOEAD"), "MOEAD", file);
        writePAtoTxt(results.get("MODE"), "MODE", file);
        writePAtoTxt(results.get("MOTLBO"), "MOTLBO", file);
        writePAtoTxt(results.get("NSGAII"), "NSGAII", file);
    }

    public static void createDataRT(Map<String, Double> results, Sheet sheet, int rowIdx, File file) {
        Row row = sheet.createRow(rowIdx);

        // Cột filename
        row.createCell(0).setCellValue(file.getName());

        // Cột HV: nếu HV là số (Double), ghi trực tiếp; nếu object/array, chuyển sang String
        writeCell(row, 1, results.get("MOEAD"), df);
        writeCell(row, 2, results.get("MODE"), df);
        writeCell(row, 3, results.get("MOTLBO"), df);
        writeCell(row, 4, results.get("NSGAII"), df);
    }

    public static void createDataGC(Sheet sheet, int rowIdx, File file) {
        Row row = sheet.createRow(rowIdx);

        // Cột filename
        row.createCell(0).setCellValue(file.getName());

        // Cột HV: nếu HV là số (Double), ghi trực tiếp; nếu object/array, chuyển sang String
        writeCell(row, 1, CountGeneration.MOEAD, df);
        writeCell(row, 2, CountGeneration.MODE, df);
        writeCell(row, 3, CountGeneration.MOTLBO, df);
        writeCell(row, 4, CountGeneration.NSGA_II, df);
    }

    public static void createDataIGD(Map<String, List<Solution>> results, Sheet sheet, int rowIdx, File file) {
        Row row = sheet.createRow(rowIdx);

        // Cột filename
        row.createCell(0).setCellValue(file.getName());

        var igdMap = IGDUtil.computeIGDs(results, true);
        writeCell(row, 1, igdMap.get("MOEAD"), df);
        writeCell(row, 2, igdMap.get("MODE"), df);
        writeCell(row, 3, igdMap.get("MOTLBO"), df);
        writeCell(row, 4, igdMap.get("NSGAII"), df);
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
        if (value.toString().length() > 32767) {
            System.out.println("Too long: " + value.toString());
        } else {
            cell.setCellValue(value.toString());
        }
    }

    public static void writePAtoTxt(
            List<Solution> solutions,
            String algorithm,
            File dataFile
    ) {
        try {
            Path outDir = Path.of("src/output/PA");
            Files.createDirectories(outDir);

            String fileName = dataFile.getName().replace(".csv", "")
                    + "_" + algorithm + ".txt";

            Path outFile = outDir.resolve(fileName);

            try (BufferedWriter bw = new BufferedWriter(new FileWriter(outFile.toFile()))) {
                for (Solution s : solutions) {
                    bw.write(s.toString()); // hoặc custom format
                    bw.newLine();
                }
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

}
