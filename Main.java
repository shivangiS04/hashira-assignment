import java.io.*;
import java.util.*;
import org.json.*;

public class Main {
    
    public static long baseToDecimal(String value, int base) {
        long result = 0;
        long power = 1;
        
        for (int i = value.length() - 1; i >= 0; i--) {
            char digit = Character.toLowerCase(value.charAt(i));
            int digitValue;
            
            if (digit >= '0' && digit <= '9') {
                digitValue = digit - '0';
            } else if (digit >= 'a' && digit <= 'z') {
                digitValue = digit - 'a' + 10;
            } else {
                throw new IllegalArgumentException("Invalid digit: " + digit);
            }
            
            if (digitValue >= base) {
                throw new IllegalArgumentException("Invalid digit for base " + base + ": " + digit);
            }
            
            result += digitValue * power;
            power *= base;
        }
        
        return result;
    }
    
    // Method 1: Lagrange Interpolation
    public static long lagrangeInterpolation(List<long[]> points, int k) {
        double constantTerm = 0.0;
        
        for (int i = 0; i < k; i++) {
            long xi = points.get(i)[0];
            long yi = points.get(i)[1];
            
            double li0 = 1.0;
            for (int j = 0; j < k; j++) {
                if (i != j) {
                    long xj = points.get(j)[0];
                    li0 *= (0.0 - xj) / (xi - xj);
                }
            }
            
            constantTerm += yi * li0;
        }
        
        return Math.round(constantTerm);
    }
    
    // Method 2: Newton's Divided Differences
    public static long newtonDividedDifferences(List<long[]> points, int k) {
        double[][] dividedDiff = new double[k][k];
        
        for (int i = 0; i < k; i++) {
            dividedDiff[i][0] = points.get(i)[1];
        }
        
        for (int j = 1; j < k; j++) {
            for (int i = 0; i < k - j; i++) {
                long xi = points.get(i)[0];
                long xi_plus_j = points.get(i + j)[0];
                dividedDiff[i][j] = (dividedDiff[i + 1][j - 1] - dividedDiff[i][j - 1]) / (xi_plus_j - xi);
            }
        }
        
        double result = dividedDiff[0][0];
        double product = 1;
        
        for (int i = 1; i < k; i++) {
            product *= (0 - points.get(i - 1)[0]);
            result += dividedDiff[0][i] * product;
        }
        
        return Math.round(result);
    }
    
    // Method 3: Cramer's Rule
    public static long cramersRule(List<long[]> points, int k) {
        double[][] matrix = new double[k][k];
        double[] constants = new double[k];
        
        for (int i = 0; i < k; i++) {
            long x = points.get(i)[0];
            long y = points.get(i)[1];
            
            for (int j = 0; j < k; j++) {
                matrix[i][j] = Math.pow(x, j);
            }
            constants[i] = y;
        }
        
        double det = determinant(matrix, k);
        if (Math.abs(det) < 1e-10) {
            throw new RuntimeException("Matrix is singular");
        }
        
        double[][] modMatrix = new double[k][k];
        for (int i = 0; i < k; i++) {
            for (int j = 0; j < k; j++) {
                if (j == 0) {
                    modMatrix[i][j] = constants[i];
                } else {
                    modMatrix[i][j] = matrix[i][j];
                }
            }
        }
        
        double detMod = determinant(modMatrix, k);
        return Math.round(detMod / det);
    }
    
    // Method 4: Matrix Inversion (LU Decomposition)
    public static long luDecomposition(List<long[]> points, int k) {
        double[][] matrix = new double[k][k];
        double[] constants = new double[k];
        
        for (int i = 0; i < k; i++) {
            long x = points.get(i)[0];
            long y = points.get(i)[1];
            
            for (int j = 0; j < k; j++) {
                matrix[i][j] = Math.pow(x, j);
            }
            constants[i] = y;
        }
        
        double[][] L = new double[k][k];
        double[][] U = new double[k][k];
        
        for (int i = 0; i < k; i++) {
            for (int j = i; j < k; j++) {
                double sum = 0;
                for (int p = 0; p < i; p++) {
                    sum += L[i][p] * U[p][j];
                }
                U[i][j] = matrix[i][j] - sum;
            }
            
            for (int j = i; j < k; j++) {
                if (i == j) {
                    L[i][i] = 1;
                } else {
                    double sum = 0;
                    for (int p = 0; p < i; p++) {
                        sum += L[j][p] * U[p][i];
                    }
                    L[j][i] = (matrix[j][i] - sum) / U[i][i];
                }
            }
        }
        
        double[] y = new double[k];
        for (int i = 0; i < k; i++) {
            double sum = 0;
            for (int j = 0; j < i; j++) {
                sum += L[i][j] * y[j];
            }
            y[i] = constants[i] - sum;
        }
        
        double[] x = new double[k];
        for (int i = k - 1; i >= 0; i--) {
            double sum = 0;
            for (int j = i + 1; j < k; j++) {
                sum += U[i][j] * x[j];
            }
            x[i] = (y[i] - sum) / U[i][i];
        }
        
        return Math.round(x[0]);
    }
    
    // Method 5: Barycentric Interpolation
    public static long barycentricInterpolation(List<long[]> points, int k) {
        double[] weights = new double[k];
        
        for (int i = 0; i < k; i++) {
            weights[i] = 1.0;
            for (int j = 0; j < k; j++) {
                if (i != j) {
                    weights[i] /= (points.get(i)[0] - points.get(j)[0]);
                }
            }
        }
        
        double numerator = 0;
        double denominator = 0;
        
        for (int i = 0; i < k; i++) {
            long xi = points.get(i)[0];
            long yi = points.get(i)[1];
            
            if (xi == 0) {
                return yi;
            }
            
            double term = weights[i] / (0 - xi);
            numerator += term * yi;
            denominator += term;
        }
        
        return Math.round(numerator / denominator);
    }
    
    // Helper method for determinant calculation
    private static double determinant(double[][] matrix, int n) {
        if (n == 1) {
            return matrix[0][0];
        }
        
        double det = 0;
        double[][] temp = new double[n - 1][n - 1];
        
        for (int i = 0; i < n; i++) {
            int row = 0;
            for (int j = 1; j < n; j++) {
                int col = 0;
                for (int k = 0; k < n; k++) {
                    if (k != i) {
                        temp[row][col++] = matrix[j][k];
                    }
                }
                row++;
            }
            det += Math.pow(-1, i) * matrix[0][i] * determinant(temp, n - 1);
        }
        
        return det;
    }
    
    public static void main(String[] args) throws Exception {
        BufferedReader br = new BufferedReader(new InputStreamReader(System.in));
        StringBuilder sb = new StringBuilder();
        String line;
        
        while ((line = br.readLine()) != null) {
            sb.append(line);
        }
        
        JSONObject obj = new JSONObject(sb.toString());
        JSONObject keys = obj.getJSONObject("keys");
        int k = keys.getInt("k");
        
        System.out.println("=== Shamir's Secret Sharing - Multiple Methods ===");
        System.out.println("k (minimum points needed): " + k);
        
        List<long[]> points = new ArrayList<>();
        
        for (String key : obj.keySet()) {
            if (key.equals("keys")) continue;
            
            JSONObject pointObj = obj.getJSONObject(key);
            int base = Integer.parseInt(pointObj.getString("base"));
            String value = pointObj.getString("value");
            
            long x = Long.parseLong(key.replaceAll("\\D", ""));
            long y = baseToDecimal(value, base);
            
            points.add(new long[]{x, y});
            System.out.println("Point: (" + x + ", " + y + ") [Base " + base + ": " + value + "]");
        }
        
        points.sort((a, b) -> Long.compare(a[0], b[0]));
        
        if (points.size() < k) {
            System.err.println("Error: Not enough points. Need " + k + ", got " + points.size());
            return;
        }
        
        System.out.println("\n" + "=".repeat(50));
        System.out.println("TESTING ALL METHODS:");
        System.out.println("=".repeat(50));
        
        try {
            long result1 = lagrangeInterpolation(points, k);
            System.out.println("Method 1 - Lagrange Interpolation: " + result1);
        } catch (Exception e) {
            System.out.println("Method 1 - Lagrange Interpolation: ERROR - " + e.getMessage());
        }
        
        try {
            long result2 = newtonDividedDifferences(points, k);
            System.out.println("Method 2 - Newton's Divided Differences: " + result2);
        } catch (Exception e) {
            System.out.println("Method 2 - Newton's Divided Differences: ERROR - " + e.getMessage());
        }
        
        try {
            long result3 = cramersRule(points, k);
            System.out.println("Method 3 - Cramer's Rule: " + result3);
        } catch (Exception e) {
            System.out.println("Method 3 - Cramer's Rule: ERROR - " + e.getMessage());
        }
        
        try {
            long result4 = luDecomposition(points, k);
            System.out.println("Method 4 - LU Decomposition: " + result4);
        } catch (Exception e) {
            System.out.println("Method 4 - LU Decomposition: ERROR - " + e.getMessage());
        }
        
        try {
            long result5 = barycentricInterpolation(points, k);
            System.out.println("Method 5 - Barycentric Interpolation: " + result5);
        } catch (Exception e) {
            System.out.println("Method 5 - Barycentric Interpolation: ERROR - " + e.getMessage());
        }
        
        long finalAnswer = lagrangeInterpolation(points, k);
        System.out.println("\n" + "=".repeat(30));
        System.out.println("FINAL ANSWER: " + finalAnswer);
        System.out.println("=".repeat(30));
    }
}