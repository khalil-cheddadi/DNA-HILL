import java.util.HashMap;
import java.util.Random;
import java.util.Scanner;
import org.apache.commons.math4.legacy.linear.Array2DRowRealMatrix;
import org.apache.commons.math4.legacy.linear.MatrixUtils;
import org.apache.commons.math4.legacy.linear.RealMatrix;
import org.apache.commons.math4.legacy.linear.LUDecomposition;

public class App {
    public static void main(String[] args) {
        Scanner input = new Scanner(System.in);
        System.out.print("enter message:");
        String userInput = input.nextLine();
        System.out.print("enter key:");
        String key = input.nextLine();
        System.out.println("\nInput:\n" + userInput);

        String preprocessed = prepare(userInput);
        System.out.println("Preprocessing:\n" + preprocessed);

        String binaryForm = convertToBinary(preprocessed);
        System.out.println("Binary Form:\n" + binaryForm);

        String dnaForm = convertToDNA(binaryForm);
        System.out.println("DNA Form:\n" + dnaForm);

        String[] result = processDNASequence(dnaForm);
        System.out.println("Output:\n" + result[0]);
        System.out.println("Ambiguity:\n" + result[1]);

        RealMatrix keyMatrix = prepareKey(key);
        System.out.println("Key Matrix :\n");
        printMatrix(keyMatrix);

        int dim = keyMatrix.getColumnDimension();
        RealMatrix messageMatrix = prepareMassage(result[0], dim);
        System.out.println("Message Matrix :\n");
        printMatrix(messageMatrix);

        RealMatrix cipherMatrix = hillEncryption(messageMatrix, keyMatrix);
        System.out.println("Cipher Matrix :\n");
        printMatrix(cipherMatrix);

        RealMatrix inverseMatrix = gererateDecryptionMatrix(keyMatrix);
        System.out.println("Inverse Matrix :\n");
        printMatrix(inverseMatrix);

        RealMatrix plaintext = hillDecreption(cipherMatrix, inverseMatrix);
        System.out.println("plaintext Matrix :\n");
        printMatrix(plaintext);

        System.out.println(" final :\n");
        printMatrix(plaintext.subtract(messageMatrix));

        input.close();
    }

    public static String prepare(String input) {
        String result = input.replaceAll(" ", "");
        return result;
    }

    private static String convertToBinary(String input) {
        StringBuilder binaryForm = new StringBuilder();
        for (char character : input.toCharArray()) {
            String buffer = Integer.toBinaryString(character);
            binaryForm.append("0".repeat(8 - buffer.length())).append(buffer);
        }
        return binaryForm.toString();
    }

    private static String convertToDNA(String binaryForm) {
        StringBuilder dnaForm = new StringBuilder();
        String[] table = { "A", "C", "G", "U" };
        Random random = new Random();
        for (int i = 0; i < binaryForm.length(); i += 2) {
            String buffer = binaryForm.substring(i, i + 2);
            dnaForm.append(table[Integer.parseInt(buffer, 2)]);
        }
        while (dnaForm.length() % 3 != 0) {
            int rand = random.nextInt(4);
            dnaForm.insert(0, table[rand]);
        }
        return dnaForm.toString();
    }

    private static String[] processDNASequence(String dnaForm) {
        String[] result = new String[2];
        StringBuilder output = new StringBuilder();
        StringBuilder ambiguity = new StringBuilder();
        HashMap<String, String> codonsTable = initCodonsTable();
        for (int i = 0; i < dnaForm.length(); i += 3) {
            String buffer = dnaForm.substring(i, i + 3);
            String[] parts = codonsTable.get(buffer).split("-");
            output.append(parts[0]);
            ambiguity.append(parts[1]);
        }
        result[0] = output.toString();
        result[1] = ambiguity.toString();
        return result;
    }

    private static RealMatrix prepareKey(String key) {
        String processedKey = key.replaceAll("[^a-zA-Z]", "").toUpperCase();

        int n = (int) Math.sqrt(processedKey.length());

        if (!isPerfectSquare(processedKey.length())) {
            throw new Error("key is unusable");
        }

        double[][] matrixData = new double[n][n];

        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                matrixData[i][j] = processedKey.charAt(i * n + j) - 65;
                if (matrixData[i][j] >= 9) {
                    matrixData[i][j]--;
                }
            }
        }

        RealMatrix matrix = new Array2DRowRealMatrix(matrixData);
        LUDecomposition luDecomposition = new LUDecomposition(matrix);
        double determinant = luDecomposition.getDeterminant();

        if (Math.round(determinant) == 0) {
            throw new Error("key is unusable bc the resualting key matrix is irreversible");

        }

        if (!areRelativelyPrime((int) Math.abs(determinant), 25)) {
            throw new Error("key is unusable bc determinant is not prime with 25");
        }

        return matrix;
    }

    private static RealMatrix prepareMassage(String message, int cols) {

        int toFill = message.length() % cols;
        if (toFill > 0) {
            message += "X".repeat(cols - toFill);
        }

        int rows = message.length() / cols;
        double[][] matrixData = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                matrixData[i][j] = message.charAt(i * cols + j) - 65;
                if (matrixData[i][j] >= 9) {
                    matrixData[i][j]--;
                }
            }
        }

        RealMatrix matrix = new Array2DRowRealMatrix(matrixData);
        return matrix;
    }

    private static RealMatrix hillEncryption(RealMatrix messageMatrix, RealMatrix keyMatrix) {
        RealMatrix result = messageMatrix.multiply(keyMatrix);
        int rows = result.getRowDimension();
        int cols = result.getColumnDimension();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                while (result.getEntry(i, j) >= 25) {
                    result.setEntry(i, j, result.getEntry(i, j) - 25);
                }
                while (result.getEntry(i, j) < 0) {
                    result.setEntry(i, j, result.getEntry(i, j) + 25);
                }
            }
        }
        return result;
    }

    private static RealMatrix gererateDecryptionMatrix(RealMatrix keyMatrix) {
        int rows = keyMatrix.getRowDimension();
        int cols = keyMatrix.getColumnDimension();
        RealMatrix inverseMatrix = MatrixUtils.inverse(keyMatrix);
        LUDecomposition luDecomposition = new LUDecomposition(keyMatrix);
        double determinant = Math.abs(luDecomposition.getDeterminant());
        int num = 0;
        while (determinant * num % 25 != 1) {
            num++;
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                inverseMatrix.setEntry(i, j, Math.round(inverseMatrix.getEntry(i, j) * determinant * num) % 25);
            }
        }

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                while (inverseMatrix.getEntry(i, j) >= 25) {
                    inverseMatrix.setEntry(i, j, inverseMatrix.getEntry(i, j) - 25);
                }
                while (inverseMatrix.getEntry(i, j) < 0) {
                    inverseMatrix.setEntry(i, j, inverseMatrix.getEntry(i, j) + 25);
                }
            }
        }

        return inverseMatrix;
    }

    static RealMatrix hillDecreption(RealMatrix cipherMatrix, RealMatrix inverseMatrix) {
        RealMatrix plaintext = cipherMatrix.multiply(inverseMatrix);
        for (int i = 0; i < plaintext.getRowDimension(); i++) {
            for (int j = 0; j < plaintext.getColumnDimension(); j++) {
                while (plaintext.getEntry(i, j) >= 25) {
                    plaintext.setEntry(i, j, plaintext.getEntry(i, j) - 25);
                }
                while (plaintext.getEntry(i, j) < 0) {
                    plaintext.setEntry(i, j, plaintext.getEntry(i, j) + 25);
                }
            }
        }
        return plaintext;
    }

    private static boolean isPerfectSquare(int toTest) {
        int sqrt = (int) Math.sqrt(toTest);
        return sqrt * sqrt == toTest;
    }

    private static void printMatrix(RealMatrix matrix) {
        int cols = matrix.getColumnDimension();
        int rows = matrix.getRowDimension();
        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                System.out.print("\t" + (int) matrix.getEntry(i, j) + "\t");
            }
            System.out.println();
            System.out.println();
        }
    }

    private static int calculateGCD(int a, int b) {
        while (b != 0) {
            int temp = b;
            b = a % b;
            a = temp;
        }
        return a;
    }

    private static boolean areRelativelyPrime(int a, int b) {
        return calculateGCD(a, b) == 1;
    }

    private static HashMap<String, String> initCodonsTable() {
        HashMap<String, String> codonsTable = new HashMap<>();
        codonsTable.put("GCU", "A-1");
        codonsTable.put("GCC", "A-2");
        codonsTable.put("GCA", "A-3");
        codonsTable.put("GCG", "A-4");
        codonsTable.put("UAA", "B-1");
        codonsTable.put("UAG", "B-2");
        codonsTable.put("UGA", "B-3");
        codonsTable.put("UGU", "C-1");
        codonsTable.put("UGC", "C-2");
        codonsTable.put("GAU", "D-1");
        codonsTable.put("GAC", "D-2");
        codonsTable.put("GAA", "E-1");
        codonsTable.put("GAG", "E-2");
        codonsTable.put("UUU", "F-1");
        codonsTable.put("UUC", "F-2");
        codonsTable.put("GGU", "G-1");
        codonsTable.put("GGC", "G-2");
        codonsTable.put("GGA", "G-3");
        codonsTable.put("GGG", "G-4");
        codonsTable.put("CAU", "H-1");
        codonsTable.put("CAC", "H-2");
        codonsTable.put("AUU", "I-1");
        codonsTable.put("AUC", "I-2");
        codonsTable.put("AUA", "I-3");
        codonsTable.put("AAA", "K-1");
        codonsTable.put("AAG", "K-2");
        codonsTable.put("CUU", "L-1");
        codonsTable.put("CUC", "L-2");
        codonsTable.put("CUA", "L-3");
        codonsTable.put("CUG", "L-4");
        codonsTable.put("AUG", "M-1");
        codonsTable.put("AAU", "N-1");
        codonsTable.put("AAC", "N-2");
        codonsTable.put("UUA", "O-1");
        codonsTable.put("UUG", "O-2");
        codonsTable.put("CCU", "P-1");
        codonsTable.put("CCC", "P-2");
        codonsTable.put("CCA", "P-3");
        codonsTable.put("CCG", "P-4");
        codonsTable.put("CAA", "Q-1");
        codonsTable.put("CAG", "Q-2");
        codonsTable.put("CGU", "R-1");
        codonsTable.put("CGC", "R-2");
        codonsTable.put("CGA", "R-3");
        codonsTable.put("CGG", "R-4");
        codonsTable.put("UCU", "S-1");
        codonsTable.put("UCC", "S-2");
        codonsTable.put("UCA", "S-3");
        codonsTable.put("UCG", "S-4");
        codonsTable.put("ACU", "T-1");
        codonsTable.put("ACC", "T-2");
        codonsTable.put("ACA", "T-3");
        codonsTable.put("ACG", "T-4");
        codonsTable.put("AGA", "U-4");
        codonsTable.put("AGG", "U-4");
        codonsTable.put("GUU", "V-1");
        codonsTable.put("GUC", "V-2");
        codonsTable.put("GUA", "V-3");
        codonsTable.put("GUG", "V-4");
        codonsTable.put("UGG", "W-1");
        codonsTable.put("AGU", "X-1");
        codonsTable.put("AGC", "X-2");
        codonsTable.put("UAU", "Y-1");
        codonsTable.put("UAC", "Z-1");
        return codonsTable;
    }
}