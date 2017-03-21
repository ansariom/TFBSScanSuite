package osu.megraw.llscan;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.DecimalFormat;

/**
 * Copyright (C) 2010  Molly Megraw
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: molly.megraw@duke.edu
 *
 */
public class Print {

    public static DecimalFormat df0 = new DecimalFormat("####0");
    public static DecimalFormat df2 = new DecimalFormat("####0.00");
    public static DecimalFormat df3 = new DecimalFormat("####0.000");
    public static DecimalFormat df4 = new DecimalFormat("####0.0000");
    public static DecimalFormat df6 = new DecimalFormat("####0.000000");
    public static DecimalFormat df10 = new DecimalFormat("####0.0000000000");
    public static DecimalFormat df14 = new DecimalFormat("####0.00000000000000");

    public Print() {
    }

    public static void printArray(int[] in_array, String str) {

        System.out.println(str);
        for (int i = 0; i < in_array.length; i++) {
            System.out.println(in_array[i]);
        }
        System.out.println();

    }
    public static void printArray(double[] in_array, String str) {

        System.out.println(str);
        for (int i = 0; i < in_array.length; i++) {
            System.out.println(df4.format(in_array[i]));
        }
        System.out.println();

    }

    public static void printArray(String[] in_array, String str) {

        System.out.println(str);
        for (int i = 0; i < in_array.length; i++) {
            System.out.println(in_array[i]);
        }
        System.out.println();

    }

    public static void printArray(char[] in_array, String str) {

        System.out.println(str);
        for (int i = 0; i < in_array.length; i++) {
            System.out.println(in_array[i]);
        }
        System.out.println();

    }

    public static void printMatrix(int[][] in_mat, String str) {

        System.out.println(str);
        for (int i = 0; i < in_mat.length; i++) {
            System.out.println("Row " + i + ":");
            for (int j = 0; j < in_mat[i].length; j++) {
                System.out.print(in_mat[i][j] + " ");
            }
            System.out.println();
            System.out.println();
        }
        System.out.println();

    }


    public static void printMatrix(double[][] in_mat, String str) {

        System.out.println(str);
        for (int i = 0; i < in_mat.length; i++) {
            System.out.println("Row " + i + ":");
            for (int j = 0; j < in_mat[i].length; j++) {
                System.out.print(df4.format(in_mat[i][j]) + " ");
            }
            System.out.println();
            System.out.println();
        }
        System.out.println();

    }

    public static void printMatrix(char[][] in_mat, String str) {

        System.out.println(str);
        for (int i = 0; i < in_mat.length; i++) {
            System.out.println("Row " + i + ":");
            for (int j = 0; j < in_mat[i].length; j++) {
                System.out.print(in_mat[i][j] + " ");
            }
            System.out.println();
            System.out.println();
        }
        System.out.println();

    }

    public static void printMatrixToFile(double[][] in_mat, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_mat.length; i++) {
                for (int j = 0; j < in_mat[i].length; j++) {
                    outputFile.print(Print.df4.format(in_mat[i][j]) + " ");
                }
                outputFile.println();
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printMatrixToFile(double[][] in_mat, String filename, DecimalFormat in_df) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_mat.length; i++) {
                for (int j = 0; j < in_mat[i].length; j++) {
                    outputFile.print(in_df.format(in_mat[i][j]) + " ");
                }
                outputFile.println();
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printMatrixToFile(int[][] in_mat, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_mat.length; i++) {
                for (int j = 0; j < in_mat[i].length; j++) {
                    outputFile.print(in_mat[i][j] + " ");
                }
                outputFile.println();
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printMatrixToFile(char[][] in_mat, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_mat.length; i++) {
                outputFile.println("Row " + i + ":");
                for (int j = 0; j < in_mat[i].length; j++) {
                    outputFile.print(in_mat[i][j]);
                }
                outputFile.println();
                outputFile.println();
            }
            outputFile.println();
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printCharMatrix(char[][] in_mat) {
        for (int i = 0; i < in_mat.length; i++) {
                System.out.println("Row " + i + ":");
                for (int j = 0; j < in_mat[i].length; j++) {
                    System.out.print(in_mat[i][j]);
                }
                System.out.println();
                System.out.println();
            }
            System.out.println();
    }

    public static void printCharMatrixToFile(char[][] in_mat, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_mat.length; i++) {
                outputFile.println("Row " + i + ":");
                for (int j = 0; j < in_mat[i].length; j++) {
                    outputFile.print(in_mat[i][j]);
                }
                outputFile.println();
                outputFile.println();
            }
            outputFile.println();
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printLabeledMatrixToFile(int[][] in_mat, String[] labels, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_mat.length; i++) {
                outputFile.print(labels[i] + "\t");
                for (int j = 0; j < in_mat[i].length; j++) {
                    outputFile.print(in_mat[i][j] + " ");
                }
                outputFile.println();
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printLabeledMatrixToFile(String[][] in_mat, String[] labels, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_mat.length; i++) {
                outputFile.print(labels[i] + "\t");
                for (int j = 0; j < in_mat[i].length; j++) {
                    outputFile.print(in_mat[i][j] + " ");
                }
                outputFile.println();
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void appendMatrixToFile(PrintWriter pw, double[][] in_mat, String str, DecimalFormat in_df) {

        pw.println(str);
        for (int i = 0; i < in_mat.length; i++) {
            for (int j = 0; j < in_mat[i].length; j++) {
                pw.print(in_df.format(in_mat[i][j]) + " ");
            }
            pw.println();
            pw.println();
        }
        pw.println();
        pw.flush();
    }

    public static void appendPWMToFile(PrintWriter pw, double[][] in_mat, String str, DecimalFormat in_df) {

        pw.println(str);
        for (int i = 0; i < in_mat.length; i++) {
            pw.print("=");
            for (int j = 0; j < in_mat[i].length; j++) {
                pw.print("\t" + in_df.format(in_mat[i][j]));
            }
            pw.println();
        }
        pw.println();
        pw.flush();
    }

    public static void printMatrixTranspose(int[][] in_mat, String str) {

        System.out.println(str);
        for (int j = 0; j < in_mat[0].length; j++) {
            for (int i = 0; i < in_mat.length; i++) {
                System.out.print(in_mat[i][j] + " ");
            }
            System.out.println();
        }
        System.out.println();

    }

    public static void printArrayToFile(int[] in_array, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_array.length; i++) {
                outputFile.print(in_array[i] + "\n");
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void appendArrayToFile( PrintWriter pw, double[] in_array, String str, DecimalFormat in_df) {

        pw.println(str);
        for (int i = 0; i < in_array.length; i++) {
            pw.println(in_df.format(in_array[i]));
        }
        pw.println();
        pw.flush();
    }

    public static void printLabeledArrayToFile(int[] in_array, String[] in_labels, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_array.length; i++) {
                outputFile.print(in_labels[i] + "\t" + in_array[i] + "\n");
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printArrayToFile(double[] in_array, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_array.length; i++) {
                outputFile.print(df4.format(in_array[i]) + "\n");
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printArrayToFile(double[] in_array, String filename, DecimalFormat in_df) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_array.length; i++) {
                outputFile.print(in_df.format(in_array[i]) + "\n");
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printAsFastaFile(String[] in_headers, char[][] in_seqs, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_headers.length; i++) {
                outputFile.println("> " + in_headers[i]);
                for (int j = 0; j < in_seqs[i].length; j++) {
                    outputFile.print(in_seqs[i][j]);
                }
                outputFile.print("\n");
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printAsFastaFile(String[] in_headers, char[][] in_seqs, File in_file) {
        try {
            FileWriter outFile = new FileWriter(in_file);
            PrintWriter outputFile = new PrintWriter(outFile);

            for (int i = 0; i < in_headers.length; i++) {
                outputFile.println("> " + in_headers[i]);
                for (int j = 0; j < in_seqs[i].length; j++) {
                    outputFile.print(in_seqs[i][j]);
                }
                outputFile.print("\n");
            }
            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printPWMsToFile(PWMReturn pwms, String filename, DecimalFormat in_df) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            double[][] in_mat;
            for (int nmat = 0; nmat < pwms.pwms.length; nmat++) {
                in_mat = pwms.pwms[nmat];

                outputFile.print("> " + pwms.labels[nmat]);
                outputFile.println();
                for (int i = 0; i < in_mat.length; i++) {
                    outputFile.print("=");
                    for (int j = 0; j < in_mat[i].length; j++) {
                        /*
                        if (j > 0) {
                            outputFile.print("\t");
                        }
                        */
                        outputFile.print("\t");
                        outputFile.print(in_df.format(in_mat[i][j]));
                    }
                    outputFile.println();
                }
                outputFile.println();
            }

            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

    public static void printFreqsToFile(FreqReturn freqs, String filename) {
        try {
            File selectedFile = new File(filename);
            FileWriter outFile = new FileWriter(selectedFile);
            PrintWriter outputFile = new PrintWriter(outFile);

            int[][] in_mat;
            for (int nmat = 0; nmat < freqs.freqs.length; nmat++) {
                in_mat = freqs.freqs[nmat];

                outputFile.print(freqs.headers[nmat]);
                outputFile.println();
                for (int i = 0; i < in_mat.length; i++) {
                    outputFile.print("=");
                    for (int j = 0; j < in_mat[i].length; j++) {
                        outputFile.print("\t" + in_mat[i][j]);
                    }
                    outputFile.println();
                }
                outputFile.println();
            }

            outputFile.close();
        }
        catch (IOException ioe) {
            System.out.println("Could not print output file " + "filename");
        }
    }

}
