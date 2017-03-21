package osu.megraw.llscan;

import java.io.File;

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
public class Background {
    // Stores sequnece identifer and background models for local sequence background distributions 
    public String seqLabel;  // Not used especially - just added here in case needed for future code updates
    public double[][] B; //0th order background model
    public double[][][] B_M1; // 1st order background model;
    public boolean useLocal;
    public boolean useEqual;
    public boolean useM1;

    // No sequence info was sent in - assume we are using a sequence expecting an 
    // equal background distribution for all bases
    public Background() {
        this.useLocal = false;
        this.useEqual = true;
        this.useM1 = false; // Since everything is equal, we don't need to worry about 1st order
    }

    // Background sequence information was sent in - assume non-equal distribution of bases
    public Background (String seqLabel, double[][] B, double[][][] B_M1, boolean useLocal) {
        this.seqLabel = seqLabel;
        this.B = B;
        this.B_M1 = B_M1;
        this.useLocal = useLocal;
        this.useEqual = false;

        // user did not supply 1st order background model
        if (this.B_M1 == null) {
            this.useM1 = false;

        // user did supply 1st order background model - use it in loglikelihood calculations
        } else {
            this.useM1 = true;
        }
    }

    public Background (String seqLabel, double[][] B, double[][][] B_M1) {
       this(seqLabel, B, B_M1, true); // By default, object assumes you are using local distribution
    }

    public static void main(String[] args) throws java.io.IOException {
        Background background = new Background();

        String usage = "usage: java -jar tfbs_scan.jar Background <FASTA FILE> <OUTFILE1> <OUTFILE2>\n\n" +
                       "The Background utility generates 0th and 1st order background models base\n" +
                       "on the <FASTA FILE> input sequences. It outputs the 0th and 1st bacgkround\n" +
                       "models to <OUTFILE1> and <OUTFILE2> respectively.  These background models\n" +
                       "can then be used for the Scan, ROEFinder, and GenFeatures scanning programs.\n";

        if (args.length < 3) {
            System.out.println(usage);
            return;
        }

        // Define input files
        String inFile = args[0];

        // Define output files
        String outFile = args[1];
        String outFile2 = args[2];

        // Read input file in Fasta format
        FastaReturn fa = Load.loadFastaFile(new File(inFile));
        char[][] S = new char[fa.lines.length][];
        for (int i = 0; i < fa.lines.length; i++) {
            S[i] = fa.lines[i].toCharArray();
        }

        double bg_m0[] = getBackground_M0(S);
        Print.printArrayToFile(bg_m0, outFile);

        double bg_m1[][] = getBackground_M1(S);
        Print.printMatrixToFile(bg_m1, outFile2);
    }

    public static double[] getBackground_M0(char[][] S) {
        double[] Theta_o = new double[Setup.NCHARS];
        int chrInd;
        double sum;

        // Initialize Theta_o
        for (int k = 0; k < Theta_o.length; k++) {
            Theta_o[k] = 0.0;
        }

        // Tally background bases
        for (int i = 0; i < S.length; i++) {
            for (int j = 0; j < S[i].length; j++) {
                chrInd = Setup.charInd(S[i][j]);
                if (!(chrInd < 0)) { // recognized character
                    Theta_o[chrInd] += 1.0;
                }
            }
        }

        // Divide by total
        sum = 0.0;
        for (int k = 0; k < Theta_o.length; k++) {
            sum += Theta_o[k];
        }
        for (int k = 0; k < Theta_o.length; k++) {
            Theta_o[k] /= sum;
        }

        return Theta_o;
    }

    public static double[][] getBackground_M1(char[][] S) {
        double[][] Theta_o_M1 = new double[Setup.NCHARS][Setup.NCHARS];

        int chrInd1, chrInd2;
        double[] sum = new double[Setup.NCHARS];

        // Initialize Theta_o_M1
        for (int k = 0; k < Theta_o_M1.length; k++) {
            for (int k2 = 0; k2 < Theta_o_M1[k].length; k2++) {
                Theta_o_M1[k][k2] = 0.0;
            }
        }

        // Tally background di-nucleotide pairs
        for (int i = 0; i < S.length; i++) {
            for (int j = 0; j < S[i].length - 1; j++) {
                chrInd1 = Setup.charInd(S[i][j]);
                chrInd2 = Setup.charInd(S[i][j+1]);
                if (!(chrInd1 < 0) && !(chrInd2 < 0))  { // recognized characters
                    Theta_o_M1[chrInd1][chrInd2] += 1.0;
                }
            }
        }

        // Divide by row totals
        for (int k = 0; k < Theta_o_M1.length; k++) {
            sum[k] = 0.0;
            for (int k2 = 0; k2 < Theta_o_M1[k].length; k2++) {
                sum[k] += Theta_o_M1[k][k2];
            }
        }
        for (int k = 0; k < Theta_o_M1.length; k++) {
            for (int k2 = 0; k2 < Theta_o_M1[k].length; k2++) {
                Theta_o_M1[k][k2] /= sum[k];
            }
        }

        return Theta_o_M1;
    }
}
