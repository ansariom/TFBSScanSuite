package osu.megraw.llscan;

import java.util.Vector;
import java.util.concurrent.Callable;

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

public class ScanRunner implements Callable<ScanResult> {
    char[] S;
    double[][] in_ThetaT;
    String strand;
    Background BG;
    double minScore;
    int beyondTSS;
    String pwmLabel;
    String seqLabel;

    public ScanRunner (char[] S, double[][] in_ThetaT, String strand, Background BG, double minScore, int beyondTSS,
                       String pwmLabel, String seqLabel) {
        this.S = S;
        this.in_ThetaT = in_ThetaT;
        this.strand = strand;
        this.BG = BG;
        this.minScore = minScore;
        this.beyondTSS = beyondTSS;
        this.pwmLabel = pwmLabel;
        this.seqLabel = seqLabel;
    }
    
    public ScanRunner (char[] S, double[][] in_ThetaT, String strand, double minScore, int beyondTSS,
            String pwmLabel, String seqLabel, int BG_WIN) {
		this.S = S;
		this.in_ThetaT = in_ThetaT;
		this.strand = strand;
		double[][] B = Utils.getWholeSeqLocalBackground(S, BG_WIN);
        double[][][] B_M1 = Utils.getWholeSeqLocalM1Background(S, BG_WIN);
        if (BG_WIN > 0) {
        	this.BG = new Background(seqLabel, B, B_M1);
        } else {
        	this.BG = new Background();
        }
		this.minScore = minScore;
		this.beyondTSS = beyondTSS;
		this.pwmLabel = pwmLabel;
		this.seqLabel = seqLabel;
    }
    

    public ScanResult call() {
    	System.out.println("ScanRunner is called for : " + pwmLabel + ", " + seqLabel);
        return this.scan();
    }

    public ScanResult scan() {
        double[][] ThetaT;
        if (strand.equals("FWD")) {
            ThetaT = in_ThetaT;
        }
        else { // REV
            ThetaT = Utils.PWMRevCmp(in_ThetaT);
        }

        int w = ThetaT.length;
        int chrInd, chrInd_prev;
        double theta_freq, bg_freq, score;
        double numLogSum, denLogSum;
        double prop;
        int[] hitLocs;
        double[] hitScores;
        Vector locVec = new Vector();
        Vector scoreVec = new Vector();
        boolean acceptMotif;

        // When bgIndex remains 0 indicates we are *NOT* using a local background sequence
        int bgIndex = 0;
        for (int j = 0; j < (S.length - w + 1); j++) {
            // Take product of freq/bg probabilities over motif width
            numLogSum = 0.0;
            denLogSum = 0.0;
            acceptMotif = true;

            // Update the background index to change to current background
            // sequence location if we are doing local background distributions,
            // otherwise assume we are only using one background distribution, located
            // at index 0
            if (BG.useLocal) {
                bgIndex = j;
            }

            for (int mp = 0; mp < w; mp++) {
                chrInd = Setup.charInd(S[j+mp]);
                if (chrInd < 0) { // encountered an unrecognized character
                    acceptMotif = false;
                    break;
                }

                theta_freq = ThetaT[mp][chrInd];

                if (!BG.useEqual) {
                    if (BG.useM1) {
                        if (mp == 0) {
                            bg_freq = BG.B[bgIndex][chrInd];
                        } else {
                            chrInd_prev = Setup.charInd(S[j + mp - 1]);
                            bg_freq = BG.B_M1[bgIndex][chrInd_prev][chrInd];
                        }
                    }
                    else {
                        bg_freq = BG.B[bgIndex][chrInd];
                    }
                } else {
                    bg_freq = 0.25;
                }

                numLogSum += Math.log(theta_freq);
                denLogSum += Math.log(bg_freq);
            }

            // Only consider motif if all characters are recognized
            if (acceptMotif) {
                score = numLogSum - denLogSum;
                if (score > minScore) {
                    int pos = j;
                    int tssLoc = S.length - beyondTSS;
                    int hitLoc = pos - tssLoc;

                    // Record hit location and score
                    locVec.addElement(new Integer(hitLoc));
                    scoreVec.addElement(new Double(score));
                }
                if (Double.isInfinite(score)) {
                    System.out.println("WARNING: Encountered an infinite score at position (" + j + ")");
                }
            }
        }

        hitLocs = new int[locVec.size()];
        hitScores = new double[scoreVec.size()];
        for (int i = 0; i < locVec.size(); i++) {
            hitLocs[i] = ((Integer) locVec.elementAt(i)).intValue();
            hitScores[i] = ((Double) scoreVec.elementAt(i)).doubleValue();
        }

        return new ScanResult(hitLocs, hitScores, pwmLabel, seqLabel, strand);
    }
}
