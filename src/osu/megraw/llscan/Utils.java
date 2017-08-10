package osu.megraw.llscan;

import java.io.BufferedReader;
import java.io.FileInputStream;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.List;

/* Copyright (c) 2013 Oregon State University
 *
 * This program is distributed under the terms listed in the
 * LICENSE file included with this software and online at
 * http://oregonstate.edu/research/occd/sites/default/files/permissionstatement.pdf.
 *
 * IN NO EVENT SHALL OREGON STATE UNIVERSITY BE LIABLE TO ANY PARTY FOR DIRECT,
 * INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES, INCLUDING LOST PROFITS,
 * ARISING OUT OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF OREGON
 * STATE UNIVERSITY HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE. OREGON STATE
 * UNIVERSITY SPECIFICALLY DISCLAIMS ANY WARRANTIES, INCLUDING, BUT NOT LIMITED TO,
 * THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * AND ANY STATUTORY WARRANTY OF NON-INFRINGEMENT. THE SOFTWARE PROVIDED HEREUNDER
 * IS ON AN "AS IS" BASIS, AND OREGON STATE UNIVERSITY HAS NO OBLIGATIONS TO
 * PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS. 
 *
 * Contact: megrawm@science.oregonstate.edu
 * */

public class Utils {
	private static void writePWMtoWMMFile(String outfile, double [][]pwm) throws IOException {
        PrintWriter writer = new PrintWriter(outfile);

        // Write out the header of the WMM file
        writer.println("WMM");
        writer.println("PROMOTER");
        // 0 = cutoff (doesn't appear to be used in get-WMM-cutoff)
        // pwm.length = the length of the pwm sequence
        // 5 = number of columns represented as A,C,G,N, and T (in that order) - N isn't found in the PWM matrix format, only four columns exist
        writer.println("0 " + pwm.length + " 5");

        // pwm[0].length = consensus window length which is equal to the pwm sequence length
        // 0 and 0 refer to the consensus offset (which we aren't setting, hence needs to be 0),
        // and the consensus sequence length, respectively - neither are needed so just set to 0
        // + refers to the strand, and all threshold calculations are made with respect to the FWD strand
        writer.println(pwm.length + " 0 0 +");

        for (int j = 0; j < pwm.length; j++) {
            String logPs = "";
            for (int i = 0; i < 5; i++) {
                // Write out 'N' value as '-inf' since we don't need it
                if (i == 3) {
                    logPs += "\t-inf";
                } else {
                    int nt_index = i;

                    // Adjust nt_index to account for the '-inf' value needed for WMM matrix which isn't present in PWM matrix
                    if (nt_index == 4) {
                        nt_index = 3;
                    }
                    // Only prepend a tab to values betyond the first value
                    if (nt_index != 0) {
                        logPs += "\t";
                    }
                    logPs += Double.toString(Math.log(pwm[j][nt_index]));
                }
            }
            writer.println(logPs);
        }
        writer.close();
    }

    public static double getSeqContent(int L, int R, char[] S, int beyondTSS, char[] bases) {

        // Convert window positions to array indices
        int tssLoc = S.length - beyondTSS;
        int winL = L + tssLoc;
        int winR = R + tssLoc;

        if (winL < 0) { winL = 0; }
        if (winR > S.length) { winR = S.length; }

        // Compute gcContent in window
        double seqContent = Utils.seqPercent(String.valueOf(S), winL, winR, bases);

        return seqContent;

    }

    public static double getGCContent(int L, int R, char[] S, int beyondTSS) {
    	
        // Convert window positions to array indices
        int tssLoc = S.length - beyondTSS;
        int winL = L + tssLoc;
        int winR = R + tssLoc;

        if (winL < 0) { winL = 0; }
        if (winR > S.length) { winR = S.length; }

        // Compute gcContent in window
        double gcContent = Utils.gcPercent(String.valueOf(S), winL, winR);

        return gcContent;

    }
	
	public static double getCumScoreInWindow(int L, int R, char[] seq, double[][] pwm, String strand, double[][] B, double[][][] B_M1, double scoreCutOffValue, int nucsAfterTSS) throws BadCharException {
        int w = pwm.length;
        int chrInd, chrInd_prev;
        double theta_freq, bg_freq, score;
        double numLogSum, denLogSum;

        double cumScore = 0.0;

        boolean hitMaskedChar;
        boolean USE_MARKOV_1_BG = false;

        if (B_M1 != null) {
            USE_MARKOV_1_BG = true;
        }

        // Convert window positions to array indices
        int tssLoc = seq.length - nucsAfterTSS;
        int winL = L + tssLoc;
        int winR = R + tssLoc;

        if (winL < 0) { winL = 0; }
        if (winR > seq.length - w) { winR = seq.length - w; }

        // Compute scores for each possible motif site in window
        for (int j = winL; j <= winR; j++) {
            numLogSum = 0.0;
            denLogSum = 0.0;
            // Take product of freq/bg probabilities over motif width
            hitMaskedChar = false;
            for (int mp = 0; mp < w; mp++) {
                chrInd = Setup.charInd(seq[j+mp]);
                if (chrInd < 0) { // encountered a masked character
                    hitMaskedChar = true;
                    break;
                }
                theta_freq = pwm[mp][chrInd];
                if (USE_MARKOV_1_BG) {
                    if (mp == 0) {
                        bg_freq = B[j][chrInd];
                    } else {
                        chrInd_prev = Setup.charInd(seq[j + mp - 1]);
                        bg_freq = B_M1[j][chrInd_prev][chrInd];
                    }
                }
                else {
                    bg_freq = B[j][chrInd];
                }
                numLogSum += Math.log(theta_freq);
                denLogSum += Math.log(bg_freq);
            }
            if (!hitMaskedChar) {
                score = numLogSum - denLogSum;
                if (score > scoreCutOffValue) {
                    // Record cumulative score
                    cumScore += score;
                }
                if (Double.isInfinite(score)) {
                    System.out.println(
                            "WARNING: Encountered an infinite score at position " + j);
                }
            }
        }
        return cumScore;

    }

    public static String reverseComplement(String sequence) {
        StringBuffer b= new StringBuffer();
        for(int i=sequence.length()-1;i>=0;--i) {
            switch(sequence.charAt(i)) {
            case 'A': b.append('T'); break;
            case 'T': b.append('A'); break;
            case 'G': b.append('C'); break;
            case 'C': b.append('G'); break;
            case 'a': b.append('t'); break;
            case 't': b.append('a'); break;
            case 'g': b.append('c'); break;
            case 'c': b.append('g'); break;
            default: b.append(sequence.charAt(i)); break;
            }
        }

        return b.toString();
    }

    public static double seqPercent(String in_sequence, char[] bases) {
        int n = 0;
        for(int i = 0; i < in_sequence.length(); ++i) {
            char base = Character.toUpperCase(in_sequence.charAt(i));
            for (int j = 0; j < bases.length; j++) {
                if (base == bases[j]) {
                    n++;
                    break;
                }
            }
        }
        return (double)(n/(double)in_sequence.length());
    }

    public static double seqPercent(String in_sequence, int L, int R, char[] bases) {
        String sequence = in_sequence.substring(L, R);
        return seqPercent(sequence, bases);
    }

    public static double gcPercent(String in_sequence) {
        int n = 0;
        for(int i = 0; i < in_sequence.length(); ++i) {
            char base = Character.toUpperCase(in_sequence.charAt(i));
            n += (base == 'G' || base == 'C'?1:0);
        }
        return (double)(n/(double)in_sequence.length());
    }

    public static double gcPercent(String in_sequence, int L, int R) {
        String sequence = in_sequence.substring(L, R);
        return gcPercent(sequence);
    }

    public static double[] getBackground(char[][] S, int start, int end) {
        double[] Theta_o = new double[Setup.NCHARS];
        int chrInd;
        double sum;

        // Initialize Theta_o
        for (int k = 0; k < Theta_o.length; k++) {
            Theta_o[k] = 0.0;
        }

        if (start > end || end < 0) {
            return Theta_o;
        }
        if (start < 0) { start = 0; }

        // Tally background bases
        for (int i = 0; i < S.length; i++) {
            if (start > S[i].length) { continue; }
            if (end > S[i].length) { end = S[i].length; }
            for (int j = start; j < end; j++) {
                chrInd = Setup.charInd(S[i][j]);
                if (!(chrInd < 0)) { // base recognized
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

    public static double[][] getM1Background(char[][] S, int start, int end) {
        double[][] Theta_o_M1 = new double[Setup.NCHARS][Setup.NCHARS];
        int chrInd1, chrInd2;
        double[] sum = new double[Setup.NCHARS];

        // Initialize Theta_o_M1
         for (int k = 0; k < Theta_o_M1.length; k++) {
             for (int k2 = 0; k2 < Theta_o_M1[k].length; k2++) {
                 Theta_o_M1[k][k2] = 0.0;
             }
         }

         if (start > end || end < 0) {
             return Theta_o_M1;
         }
         if (start < 0) { start = 0; }

         // Tally background di-nucleotide pairs
         for (int i = 0; i < S.length; i++) {
             if (start > S[i].length) { continue; }
             if (end > S[i].length) { end = S[i].length; }
             for (int j = start; j < end - 1; j++) {
                 chrInd1 = Setup.charInd(S[i][j]);
                 chrInd2 = Setup.charInd(S[i][j+1]);
                 if (!(chrInd1 < 0) && !(chrInd2 < 0))  { // base recognized
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

    public static double[] getBackground(char[] S, int start, int end) {
        double[] Theta_o = new double[Setup.NCHARS];
        int chrInd;
        double sum;

        // Initialize Theta_o
        for (int k = 0; k < Theta_o.length; k++) {
            Theta_o[k] = 0.0;
        }

        if (start > end || end < 0) {
            return Theta_o;
        }
        if (start < 0) { start = 0; }

        // Tally background bases
        if (start > S.length) { return Theta_o; }
        if (end > S.length) { end = S.length; }
        for (int j = start; j < end; j++) {
            chrInd = Setup.charInd(S[j]);
            if (!(chrInd < 0)) { // base recognized
                Theta_o[chrInd] += 1.0;
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

    /**
     * Mitra
     * @param S : one sequence line for a promoter
     * @param BG_WIN : length of the window
     * @return B[S.length][4]: frequencies of A,C,G, and T corresponding to each location in the S
     */
    public static double[][] getWholeSeqLocalBackground(char[] S, int BG_WIN) {
        double[][] B = new double[S.length][Setup.NCHARS];
        if (S.length == 0 || S.length < 2*BG_WIN + 1) { return B; }

        // First calculate and store counts

        // first vector : 
        // Mitra - Count the frequency of each nt (A,C,G,T) within first window of length BG_WIN
        for (int j = 0; j <= BG_WIN; j++) {
            int chrInd = Setup.charInd(S[j]);
            if (!(chrInd < 0)) { // base recognized
                B[0][chrInd] += 1.0;
            }
        }

        // up through first full window
        // Mitra: Shift first window to the right by 1bp and look at the new nt at the right and add 1 to the corresponding bucket (A,C,G,T)
        int firstFullWinIdx = BG_WIN;
        for (int p = 1; p <= firstFullWinIdx; p++) {
            System.arraycopy(B[p-1], 0, B[p], 0, B[p-1].length);
            int endOfWinIdx = p + BG_WIN;
            int chrInd = Setup.charInd(S[endOfWinIdx]);
            if (!(chrInd < 0)) { // base recognized
                B[p][chrInd] += 1.0;
            }
        }

        // all full windows
        // Mitra: Go through all overlapping wins and sum up counts
        int lastFullWinIdx = S.length - 1 - BG_WIN;
        for (int p = firstFullWinIdx + 1; p <= lastFullWinIdx; p++) {
            System.arraycopy(B[p-1], 0, B[p], 0, B[p-1].length);
            int endOfWinIdx = p + BG_WIN;
            int prevStartOfWinIdx = p - BG_WIN - 1;
            int chrInd = Setup.charInd(S[endOfWinIdx]);
            if (!(chrInd < 0)) { // base recognized
                B[p][chrInd] += 1.0;
            }
            chrInd = Setup.charInd(S[prevStartOfWinIdx]);
            if (!(chrInd < 0)) { // base recognized
                B[p][chrInd] -= 1.0;
            }
        }

        // remaining positions with partial windows at end of sequence
        for (int p = lastFullWinIdx + 1; p < S.length; p++) {
            System.arraycopy(B[p-1], 0, B[p], 0, B[p-1].length);
            int prevStartOfWinIdx = p - BG_WIN - 1;
            int chrInd = Setup.charInd(S[prevStartOfWinIdx]);
            if (!(chrInd < 0)) { // base recognized
                B[p][chrInd] -= 1.0;
            }
        }

        // Divide by totals to get frequencies
       for (int i = 0; i < S.length; i++) {
            double sum = 0.0;
            for (int k = 0; k < B[i].length; k++) {
                sum += B[i][k];
            }
            for (int k = 0; k < B[i].length; k++) {
                B[i][k] /= sum;
            }
        }

        return B;
    }

    public static double[][] getM1Background(char[] S, int start, int end) {
        double[][] Theta_o_M1 = new double[Setup.NCHARS][Setup.NCHARS];
        int chrInd1, chrInd2;
        double[] sum = new double[Setup.NCHARS];

        // Initialize Theta_o_M1
         for (int k = 0; k < Theta_o_M1.length; k++) {
             for (int k2 = 0; k2 < Theta_o_M1[k].length; k2++) {
                 Theta_o_M1[k][k2] = 0.0;
             }
         }

         if (start > end || end < 0) {
             return Theta_o_M1;
         }
         if (start < 0) { start = 0; }

         // Tally background di-nucleotide pairs
         if (start > S.length) { return Theta_o_M1; }
         if (end > S.length) { end = S.length; }
         for (int j = start; j < end - 1; j++) {
             chrInd1 = Setup.charInd(S[j]);
             chrInd2 = Setup.charInd(S[j+1]);
             if (!(chrInd1 < 0) && !(chrInd2 < 0))  { // base recognized
                 Theta_o_M1[chrInd1][chrInd2] += 1.0;
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

    public static double[][][] getWholeSeqLocalM1Background(char[] S, int BG_WIN) {
        double[][][] B_M1 = new double[S.length][Setup.NCHARS][Setup.NCHARS];
        if (S.length == 0 || S.length < 2*BG_WIN + 2) { return B_M1; }

        // First calculate and store counts

        // first vector
        for (int j = 0; j <= BG_WIN; j++) {
            int chrInd1 = Setup.charInd(S[j]);
            int chrInd2 = Setup.charInd(S[j+1]);
            if (!(chrInd1 < 0) && !(chrInd2 < 0))  { // bases recognized
                B_M1[0][chrInd1][chrInd2] += 1.0;
            }
        }

        // up through first full window
        int firstFullWinIdx = BG_WIN;
        for (int p = 1; p <= firstFullWinIdx; p++) {
            for (int q = 0; q < B_M1[p].length; q++) {
                System.arraycopy(B_M1[p - 1][q], 0, B_M1[p][q], 0, B_M1[p - 1][q].length);
            }
            int endOfWinIdx = p + BG_WIN;
            int chrInd1 = Setup.charInd(S[endOfWinIdx]);
            int chrInd2 = Setup.charInd(S[endOfWinIdx+1]);
            if (!(chrInd1 < 0) && !(chrInd2 < 0))  { // bases recognized
                B_M1[p][chrInd1][chrInd2] += 1.0;
            }
        }

        // all full windows
        int lastFullWinIdx = S.length - 2 - BG_WIN;
        for (int p = firstFullWinIdx + 1; p <= lastFullWinIdx; p++) {
            for (int q = 0; q < B_M1[p].length; q++) {
                System.arraycopy(B_M1[p - 1][q], 0, B_M1[p][q], 0, B_M1[p - 1][q].length);
            }
            int endOfWinIdx = p + BG_WIN;
            int prevStartOfWinIdx = p - BG_WIN - 1;
            int chrInd1 = Setup.charInd(S[endOfWinIdx]);
            int chrInd2 = Setup.charInd(S[endOfWinIdx+1]);
            if (!(chrInd1 < 0) && !(chrInd2 < 0))  { // bases recognized
                B_M1[p][chrInd1][chrInd2] += 1.0;
            }
            chrInd1 = Setup.charInd(S[prevStartOfWinIdx]);
            chrInd2 = Setup.charInd(S[prevStartOfWinIdx+1]);
            if (!(chrInd1 < 0) && !(chrInd2 < 0))  { // bases recognized
                B_M1[p][chrInd1][chrInd2] -= 1.0;
            }
        }

        // remaining positions with partial windows at end of sequence
        for (int p = lastFullWinIdx + 1; p < S.length; p++) {
            for (int q = 0; q < B_M1[p].length; q++) {
                System.arraycopy(B_M1[p - 1][q], 0, B_M1[p][q], 0, B_M1[p - 1][q].length);
            }
            int prevStartOfWinIdx = p - BG_WIN - 1;
            int chrInd1 = Setup.charInd(S[prevStartOfWinIdx]);
            int chrInd2 = Setup.charInd(S[prevStartOfWinIdx+1]);
            if (!(chrInd1 < 0) && !(chrInd2 < 0))  { // bases recognized
                B_M1[p][chrInd1][chrInd2] -= 1.0;
            }
        }

        // Divide by row totals to get frequencies
        for (int i = 0; i < S.length; i++) {
            double[] sum = new double[Setup.NCHARS];
            for (int k = 0; k < B_M1[i].length; k++) {
                sum[k] = 0.0;
                for (int k2 = 0; k2 < B_M1[i][k].length; k2++) {
                    sum[k] += B_M1[i][k][k2];
                }
            }
            for (int k = 0; k < B_M1[i].length; k++) {
                for (int k2 = 0; k2 < B_M1[i][k].length; k2++) {
                    B_M1[i][k][k2] /= sum[k];
                }
            }
        }

        return B_M1;
    }

    public static double[][] PWMRevCmp(double[][] pwm) {
        double[][] pwm_revcmp = new double[pwm.length][Setup.NCHARS];
        double[][] pwm_tmp = new double[pwm.length][Setup.NCHARS];

        // Complement: flip columns representing A,T and C,G
        for (int mp = 0; mp < pwm.length; mp++) {
            pwm_tmp[mp][Setup.charInd('A')] = pwm[mp][Setup.charInd('T')];
            pwm_tmp[mp][Setup.charInd('T')] = pwm[mp][Setup.charInd('A')];
            pwm_tmp[mp][Setup.charInd('C')] = pwm[mp][Setup.charInd('G')];
            pwm_tmp[mp][Setup.charInd('G')] = pwm[mp][Setup.charInd('C')];
        }

        // Reverse: invert row order
        for (int mp = 0; mp < pwm.length; mp++) {
            int revmp = (pwm.length - 1) - mp;
            System.arraycopy(pwm_tmp[mp], 0, pwm_revcmp[revmp], 0, Setup.NCHARS);
        }

        return pwm_revcmp;
    }

    public static Hashtable <String, Hashtable <Integer, Double>> getCumScore(List<ScanResult> results) {
        return getCumScore(results, "BOTH");
    }

    public static Hashtable <String, Hashtable <Integer, Double>> getCumScore(List<ScanResult> results, String strand) {
        Hashtable <String, Hashtable <Integer, Double>> pwmResults = new Hashtable <String, Hashtable <Integer, Double>>();

        for (int i = 0; i < results.size(); i++) {
            ScanResult result = results.get(i);

            if (!strand.equals("BOTH") && !strand.equals(result.strand)) {
                continue;
            }

            String pwmLabel = result.pwmLabel;

            // Initialize pwmResult cumulative score if we haven't encounterd a scan with this PWM yet
            if (!pwmResults.containsKey(pwmLabel)) pwmResults.put(pwmLabel, new Hashtable <Integer, Double>());

            Hashtable <Integer, Double> pwmResult = pwmResults.get(pwmLabel);

            for (int j = 0; j < result.hitLocs.length; j++) {
                Integer hitLoc = new Integer(result.hitLocs[j]);
                Double hitScore = new Double(result.hitScores[j]);

                // Record cumulative score
                if (pwmResult.containsKey(hitLoc)) {
                    pwmResult.put(hitLoc, new Double(pwmResult.get(hitLoc) + hitScore));
                }
                else {
                    pwmResult.put(hitLoc, hitScore);
                }
            }
        }

        return pwmResults;
    } 

    public static SysCom runSystemCommand (String command) {

        SysCom syscom = new SysCom();

        try {

            Process p = Runtime.getRuntime().exec(command);

            BufferedReader stdInput = new BufferedReader(new InputStreamReader(p.getInputStream()));
            BufferedReader stdError = new BufferedReader(new InputStreamReader(p.getErrorStream()));

            String s, outs;

            // read stdout
            outs = new String();
            while ((s = stdInput.readLine()) != null) {
                outs += s + "\n";
            }
            syscom.stdout = outs.split("\n");

            // read stderr
            outs = new String();
            while ((s = stdError.readLine()) != null) {
                outs += s + "\n";
            }
            syscom.stderr = outs.split("\n");

            syscom.failed = false;
        }
        catch (IOException e) {
        	e.printStackTrace();
            syscom.errmessage = e.toString();
            syscom.failed = true;
        }

        return syscom;
    }
    
    /**
     * @author mitra
     */
	public static  double computeLogLikScoreForStringUsingM0(String seq, int windowLength, double[][] pwmFreq, double[] bgFreq) {
		double logLikScore = 0d;
		String diNts = "ACGT";
		char[] seqArr = seq.toCharArray();
		BigDecimal numLog = new BigDecimal(0);
		BigDecimal dnumLog = new BigDecimal(0);
		int currentCharIdx = diNts.indexOf(seqArr[0]);
		for (int i = 0; i < windowLength; i++) {
			BigDecimal bgFrq = null;
			bgFrq = new BigDecimal(bgFreq[currentCharIdx]);
			currentCharIdx = diNts.indexOf(seqArr[i]);
			dnumLog = dnumLog.add(new BigDecimal(Math.log(bgFrq.doubleValue())));
			
			BigDecimal fgFrq = new BigDecimal(pwmFreq[i][currentCharIdx]);
			numLog = numLog.add(new BigDecimal(Math.log(fgFrq.doubleValue())));
		}
		logLikScore = numLog.subtract(dnumLog).doubleValue();
		return logLikScore;
		
	}

    /**
     * @author mitra
     */
	public static String getRevComplement(String string) {
		String revStr = new StringBuilder(string).reverse().toString();
		 //calculate reverse complement
		revStr = revStr.replace("A", "t").replace("T", "a").
		        replace("C", "g").replace("G", "c").toUpperCase();
       return revStr;
	}
	
    /**
     * @author mitra
     */
	// Load PWM
	public static double[][] getPWM(String inputPWMsFile, String pwmName) {
		List<List<Double>> pwmList = new ArrayList<List<Double>>();
		double[][] pwm = null;
		try {
			BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(new FileInputStream(inputPWMsFile)));
			String line = null;
			boolean foundPWM = false;
			String pLine = "";
			int row = 0;
			while ((line = bufferedReader.readLine()) != null) {
				if (line.startsWith(">")) {
					if (!line.substring(2).equalsIgnoreCase(pwmName)) continue;
					foundPWM = true;
					pLine = ">";
				} else if (line.startsWith("=")){
					if (!foundPWM) continue;
					pLine = "=";
					String[] freqs = line.split("\t");
					pwmList.add(new ArrayList<Double>());
					for (int i = 1; i < freqs.length; i++) {
						pwmList.get(row).add(new BigDecimal(freqs[i]).doubleValue());
					}
					row++;
					
				} else {
					if (foundPWM && pLine.equalsIgnoreCase("="))
						break;
					pLine = "";
				}
			}
			bufferedReader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		pwm = new double[pwmList.size()][4];
		for (int i = 0; i < pwmList.size(); i++) {
			for (int j = 0; j < 4; j++) {
				pwm[i][j] = pwmList.get(i).get(j);
			}
		}
		return pwm;
	}
	
	


}
