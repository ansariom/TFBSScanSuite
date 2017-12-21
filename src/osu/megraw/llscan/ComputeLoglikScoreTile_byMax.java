package osu.megraw.llscan;

import java.util.HashMap;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.Random;
import java.util.concurrent.Callable;

public class ComputeLoglikScoreTile_byMax implements Callable<LoglikScoreResult>{

	private String sampleName;
	private int BG_WIN;
	private char[] seq;
	DoubleColReturn scoreCutoffValue;
	int nucsAfterTSS;
	private PWMReturn pwms;
	private HashMap<String, LefRightPair<Integer, Integer>> winHash;
	String[] strands = {"FWD", "REV"};
	private String seqName;
	long numberNonzeros = 0;
	private int winsWidth = 100;
	DoubleColReturn maxScores;

	public ComputeLoglikScoreTile_byMax(char[] seq, DoubleColReturn scoreCutoffValue,
			int nucsAfterTSS, int bG_WIN, String seqName, PWMReturn pwms, 
			HashMap<String, LefRightPair<Integer, Integer>> windowsHash, int winsWidth, DoubleColReturn maxScores) {
		super();
		this.seqName = seqName;
		this.seq = seq;
		this.scoreCutoffValue = scoreCutoffValue;
		this.pwms = pwms;
		this.BG_WIN = bG_WIN;
		this.winHash = windowsHash;
		this.nucsAfterTSS = nucsAfterTSS;
		this.winsWidth = winsWidth;
		this.maxScores = maxScores;
	}
	
	@Override
	public LoglikScoreResult call() throws Exception {
		// Compute and store local background values for this sequence
		double[][] B = Utils.getWholeSeqLocalBackground(seq, BG_WIN);
		double[][][] B_M1 = Utils.getWholeSeqLocalM1Background(seq, BG_WIN);

		int pos = 0;
		int nmats = pwms.labels.length;
        int numFeatures = ((2 * winHash.size()) * pwms.labels.length) + 3; // number of windows, plus GC, CA, GA contents
        double L = Double.valueOf(winsWidth).doubleValue();
        // Stores featureWinId and corresponding cumScore value
        HashMap<String, Double> featureHash = new HashMap<>();
        
        for (String strand : strands) {
	        for (int winNo = 1; winNo <= winHash.size(); winNo++) {
	        	int left = winHash.get(String.valueOf(winNo)).getL();
	        	int right = winHash.get(String.valueOf(winNo)).getR();
	        	
	        	for (int i = 0; i < nmats; i++) {
	        		String pwmID = pwms.labels[i];
					String featureId = pwmID + "_" + strand + "_" + winNo + "_tile" + winsWidth ;
					
					int pwmIdx = pwms.labelIndex.get(pwmID);
					double[][] ThetaT;
				    if (strand.equals("FWD")) {
				        ThetaT = pwms.pwms[pwmIdx];
				    } else { // REV
				        ThetaT = pwms.revcmp_pwms[pwmIdx];
				    }
				    
				    //Compute score
				    double score = Utils.getCumScoreInWindow(left, right, seq, ThetaT, strand,
                            B, B_M1, scoreCutoffValue.values[pwmIdx], nucsAfterTSS);
				    if (score != 0) numberNonzeros++;
				    score = (score / (L * maxScores.values[pwmIdx])) * 10;
				    featureHash.put(featureId, score);
				}
			}
        }
        
        // Calculate GC Content (-100 to +100 related to TSS)
		int GC_WIN = 100;

        double gccontent = Utils.getGCContent( -1 * GC_WIN, GC_WIN, seq, nucsAfterTSS);
        String featureId = "GCcontent";
        featureHash.put(featureId, gccontent);
        if (gccontent != 0) numberNonzeros++;

        double CAcontent = Utils.getSeqContent(-1 * GC_WIN, GC_WIN, seq, nucsAfterTSS, new char[]{'C', 'A'});
        featureId = "CAcontent";
        featureHash.put(featureId, CAcontent);
        if (gccontent != 0) numberNonzeros++;

        double GAcontent = Utils.getSeqContent(-1 * GC_WIN, GC_WIN, seq, nucsAfterTSS, new char[]{'G', 'A'});
        featureId = "GAcontent";
        featureHash.put(featureId, GAcontent);
        if (gccontent != 0) numberNonzeros++;
            
        sampleName = seqName + "_" + pos;
        LoglikScoreResult loglikScoreResult = new LoglikScoreResult(featureHash, sampleName, numberNonzeros);
		return loglikScoreResult;
	}

}
