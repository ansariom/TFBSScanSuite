package osu.megraw.llscan;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.Callable;

public class ComputeLoglikScoreROEROC implements Callable<LoglikScoreResult>{

	private String sampleName;
	private int BG_WIN;
	private char[] seq;
	DoubleColReturn scoreCutoffValue;
	int nucsAfterTSS;
	private PWMReturn pwms;
	String[] strands = {"FWD", "REV"};
	private String seqName;
	long numberNonzeros = 0;
	private HashMap<String, List<Coordinate>> openCoordinatesLeaf;
	private HashMap<String, List<Coordinate>> openCoordinatesRoot;
	
	private String[] winVarNames;
	private int[] winVarL;
	private int[] winVarR;
	private int[] winVarPWMIndex;
	private String[] winVarStrand;
	private int nWinVars;


	
	public ComputeLoglikScoreROEROC(char[] seq, DoubleColReturn scoreCutOffs, int nucsAfterTSS, int bG_WIN,
			String seqName, PWMReturn pwms, String[] winVarNames, int[] winVarL, int[] winVarR, int[] winVarPWMIndex,
			String[] winVarStrand, HashMap<String, List<Coordinate>> leafCoordsHash,
			HashMap<String, List<Coordinate>> rootCoordsHash) {
		super();
		this.seqName = seqName;
		this.seq = seq;
		this.scoreCutoffValue = scoreCutOffs;
		this.pwms = pwms;
		this.BG_WIN = bG_WIN;
		this.nucsAfterTSS = nucsAfterTSS;
		this.openCoordinatesLeaf = leafCoordsHash;
		this.openCoordinatesRoot = rootCoordsHash;
		this.winVarNames = winVarNames;
		this.winVarL = winVarL;
		this.winVarR = winVarR;
		this.winVarPWMIndex = winVarPWMIndex;

	}

	@Override
	public LoglikScoreResult call() throws Exception {
		// Compute and store local background values for this sequence
		double[][] B = Utils.getWholeSeqLocalBackground(seq, BG_WIN);
		double[][][] B_M1 = Utils.getWholeSeqLocalM1Background(seq, BG_WIN);

		int pos = 0;
        
        // Stores featureWinId and corresponding cumScore value
        HashMap<String, Double> featureHash = new HashMap<>();

        if (openCoordinatesLeaf == null) {
        	System.out.println("Open Chromatin is Null " + seqName);
        	openCoordinatesLeaf = new HashMap<>();
        }
        
        if (openCoordinatesRoot == null) {
        	System.out.println("Open Chromatin is Null " + seqName);
        	openCoordinatesRoot = new HashMap<>();
        }
        
		for (int wv = 0; wv < nWinVars; wv++) {
		    int pwmIdx = winVarPWMIndex[wv];
		    String strand = winVarStrand[wv];

		    double[][] ThetaT;
		    // SHAWN
		   	//System.out.print("Label: " + pwms.labels[pwmIdx] + " " + seqName + " ");
		    if (strand.equals("FWD")) {
		        ThetaT = pwms.pwms[pwmIdx];
		    } else { // REV
		        ThetaT = pwms.revcmp_pwms[pwmIdx];
		    }
		    
		    List<Coordinate> leafCoorList = openCoordinatesLeaf.get(winVarNames[wv]);
		    List<Coordinate> rootCoorList = openCoordinatesRoot.get(winVarNames[wv]);

		    if (leafCoorList == null) {
		    	leafCoorList = new ArrayList<>();
		    	System.out.println("chromatin is closed < " + seqName + " - " + winVarNames[wv]);
		    }
		    
		    if (rootCoorList == null) {
		    	rootCoorList = new ArrayList<>();
			    System.out.println("chromatin is closed < " + seqName + " - " + winVarNames[wv]);
		    }
		    
		    Collections.sort(leafCoorList);
		    Collections.sort(rootCoorList);
		    
			String leafFeatureId = winVarNames[wv] + "_LEAF";
			double leafWinSum = 0.0d;
			for (Coordinate coordinate : leafCoorList) {
			    double score = Utils.getCumScoreInWindow(coordinate.getStart(), coordinate.getEnd(), seq, ThetaT, strand,
                        B, B_M1, scoreCutoffValue.values[pwmIdx], nucsAfterTSS);
			    leafWinSum += score;
			}
        	featureHash.put(leafFeatureId, leafWinSum);

        	String rootFeatureId = winVarNames[wv] + "_ROOT";
        	double rootWinSum = 0.0d;
        	for (Coordinate coordinate : rootCoorList) {
        		double score = Utils.getCumScoreInWindow(coordinate.getStart(), coordinate.getEnd(), seq, ThetaT, strand,
        				B, B_M1, scoreCutoffValue.values[pwmIdx], nucsAfterTSS);
        		rootWinSum += score;
        	}
        	featureHash.put(rootFeatureId, rootWinSum);
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
            
        sampleName = seqName + "_0";
        LoglikScoreResult loglikScoreResult = new LoglikScoreResult(featureHash, sampleName, numberNonzeros);
		return loglikScoreResult;
	}

}
