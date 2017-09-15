package osu.megraw.llscan;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.Callable;

public class ComputeLoglikScoreROC implements Callable<LoglikScoreResult>{

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
	private int nucsUps = 2000;
	private List<Coordinate> openCoordinatesLeaf;
	private List<Coordinate> openCoordinatesRoot;

	public ComputeLoglikScoreROC(char[] seq, DoubleColReturn scoreCutoffValue,
			int nucsAfterTSS, int bG_WIN, String seqName, PWMReturn pwms, 
			HashMap<String, LefRightPair<Integer, Integer>> windowsHash, int winsWidth, List<Coordinate> openCoordinatesLeaf, List<Coordinate> openCoordinatesRoot, int nucsUpstream) {
		super();
		this.seqName = seqName;
		this.seq = seq;
		this.scoreCutoffValue = scoreCutoffValue;
		this.pwms = pwms;
		this.BG_WIN = bG_WIN;
		this.winHash = windowsHash;
		this.nucsAfterTSS = nucsAfterTSS;
		this.winsWidth = winsWidth;
		this.openCoordinatesLeaf = openCoordinatesLeaf;
		this.openCoordinatesRoot = openCoordinatesRoot;
		nucsUps = nucsUpstream;
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
        	openCoordinatesLeaf = new ArrayList<>();
        }
        
        if (openCoordinatesRoot == null) {
        	System.out.println("Open Chromatin is Null " + seqName);
        	openCoordinatesRoot = new ArrayList<>();
        }
		Collections.sort(openCoordinatesLeaf);
		Collections.sort(openCoordinatesRoot);
		
		HashMap<Integer, List<Coordinate>> fullWinsLeaf = new HashMap<>();
		
		for (Coordinate coordinate : openCoordinatesLeaf) {
//			System.out.println(coordinate);
			int wlow = ((coordinate.getStart() + nucsUps)/winsWidth) + 1;
			int wNext = ((coordinate.getEnd() + nucsUps)/winsWidth);
			if ((coordinate.getEnd() + nucsUps) % winsWidth > 0)
				wNext += 1;
			
			for (int w = wlow; w <= wNext; w++) {
				int left = Math.max(winHash.get(String.valueOf(w)).getL(), coordinate.getStart());
				int right = Math.min(winHash.get(String.valueOf(w)).getR(), coordinate.getEnd());
				
				Coordinate ocCoord = new Coordinate("", left, right);
				List<Coordinate> cList = fullWinsLeaf.get(w);
				if (cList == null) {
					cList = new ArrayList<>();
					cList.add(ocCoord);
				} else 
					cList.add(ocCoord);
				fullWinsLeaf.put(w, cList);
			}
		}
		
		HashMap<Integer, List<Coordinate>> fullWinsRoot = new HashMap<>();
		
		for (Coordinate coordinate : openCoordinatesRoot) {
			int wlow = ((coordinate.getStart() + nucsUps)/winsWidth) + 1;
			int wNext = ((coordinate.getEnd() + nucsUps)/winsWidth);
			if ((coordinate.getEnd() + nucsUps) % winsWidth > 0)
				wNext += 1;
			
			for (int w = wlow; w <= wNext; w++) {
				int left = Math.max(winHash.get(String.valueOf(w)).getL(), coordinate.getStart());
				int right = Math.min(winHash.get(String.valueOf(w)).getR(), coordinate.getEnd());
				
				Coordinate ocCoord = new Coordinate("", left, right);
				List<Coordinate> cList = fullWinsRoot.get(w);
				if (cList == null) {
					cList = new ArrayList<>();
					cList.add(ocCoord);
				} else 
					cList.add(ocCoord);
				fullWinsRoot.put(w, cList);
			}
		}
        
        for (String strand : strands) {
	        for (int winNo = 1; winNo <= winHash.size(); winNo++) {
				for (String pwmID : pwms.labels) {
					String featureIdLeaf = pwmID + "_" + strand + "_" + winNo + "_ROC_LEAF_" + winsWidth ;
					String featureIdRoot = pwmID + "_" + strand + "_" + winNo + "_ROC_ROOT_" + winsWidth ;
					
					int pwmIdx = pwms.labelIndex.get(pwmID);
					double[][] ThetaT;
				    if (strand.equals("FWD")) {
				        ThetaT = pwms.pwms[pwmIdx];
				    } else { // REV
				        ThetaT = pwms.revcmp_pwms[pwmIdx];
				    }

					List<Coordinate> cList = fullWinsLeaf.get(winNo);
					double winScore = 0.0d;
		        	if (cList != null) {
		        		for (Coordinate coordinate : cList) {
						    double score = Utils.getCumScoreInWindow(coordinate.getStart(), coordinate.getEnd(), seq, ThetaT, strand,
		                            B, B_M1, scoreCutoffValue.values[pwmIdx], nucsAfterTSS);
						    winScore += score;
						}
		        	}
		        	featureHash.put(featureIdLeaf, winScore);
		        	
		        	cList = fullWinsRoot.get(winNo);
		        	winScore = 0.0d;
		        	if (cList != null) {
		        		for (Coordinate coordinate : cList) {
		        			double score = Utils.getCumScoreInWindow(coordinate.getStart(), coordinate.getEnd(), seq, ThetaT, strand,
		        					B, B_M1, scoreCutoffValue.values[pwmIdx], nucsAfterTSS);
		        			winScore += score;
		        		}
		        	}
		        	featureHash.put(featureIdRoot, winScore);
				}
			}
        }
        
        // Calculate GC Content (-100 to +100 related to TSS)
//		int GC_WIN = 100;
//
//        double gccontent = Utils.getGCContent( -1 * GC_WIN, GC_WIN, seq, nucsAfterTSS);
//        String featureId = "GCcontent";
//        featureHash.put(featureId, gccontent);
//        if (gccontent != 0) numberNonzeros++;
//
//        double CAcontent = Utils.getSeqContent(-1 * GC_WIN, GC_WIN, seq, nucsAfterTSS, new char[]{'C', 'A'});
//        featureId = "CAcontent";
//        featureHash.put(featureId, CAcontent);
//        if (gccontent != 0) numberNonzeros++;
//
//        double GAcontent = Utils.getSeqContent(-1 * GC_WIN, GC_WIN, seq, nucsAfterTSS, new char[]{'G', 'A'});
//        featureId = "GAcontent";
//        featureHash.put(featureId, GAcontent);
//        if (gccontent != 0) numberNonzeros++;
            
        sampleName = seqName + "_" + pos;
        LoglikScoreResult loglikScoreResult = new LoglikScoreResult(featureHash, sampleName, numberNonzeros);
		return loglikScoreResult;
	}

}
