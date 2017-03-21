package osu.megraw.llscan;

import java.util.HashMap;

public class LoglikScoreResult {
    long numberNonzeros = 0;
	String sampleName;
	HashMap<String, Double> featureHash;
	
	public LoglikScoreResult() {
	}
	
	public LoglikScoreResult(HashMap<String, Double> featureHash, String sampleName, long numberNonzeros) {
		this.sampleName = sampleName;
        this.numberNonzeros = numberNonzeros;
        this.featureHash = featureHash;
	}

}
