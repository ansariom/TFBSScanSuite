package osu.megraw.llscan;

import java.util.concurrent.Callable;

public class BGRunner implements Callable<Background>{

	char[] S;
	String seqLabel;
	int BG_WIN = 1;
	
	public BGRunner(char[] S, int BG_WIN, String seqLabel) {
		this.S = S;
		this.BG_WIN = BG_WIN;
		this.seqLabel = seqLabel;
	}
	
	@Override
	public Background call() throws Exception {
        double[][] B = Utils.getWholeSeqLocalBackground(S, BG_WIN);
        double[][][] B_M1 = Utils.getWholeSeqLocalM1Background(S, BG_WIN);
        return new Background(seqLabel, B, B_M1);
	}
	

}
