package osu.megraw.llscan;

public class ScanResult {
    String pwmLabel;
    String seqLabel;
    String strand;
    int[] hitLocs;
    double[] hitScores;

    public ScanResult(int[] hitLocs, double[] hitScores, String pwmLabel, String seqLabel, String strand) {
        this.hitLocs = hitLocs;
        this.hitScores = hitScores;
        this.pwmLabel = pwmLabel;
        this.seqLabel = seqLabel;
        this.strand = strand;
    }
}
