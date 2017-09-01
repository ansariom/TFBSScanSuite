package osu.megraw.llscan;

public class RefSeqData {
    public String id = null; // Chromosome
    public int start = -1; // promoter start
    public String strand = "+";

    public RefSeqData(String id, int start) {
        this.id = id;
        this.start = start;
    }
    
    public RefSeqData(String id, int start, String strand) {
        this.id = id;
        this.start = start;
        this.strand = strand;
    }


}
