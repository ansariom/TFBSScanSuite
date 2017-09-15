package osu.megraw.llscan;

/**
 * @author mitra
 *
 */
public class Coordinate implements Comparable<Coordinate> {
	private String chromosome;
	private int start;
	private int end;
	public Coordinate(String chrom, int start, int end) {
		chromosome = chrom;
		this.start = start;
		this.end = end;
	}
	public String getChromosome() {
		return chromosome;
	}
	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}
	public int getStart() {
		return start;
	}
	public void setStart(int start) {
		this.start = start;
	}
	public int getEnd() {
		return end;
	}
	public void setEnd(int end) {
		this.end = end;
	}
	@Override
	public int compareTo(Coordinate c2) {
		if (this.start < c2.start)
			return -1;
		else if (start > c2.start)
			return 1;
		return 0;
	}
	
	@Override
	public String toString() {
		// TODO Auto-generated method stub
		return chromosome + "_" + start + "_" + end;
	}
	
	
}
