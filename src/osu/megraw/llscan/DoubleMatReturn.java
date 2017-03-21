package osu.megraw.llscan;

/**
 * <p>Title: Log Likelihood Scanner Project</p>
 *
 * <p>Description: Log Likelihood Scanner</p>
 *
 * <p>Copyright: Copyright (c) 2005</p>
 *
 * <p>Company: </p>
 *
 * @author Molly Megraw
 * @version 1.0
 */
public class DoubleMatReturn {

    public String[] rowLabels;
    public String[] colLabels;
    public double[][] values;
    public boolean[] flags;
    public int nflags;

    public DoubleMatReturn() {
    }

    public DoubleMatReturn(String[] in_rowLabels, String[] in_colLabels, double[][] in_values) {
        rowLabels = in_rowLabels;
        colLabels = in_colLabels;
        values = in_values;
    }

    public DoubleMatReturn(String[] in_rowLabels, double[][] in_values) {
        rowLabels = in_rowLabels;
        values = in_values;
    }

}
