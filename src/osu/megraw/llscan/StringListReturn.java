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
public class StringListReturn {

    public String[] labels;
    public String[][] values;

    public StringListReturn() {
    }

    public StringListReturn(String[] in_labels, String[][] in_values) {
        labels = in_labels;
        values = in_values;
    }
}
