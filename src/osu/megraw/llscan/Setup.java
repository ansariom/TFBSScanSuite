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
public class Setup {

    public static final int NCHARS = 4;

    /**
     * Returns defined index of characters A, C, G, T (or a, c, g, t).
     *
     * @param chr     Character
     *
     * @return        Defined character index
     *
     */
    public static int charInd(char chr) {

        int ind = -1;

        if (chr == 'A' || chr == 'a') { ind = 0; }
        else if (chr == 'C' || chr == 'c') { ind = 1; }
        else if (chr == 'G' || chr == 'g') { ind = 2; }
        else if (chr == 'T' || chr == 't') { ind = 3; }
        else {}

        return ind;
    }

    /**
     * Checks the validity of a character according to method charInd.
     *
     * @param chr     Character
     *
     * @exception     BadCharException  if input index is unrecognized
     *
     */
    public static void checkChar(char chr) throws BadCharException {
        int ind = charInd(chr);
        if (ind < 0 || ind > (NCHARS - 1)) {
            throw new BadCharException("Unrecognized character " + chr);
        }
    }

    /**
     * Returns defined character (A, C, G, or T) associated with index.
     *
     * @param ind       Character index
     *
     * @param upperCase True if returned character should be uppercase.
     *
     * @return          Defined character
     *
     * @exception BadIndexException  if input index is unrecognized
     *
     */
    public static char indChar(int ind, boolean upperCase) {

        char[] UpperCase = {'A', 'C', 'G', 'T'};
        char[] LowerCase = {'a', 'c', 'g', 't'};

        if (upperCase) {
            return UpperCase[ind];
        }
        else {
            return LowerCase[ind];
        }

    }

    /**
     * Checks the validity of a character index.
     *
     * @param ind     Character index
     *
     * @exception     BadIndexException  if input character index is unrecognized
     *
     */
    public static void checkInd(int ind) throws BadIndexException {
        if (ind < 0 || ind > (NCHARS - 1)) {
            throw new BadIndexException("Unrecognized character index " + ind);
        }
    }

}
