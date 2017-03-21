package osu.megraw.llscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Vector;

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
public class Load {
    public Load() {
        try {
            jbInit();
        } catch (Exception ex) {
            ex.printStackTrace();
        }
    }

    /**
     * Loads a FASTA file and returns String arrays of headers and lines.
     *
     * @param in_file File  Input file to be loaded.
     *
     * @return FastaReturn  Contains String arrays of headers and lines.
     *
     * @throws IOException  if file loading fails.
     *
     */
    public static FastaReturn loadFastaFile(File in_file) throws IOException {

        FastaReturn fa = new FastaReturn();

        BufferedReader inputFile = new BufferedReader(new FileReader(in_file));

        // Read each line from the input file, store headers and lines
        String line;
        int begin;
        Vector linesVec = new Vector();
        Vector headersVec = new Vector();
        StringBuilder fastaSeq = new StringBuilder();
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                // Only add in new sequences after the first FASTA header/sequence 
                // pair has been read
                if (headersVec.size() > 0) {
                    linesVec.addElement(fastaSeq.toString());
                    fastaSeq.setLength(0);
                }
                begin = 1;
                while (Character.isWhitespace(line.charAt(begin))) {
                    begin++;
                }
                headersVec.addElement(line.substring(begin, line.length()));
			} else {
                fastaSeq.append(line);
            }
        }

        // Only add remaining sequence if FASTA header was read in
        if (headersVec.size() > 0) {
            linesVec.addElement(fastaSeq.toString());
        }

        Object[] headersArr = headersVec.toArray();
        fa.headers = new String[headersArr.length];
        for (int i = 0; i < headersArr.length; i++) {
            fa.headers[i] = (String)headersArr[i];
        }

        Object[] linesArr = linesVec.toArray();
        fa.lines = new String[linesArr.length];
        for (int i = 0; i < linesArr.length; i++) {
            fa.lines[i] = (String)linesArr[i];
        }

        inputFile.close();
        if (fa.headers.length == fa.lines.length) {
            return fa;
        }
        else {
            throw new IOException("Number of headers does not match number of lines in file.");
        }
    }

    /**
     * Loads a FASTA file, reverse complements all sequences, and returns String arrays of headers and lines.
     *
     * @param in_file File  Input file to be loaded.
     *
     * @return FastaReturn  Contains String arrays of headers and lines.
     *
     * @throws IOException  if file loading fails.
     *
     */
    public static FastaReturn loadFastaFileRevcomp(File in_file) throws IOException {

        FastaReturn fa = new FastaReturn();

        BufferedReader inputFile = new BufferedReader(new FileReader(in_file));

        // Read each line from the input file, store headers and lines
        String line;
        int begin;
        Vector linesVec = new Vector();
        Vector headersVec = new Vector();
        StringBuilder fastaSeq = new StringBuilder();
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                // Only add in new sequences after the first FASTA header/sequence 
                // pair has been read
                if (headersVec.size() > 0) {
                    linesVec.addElement(fastaSeq.toString());
                    fastaSeq.setLength(0);
                }
                begin = 1;
                while (Character.isWhitespace(line.charAt(begin))) {
                    begin++;
                }
                headersVec.addElement(line.substring(begin, line.length()));
			} else {
                fastaSeq.append(line);
            }
        }

        Object[] headersArr = headersVec.toArray();
        fa.headers = new String[headersArr.length];
        for (int i = 0; i < headersArr.length; i++) {
            fa.headers[i] = (String)headersArr[i];
        }

        Object[] linesArr = linesVec.toArray();
        fa.lines = new String[linesArr.length];
        for (int i = 0; i < linesArr.length; i++) {
            fa.lines[i] = Utils.reverseComplement((String)linesArr[i]);
        }

        inputFile.close();
        if (fa.headers.length == fa.lines.length) {
            return fa;
        }
        else {
            throw new IOException("Number of headers does not match number of lines in file.");
        }
    }

    /**
     * Loads file containing lists of sequences from which to create PWMs
     *
     * @param filename String Filename to be loaded.
     *
     * @return FreqReturn object which contains the list of frequency matrices.
     *
     * @throws IOException if file loading fails.
     *
     */
    public static FreqReturn loadPWMSeqFile(String filename) throws IOException {
        BufferedReader inputFile;

        // Read each set of sequences in the file
        String[] Freq_Headers;
        int[][][] Freqs;
        String line = "";
        int chrInd;

        // First count number of sequence sets
        inputFile = new BufferedReader(new FileReader(filename));
        int nmats = 0;
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                nmats++;
            }
        }
        inputFile.close();

        // Now load all frequency matrices
        Freqs = new int[nmats][][];
        Freq_Headers = new String[nmats];
        inputFile = new BufferedReader(new FileReader(filename));
        nmats = 0;
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                Freq_Headers[nmats] = line;
                Freqs[nmats] = new int[0][Setup.NCHARS];

            }
            else if (!(line.length() == 0)) {
                if (Freqs[nmats].length == 0) {
                    Freqs[nmats] = new int[line.length()][Setup.NCHARS];
                }
                // Tally bases in this sequence
                char[] S = line.toCharArray();
                for (int i = 0; i < S.length; i++) {
                    chrInd = Setup.charInd(S[i]);
                    if (chrInd == -1) { // unrecognized character
                        System.out.println("Unrecognized char " + S[i] + " in sequences of file " + filename + " under header " + Freq_Headers[nmats]);
                    }
                    else {
                        Freqs[nmats][i][chrInd] += 1;
                    }
                }
            }
            else {
                // Done reading this sequence set
                nmats++;
            }
        }
        inputFile.close();
        return (new FreqReturn(Freq_Headers, Freqs));
    }

    /**
     * Loads file containing lists of sequences for loading as character arrays
     *
     * @param filename String Filename to be loaded.
     *
     * @return SeqArrayReturn object which contains the list of character arrays.
     *
     * @throws IOException if file loading fails.
     *
     */
    public static SeqArrayReturn loadSeqArrayFile(String filename) throws IOException {
        BufferedReader inputFile;

        // Read each set of sequences in the file
        String[] Headers;
        String[] Labels;
        char[][][] SeqArrays;
        String line = "";
        String seqlines = "";
        int chrInd;

        // First count number of sequence sets
        inputFile = new BufferedReader(new FileReader(filename));
        int nmats = 0;
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                nmats++;
            }
        }
        inputFile.close();

        // Now load all frequency matrices
        SeqArrays = new char[nmats][][];
        Headers = new String[nmats];
        Labels = new String[nmats];
        inputFile = new BufferedReader(new FileReader(filename));
        nmats = 0;
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                String label = line.substring(1, line.length());
                Labels[nmats] = label;
                Headers[nmats] = line;
                SeqArrays[nmats] = new char[0][0];
                seqlines = "";
            }
            else if (!(line.length() == 0)) {
                seqlines += line;
                seqlines += "#";
            }
            else {
                // Done reading this sequence set, process
                String[] mat = seqlines.split("#");
                SeqArrays[nmats] = new char[mat.length][];
                for (int i = 0; i < mat.length; i++) {
                    SeqArrays[nmats][i] = mat[i].toCharArray();
                }

                nmats++;
            }
        }
        inputFile.close();
        return (new SeqArrayReturn(Labels, Headers, SeqArrays));
    }

    /**
     * Loads a CpG Islands output file and returns a Hashtable of Vectors of coordinate Pairs.
     *
     * @param in_file File  Input file to be loaded.
     *
     * @return Hashtable  Contains a Hashtable of Vectors of coordinate Pairs.
     *
     * @throws IOException  if file loading fails.
     *
     */
    public static Hashtable loadCpGFile(File in_file) throws IOException {

        BufferedReader inputFile = new BufferedReader(new FileReader(in_file));

        // Read and parse each line from the input file
        Hashtable hash = new Hashtable();
        String header = inputFile.readLine();

        String name, line;
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            String[] parts = line.split("\\s+");
            name = parts[0];
            Vector vec = new Vector();
            for (int i = 1; i < parts.length; i++) {
                String[] coordstr = parts[i].split(",");
                Pair pair = new Pair(Integer.parseInt(coordstr[0]), Integer.parseInt(coordstr[1]));
                vec.addElement(pair);
            }
            hash.put(name, vec);
        }

        inputFile.close();
        return hash;
    }

    /**
     * Loads a file containing labeled positional draws of the form
     * sequenceLabel_position, for example T543F3234_-56,
     * in the first whitespace-delimited field
     * and returns a Hashtable of String Vectors by label and position set.
     *
     * @param in_file File  Input file to be loaded.
     *
     * @return Hashtable  Contains a Hashtable of String Vectors.
     *
     * @throws IOException  if file loading fails.
     *
     */
    public static Hashtable loadDrawRepFile(File in_file) throws IOException {

        BufferedReader inputFile = new BufferedReader(new FileReader(in_file));

        // Read and parse each line from the input file
        Hashtable hash = new Hashtable();
        String header = inputFile.readLine();

        String name, line, label, position;
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            String[] parts = line.split("\\s+");
            name = parts[0];
            String[] nameparts = name.split("_");
            label = nameparts[0];
            position = nameparts[1];
            if (!hash.containsKey(label)) {
                hash.put(label, new Vector());
            }
            Vector vec = (Vector) hash.get(label);
            vec.addElement(position);
            hash.put(label, vec);
        }

        inputFile.close();
        return hash;
    }

    public static DoubleColReturn loadDoubleColumnFile(File in_file) throws IOException {
        return loadDoubleColumnFile(in_file, null);
    }

    /**
     * Loads a file containing a column of label strings followed by a column of doubles.
     *
     * @param in_file File  Input file to be loaded.
     *
     * @return DoubleColReturn  Contains String arrays of headers and lines.
     *
     * @throws IOException  if file loading fails.
     *
     */
    public static DoubleColReturn loadDoubleColumnFile(File in_file, HashSet<String> allowedNames) throws IOException {

        DoubleColReturn dr = new DoubleColReturn();

        BufferedReader inputFile = new BufferedReader(new FileReader(in_file));

        String line;
        Vector valuesVec = new Vector();
        Vector labelsVec = new Vector();
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            String[] parts = line.split("\\s+");

            if (allowedNames != null && !allowedNames.contains(parts[0])) {
                continue;
            }

            labelsVec.addElement(parts[0]);
            valuesVec.addElement(parts[1]);
        }

        Object[] labelsArr = labelsVec.toArray();
        dr.labels = new String[labelsArr.length];
        for (int i = 0; i < labelsArr.length; i++) {
            dr.labels[i] = (String)labelsArr[i];
        }

        Object[] valuesArr = valuesVec.toArray();
        dr.values = new double[valuesArr.length];
        for (int i = 0; i < valuesArr.length; i++) {
            dr.values[i] = Double.parseDouble((String)valuesArr[i]);
        }

        inputFile.close();
        return dr;
    }

    /**
     * Loads a file containing a column of label strings followed by a column of strings.
     *
     * @param in_file File  Input file to be loaded.
     *
     * @return StringListReturn  Contains String arrays of headers and lines.
     *
     * @throws IOException  if file loading fails.
     *
     */
    public static StringListReturn loadStringListColumnFile(File in_file, char listDelimiterChar) throws IOException {

        StringListReturn slr = new StringListReturn();

        BufferedReader inputFile = new BufferedReader(new FileReader(in_file));

        String line;
        Vector valuesVec = new Vector();
        Vector labelsVec = new Vector();
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            String[] parts = line.split("\\t");
            labelsVec.addElement(parts[0]);
            valuesVec.addElement(parts[1]);
        }

        Object[] labelsArr = labelsVec.toArray();
        slr.labels = new String[labelsArr.length];
        for (int i = 0; i < labelsArr.length; i++) {
            slr.labels[i] = (String)labelsArr[i];
        }

        Object[] valuesArr = valuesVec.toArray();
        slr.values = new String[valuesArr.length][];
        for (int i = 0; i < valuesArr.length; i++) {
            String itemList = (String) valuesArr[i];
            if (itemList.charAt(0) == listDelimiterChar) {
                itemList = itemList.substring(1, itemList.length());
            }
            if (itemList.charAt(itemList.length() - 1) == listDelimiterChar) {
                itemList = itemList.substring(0, itemList.length() - 1);
            }
            String[] items = itemList.split("\\" + String.valueOf(listDelimiterChar));
            slr.values[i] = new String[items.length];
            for (int j = 0; j < items.length; j++) {
                slr.values[i][j] = items[j];
            }
        }

        inputFile.close();
        return slr;
    }

    /**
     * Loads a file containing a column of label strings followed by a matrix of doubles.
     *
     * @param in_file File  Input file to be loaded.
     *
     * @return DoubleColReturn  Contains String arrays of headers and lines.
     *
     * @throws IOException  if file loading fails.
     *
     */
    public static DoubleMatReturn loadDoubleMatrixFile(File in_file, boolean hasHeader) throws IOException {
        return loadDoubleMatrixFile(in_file, hasHeader, null);
    }

    public static DoubleMatReturn loadDoubleMatrixFile(File in_file, boolean hasHeader, HashSet<String> includedFeatures) throws IOException {

        DoubleMatReturn mr = new DoubleMatReturn();

        BufferedReader inputFile = new BufferedReader(new FileReader(in_file));

        String line;
        String[] headerArr;
        Vector linesVec = new Vector();

        if (hasHeader) {
            line = inputFile.readLine();
            mr.colLabels = line.split("\\s+");
        }

        while ((line = inputFile.readLine()) != null) {
            String[] parts = line.split("\\s+");
            if (parts.length > 0 && includedFeatures != null && !includedFeatures.contains(parts[0])) {
                continue;
            }
            linesVec.addElement(line.trim());
        }

        int nLines = linesVec.size();
        mr.rowLabels = new String[nLines];
        mr.values = new double[nLines][];
        mr.flags = new boolean[nLines]; // flag rows which cannot be parsed as doubles
        mr.nflags = 0; // count number of flagged rows
        for (int i = 0; i < nLines; i++) {
            line = (String)linesVec.get(i);
            String[] parts = line.split("\\s+");
            mr.rowLabels[i] = parts[0];
            mr.values[i] = new double[parts.length - 1];
            for (int j = 1; j < parts.length; j++) {
                try {
                    mr.values[i][j - 1] = Double.parseDouble(parts[j]);
                }
                catch (java.lang.NumberFormatException nfe) {
                    mr.flags[i] = true;
                }
            }
            if (mr.flags[i] == true) { mr.nflags++; }
        }

        inputFile.close();
        return mr;
    }

    /**
     * Loads file containing a list of PWMs, retrieves a specific PWM labeled
     * by PWM_Label
     *
     * @param filename String Filename to be loaded.
     *
     * @param PWM_Label String Label of PWM to be retrieved.
     *
     * @return double[][] PWM associated with PWM_Label.
     *
     * @return PseudoCountsVal Value to be added to each entry of loaded PWM counts before normalizing.
     *
     * @throws IOException if file loading fails.
     *
     */
    public static double[][] loadPWM(String filename, String PWM_Label, double PseudoCountsVal) throws IOException {
        BufferedReader inputFile = new BufferedReader(new FileReader(filename));

        // Find desired pwm in the input file
        String PWM = PWM_Label;
        String line = "";
        String header = "";
        String matlines = "";
        double[][] Theta = new double[0][0];
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                int begin = 1;
                while (!line.substring(begin, begin + 3).equals("AC="))
                {
                    begin++;
                }
                header = line.substring(begin + 3, begin + 9);

                matlines = "";
            }
            else if (line.startsWith("=")) {
                matlines += line;
                matlines += "#";
            }
            else {
                // Done reading this PWM, process if it's the desired one
                if (header.equals(PWM)) {
                    String[] mat = matlines.split("#");
                    Theta = new double[mat.length][];
                    for (int i = 0; i < mat.length; i++) {
                        String[] values = mat[i].split("\\s+");
                        Theta[i] = new double[values.length - 1];
                        for (int j = 1; j < values.length; j++) {
                            Theta[i][j-1] = Double.parseDouble(values[j]);
                        }
                    }

                    // Add pseudocounts: add PseudoCountsVal to each entry of Theta
                    for (int mp = 0; mp < Theta.length; mp++) {
                        for (int k = 0; k < Theta[mp].length; k++) {
                            Theta[mp][k] += PseudoCountsVal;
                        }
                    }

                    // Divide by row sums
                    double rowsum;
                    for (int mp = 0; mp < Theta.length; mp++) {
                        rowsum = 0.0;
                        for (int k = 0; k < Theta[mp].length; k++) {
                            rowsum += Theta[mp][k];
                        }
                        for (int k = 0; k < Theta[mp].length; k++) {
                            Theta[mp][k] /= rowsum;
                        }
                    }
                    break;
                }
            }
        }
        inputFile.close();
        return Theta;
    }

    /**
     * Loads file containing a list of PWMs.
     *
     * @param filename String Filename to be loaded.
     *
     * @param PWM_Label String Label of PWM to be retrieved.
     *
     * @return PseudoCountsVal Value to be added to each entry of loaded PWM counts before normalizing.
     *
     * @return double[][] PWM associated with PWM_Label.
     *
     * @throws IOException if file loading fails.
     *
     */
    public static PWMReturn loadPWMFile(String filename, double PseudoCountsVal) throws IOException {
        BufferedReader inputFile;

        // Read each PWM in the file
        String[] PWM_Labels;
        String[] PWM_Headers;
        double[][][] PWMs;
        String line = "";
        String matlines = "";

        // First count number of PWMs
        inputFile = new BufferedReader(new FileReader(filename));
        int nmats = 0;
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                nmats++;
            }
        }
        inputFile.close();

        // Now load all PWMs
        PWMs = new double[nmats][][];
        PWM_Labels = new String[nmats];
        PWM_Headers = new String[nmats];
        inputFile = new BufferedReader(new FileReader(filename));
        nmats = 0;
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                int begin = 1;
                while (!line.substring(begin, begin + 3).equals("AC="))
                {
                    begin++;
                }
                PWM_Labels[nmats] = line.substring(begin + 3, begin + 9);
                PWM_Headers[nmats] = line;
                PWMs[nmats] = new double[0][0];

                matlines = "";
            }
            else if (line.startsWith("=")) {
                matlines += line;
                matlines += "#";
            }
            else {
                // Done reading this PWM, process
                String[] mat = matlines.split("#");
                PWMs[nmats] = new double[mat.length][];
                for (int i = 0; i < mat.length; i++) {
                    String[] values = mat[i].split("\\s+");
                    PWMs[nmats][i] = new double[values.length - 1];
                    for (int j = 1; j < values.length; j++) {
                        PWMs[nmats][i][j-1] = Double.parseDouble(values[j]);
                    }
                }

                // Add pseudocounts: add PseudoCountsVal to each entry of Theta
                for (int mp = 0; mp < PWMs[nmats].length; mp++) {
                    for (int k = 0; k < PWMs[nmats][mp].length; k++) {
                        PWMs[nmats][mp][k] += PseudoCountsVal;
                    }
                }

                // Divide by row sums
                double rowsum;
                for (int mp = 0; mp < PWMs[nmats].length; mp++) {
                    rowsum = 0.0;
                    for (int k = 0; k < PWMs[nmats][mp].length; k++) {
                        rowsum += PWMs[nmats][mp][k];
                    }
                    for (int k = 0; k < PWMs[nmats][mp].length; k++) {
                        PWMs[nmats][mp][k] /= rowsum;
                    }
                }

                nmats++;
            }
        }
        inputFile.close();
        return (new PWMReturn(PWM_Labels, PWM_Headers, PWMs));
    }

    /**
     * Loads file containing a list of PWMs.
     *
     * @param filename String Filename to be loaded.
     *
     * @param PWM_Label String Label of PWM to be retrieved.
     *
     * @return PseudoCountsVal Value to be added to each entry of loaded PWM counts before normalizing.
     *
     * @return double[][] PWM associated with PWM_Label.
     *
     * @throws IOException if file loading fails.
     *
     */
    public static PWMReturn loadPWMFileSimpleHeader(String filename, double PseudoCountsVal) throws IOException {
        return loadPWMFileSimpleHeader(filename, PseudoCountsVal, null);
    }

    public static PWMReturn loadPWMFileSimpleHeader(String filename, double PseudoCountsVal, HashSet<String> includedFeatures) throws IOException {
        BufferedReader inputFile;

        // Read each PWM in the file
        String[] PWM_Labels;
        String[] PWM_Headers;
        double[][][] PWMs;
        String line = "";
        String matlines = "";

        // First count number of PWMs
        inputFile = new BufferedReader(new FileReader(filename));
        int nmats = 0;
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                String label = line.substring(2, line.length());
                if (includedFeatures != null && !includedFeatures.contains(label)) {
                    continue;
                }

                nmats++;
            }
        }
        inputFile.close();

        // Now load all PWMs
        PWMs = new double[nmats][][];
        PWM_Labels = new String[nmats];
        PWM_Headers = new String[nmats];
        inputFile = new BufferedReader(new FileReader(filename));
        nmats = 0;
        boolean includeThisPWM = true;
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (line.startsWith(">")) {
                String label = line.substring(2, line.length());
                if (includedFeatures != null && !includedFeatures.contains(label)) {
                    includeThisPWM = false;
                    matlines = "";
                    continue;
                }

                PWM_Labels[nmats] = label;
                PWM_Headers[nmats] = line;
                PWMs[nmats] = new double[0][0];

                matlines = "";
            }
            else if (line.startsWith("=")) {
                matlines += line;
                matlines += "#";
            }
            else {
                if (includeThisPWM) {
                    // Done reading this PWM, process
                    String[] mat = matlines.split("#");
                    PWMs[nmats] = new double[mat.length][];
                    for (int i = 0; i < mat.length; i++) {
                        String[] values = mat[i].split("\\s+");
                        PWMs[nmats][i] = new double[values.length - 1];
                        for (int j = 1; j < values.length; j++) {
                            PWMs[nmats][i][j-1] = Double.parseDouble(values[j]);
                        }
                    }

                    // Add pseudocounts: add PseudoCountsVal to each entry of Theta
                    for (int mp = 0; mp < PWMs[nmats].length; mp++) {
                        for (int k = 0; k < PWMs[nmats][mp].length; k++) {
                            PWMs[nmats][mp][k] += PseudoCountsVal;
                        }
                    }

                    // Divide by row sums
                    double rowsum;
                    for (int mp = 0; mp < PWMs[nmats].length; mp++) {
                        rowsum = 0.0;
                        for (int k = 0; k < PWMs[nmats][mp].length; k++) {
                            rowsum += PWMs[nmats][mp][k];
                        }
                        for (int k = 0; k < PWMs[nmats][mp].length; k++) {
                            PWMs[nmats][mp][k] /= rowsum;
                        }
                    }

                    nmats++;
                }

                includeThisPWM = true;
            }
        }
        inputFile.close();
        return (new PWMReturn(PWM_Labels, PWM_Headers, PWMs));
    }

    /**
     * Loads a file containing a column of doubles.
     *
     * @param in_file File  Input file to be loaded.
     *
     * @return double[]  Array of doubles.
     *
     * @throws IOException  if file loading fails.
     *
     */
    public static double[] loadVector(File in_file) throws IOException {

        double[] dr;

        BufferedReader inputFile = new BufferedReader(new FileReader(in_file));

        String line;
        Vector valuesVec = new Vector();
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (!line.equals("")) {
                valuesVec.addElement(line);
            }
        }

        Object[] valuesArr = valuesVec.toArray();
        dr = new double[valuesArr.length];
        for (int i = 0; i < valuesArr.length; i++) {
            dr[i] = Double.parseDouble((String)valuesArr[i]);
        }

        inputFile.close();
        return dr;
    }

    /**
     * Loads a file containing a matrix of doubles.
     *
     * @param in_file File  Input file to be loaded.
     *
     * @return double[][] Matrix of doubles.
     *
     * @throws IOException  if file loading fails.
     *
     */
    public static double[][] loadMatrix(File in_file) throws IOException {

        double[][] dr;

        BufferedReader inputFile = new BufferedReader(new FileReader(in_file));

        String line;
        String[] parts;
        Vector valuesVec = new Vector();
        while ((line = inputFile.readLine()) != null) {
            line.trim();
            if (!line.equals("")) {
                valuesVec.addElement(line);
            }
        }

        dr = new double[valuesVec.size()][];
        for (int i = 0; i < valuesVec.size(); i++) {
            line = valuesVec.get(i).toString();
            parts = line.split("\\s+");
            dr[i] = new double[parts.length];
            for (int j = 0; j < parts.length; j++) {
                dr[i][j] = Double.parseDouble(parts[j]);
            }
        }

        inputFile.close();
        return dr;
    }

    private void jbInit() throws Exception {
    }


}
