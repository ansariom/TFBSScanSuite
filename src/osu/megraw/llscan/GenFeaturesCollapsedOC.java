package osu.megraw.llscan;

import java.io.BufferedOutputStream;
import java.io.DataOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Hashtable;
import java.util.List;
import java.util.Random;
import java.util.Vector;
import java.util.concurrent.Callable;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

// Use Apache Commons CLI Parser
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

/**
 * <p>Title: Log Likelihood Scanner Project</p>
 *
 * <p>Description: Log Likelihood Scanner</p>
 *
 * <p>Copyright: Copyright (c) 2005-2013</p>
 *
 * <p>Company: </p>
 *
 * @author Molly Megraw
 * @version 1.0
 */
public class GenFeaturesCollapsedOC {

    public static char[][] S; // holds input sequences from input FASTA file
    public static DoubleColReturn Scores;
    public static int seqLength = -1;

    // CLI optional arguments
    public static boolean help = false;
    public static String posArrString = "0";
    public static int nucsDownStream = -1;
    public static String scoreCutoffs_Fname = null;
    public static int BG_WIN = 250;
    public static int nproc = 1;
    public static int flankingFeatureWindowLength = 0;
    public static int numFlankingWindows = 0;
    public static String processPeak = null;
    public static double PseudoCountsVal = 0.01; // SHAWN: was 0.25
    public static boolean outputBin = false;
    public static boolean byNT = false;
    public static String map_Fname = null;

    // CLI required arguments
    public static String factorInfo_FWD_Fname;
    public static String factorInfo_REV_Fname;
    public static String promoterSeqs_Fname;
    public static String pwms_Fname;
    public static String out_Fname;

    public static boolean USE_MARKOV_1_BG = true;

    public static double[] bgdefault = {0.25, 0.25, 0.25, 0.25};
    public static double[][] bgM1default = {{0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25}};

    /**
     * This code has chenged by Mitra in order to print out "true" genome coordinates of ROE sites for genes located on reverse strand
     * into map_file
     * We use Map file for OC-overlap computation!
     */
    public GenFeaturesCollapsedOC() {
    }

    public static void main(String[] args) throws java.io.IOException, java.lang.InterruptedException {
        GenFeaturesCollapsedOC ss = new GenFeaturesCollapsedOC();

        parseArgs(args);

        // Read promoter sequences file
        FastaReturn Seqs;
        Seqs = Load.loadFastaFile(new File(promoterSeqs_Fname));
        S = new char[Seqs.lines.length][];
        for (int i = 0; i < Seqs.lines.length; i++) {
            S[i] = Seqs.lines[i].toCharArray();

            if (seqLength == -1) seqLength = S[i].length;

            if (seqLength != S[i].length) {
                System.err.println("Error: input FASTA sequences of different length found with (" + Seqs.headers[i] + ")");
                return;
            }
        }

        // no TSS offset given - autoset it to midpoint of input sequence length
        if (nucsDownStream == -1) {
            // For odd length sequences, half the length of the sequence rounded down will be midpoint since
            // sequences are 0-indexed (i.e. a sequence 10001 nt long needs nucsAfterTSS=5000 to set the TSS
            // at the midpoint, i.e. character 5001 in a 1-index based setting)
            nucsDownStream = seqLength / 2;  // integer arithmetic - automatically floors
        }

        // Get labels for seqs
        HashMap <String, RefSeqData> fastaCoords = new HashMap <String, RefSeqData>();
        String[] seqLabels = new String[Seqs.headers.length];
        for (int i = 0; i < seqLabels.length; i++) {
            String header = Seqs.headers[i].substring(0, Seqs.headers[i].length());
            String[] header_parts = header.split("\\s+");
            seqLabels[i] = header_parts[0];

            // user wants to auto-produce a mapping file for features generated and the genomic coordinates of each feature
            // In order to use this feature, each sequence *must* have the header id formatted as follows
            // >fastaID_<refseq>_<start> where <refseq> is an alphanumeric value ID of the reference sequence
            // this fasta sequence came from and <start> is an integer denoting the start coordinate of this fasta sequence
            // within the reference sequence
            if (map_Fname != null) {
                String[] genomic_parts = header_parts[0].split("_");
                if (genomic_parts.length < 3) { // need *at least* FASTAID, REFSEQ, START
                    System.err.println("FASTA identifier does not have the required format to continue.  Fix header, or restart without making a mapping file.");
                    return;
                }
                try {
                    String ref = genomic_parts[1];
                    int start = Integer.parseInt(genomic_parts[2]);
                    String strand = genomic_parts[3];
                    RefSeqData refData = new RefSeqData(ref, start, strand);
                    fastaCoords.put(header_parts[0], refData);
                }
                catch (NumberFormatException e) {
                    System.out.println("Unexpected formatting problem: " + e.getMessage());
                    e.printStackTrace();
                    return;
                }
            }
        }


        // Read PWM file, store rev comp PWMs for scanning opposite strand
        // SHAWN: THIS WAS NOT HERE IN JASONS (this file), BUT WAS IN TAJ's
        // PseudoCountsVal = 0.0;
        PWMReturn pwms = Load.loadPWMFileSimpleHeader(pwms_Fname, PseudoCountsVal);
        pwms.ComputeRevCmpPWMs();

        // Read score threshold file corresponding to PWM file
        if (scoreCutoffs_Fname != null) {
            Scores = Load.loadDoubleColumnFile(new File(scoreCutoffs_Fname));
        } else {
            // No threshold was given - assume that threshold will be zero
            Scores = new DoubleColReturn(pwms.labels);
        }

        // Read factor information file corresponding to PWM file
        DoubleMatReturn facInfo_FWD = Load.loadDoubleMatrixFile(new File(factorInfo_FWD_Fname), true);
        DoubleMatReturn facInfo_REV = Load.loadDoubleMatrixFile(new File(factorInfo_REV_Fname), true);

        // Process factorInfo file into variable window definitions
        int NPI = 5; // Number of windows used to cover the peak region
        int N = NPI + 2; // total number of windows

        DoubleMatReturn[] facInfoArr = {facInfo_FWD, facInfo_REV};
        String[] strandArr = {"FWD", "REV"};

        int nWinVars = 0;
        for (int fi = 0; fi < facInfoArr.length; fi++) {
            DoubleMatReturn facInfo = facInfoArr[fi];
            nWinVars += N*(facInfo.rowLabels.length - facInfo.nflags); // N * number of factors with numeric parameters
        }
        String[] winVarNames = new String[nWinVars];
        int[] winVarL = new int[nWinVars];
        int[] winVarR = new int[nWinVars];
        int[] winVarPWMIndex = new int[nWinVars];
        String[] winVarStrand = new String[nWinVars];
        int Mcol = 0;  // facInfo data column for window center
        int Hcol = 1;  // facInfo data column for window half-width

        // Determine positions of feature windows relative to TSS using the ROE position information
        int v = 0;
        for (int fi = 0; fi < facInfoArr.length; fi++) {
            DoubleMatReturn facInfo = facInfoArr[fi];
            for (int p = 0; p < facInfo.values.length; p++) {
                if (facInfo.flags[p] == true) { continue; } // contained non-numeric data
                String varNameBase = facInfo.rowLabels[p];
                int cent = (int) facInfo.values[p][Mcol];
                double halfWidth = facInfo.values[p][Hcol];
                double entireWin = 2.0 * halfWidth;
                int Half_WIN = (int) (entireWin/(((double)NPI) + 1));
                if (Half_WIN <= 0) { Half_WIN = 1; }
                int hwv = (int) (((double) (N - 1))/2.0);
                for (int wv = -hwv; wv <= hwv; wv++) {
                    int wvn = wv + hwv + 1;
                    winVarNames[v] = varNameBase + "_" + strandArr[fi] + "_" + wvn;
                    int pcent, HW;
                    if (wv == -hwv && halfWidth < 100.0 ) { // left peak flank for narrow peaks
                        pcent = cent - (int) entireWin;
                         int lnext = cent + wv * Half_WIN;
                         HW = lnext - pcent;
                    }
                    else if (wv == hwv && halfWidth < 100.0 ) { // right peak flank for narrow peaks
                        pcent = cent + (int) entireWin;
                        int rnext = cent + wv * Half_WIN;
                        HW = pcent - rnext;
                    }
                    else { // peak
                        pcent = cent + wv * Half_WIN;
                        HW = Half_WIN;
                    }
                    winVarL[v] = pcent - HW;
                    winVarR[v] = pcent + HW;
                    winVarPWMIndex[v] = p;
                    winVarStrand[v] = strandArr[fi];
                    v++;
                }
            }
        }

        // Determine positions from which to calculate features (based on a random position pull,
        // range position pull, replicate of previous random pull, or fixed number of positions)
        String[] posArgs = posArrString.split("\\s+");
        Random rand = new Random();
        boolean drawRand = false;
        boolean drawRep = false;
        int drawRandN = 0;
        int drawRandL = 0;
        int drawRandR = 0;
        int rangeS = 0;
        int rangeL = 0;
        int rangeR = 0;
        int[] fixedPosArr = {};
        Hashtable repHash = new Hashtable();
        if (posArgs[0].equals("Rand")) {
            drawRand = true;
            drawRandN = Integer.parseInt(posArgs[1]);
            drawRandL = Integer.parseInt(posArgs[2]);
            drawRandR = Integer.parseInt(posArgs[3]);
        }
        else if (posArgs[0].equals("Range")) {
            rangeS = Integer.parseInt(posArgs[1]);
            rangeL = Integer.parseInt(posArgs[2]);
            rangeR = Integer.parseInt(posArgs[3]);
            int npts = (rangeR - rangeL + 1)/rangeS;
            fixedPosArr = new int[npts];
            for (int i = 0; i < npts; i += rangeS) {
                fixedPosArr[i] = rangeL + i;
            }
        }
        else if (posArgs[0].equals("Replicate")) {
            drawRep = true;
            String rep_Fname = posArgs[1];
            repHash = Load.loadDrawRepFile(new File(rep_Fname));
        }
        else {
            fixedPosArr = new int[posArgs.length];
            for (int i = 0 ; i < posArgs.length; i++) {
                fixedPosArr[i] = Integer.parseInt(posArgs[i]);
            }
        }

        // Process each sequence
        int[] posArr = {};
        if (drawRand) {
            posArr = new int[drawRandN];
        }
        else if (!drawRep) {
            posArr = fixedPosArr;
        }
        else {}

        ThreadPoolExecutor threadPool = new ThreadPoolExecutor(nproc, nproc, 4, TimeUnit.SECONDS, new LinkedBlockingQueue());
        List<Future<TagScoreResult>> results = new ArrayList<Future<TagScoreResult>>();
        ArrayList<String> futureSeqNames = new ArrayList<String>();

        // file to output mapped genomic coordinates to
        PrintWriter mapOutFile = null;
        if (map_Fname != null) {
            mapOutFile = new PrintWriter(new FileWriter(map_Fname));

            // header for mapping file
            mapOutFile.println("seq_id\tfeature_id\trelative_start\trelative_end\ttss_start\tabsolute_start\tabsolute_end\tref_seq\tgenomic_start\tgenomic_end");
        }
        for (int i = 0; i < S.length; i++) {
        	String seqName = seqLabels[i];
            if (processPeak != null && !seqName.equals(processPeak)) { // we are not responsible for this peak
                continue;
            }

            // Set up random drawing from positions around the sequence
			if (drawRand) {
				posArr = new int[drawRandN];
			    Hashtable hash = new Hashtable();
			    for (int d = 0; d < drawRandN; d++) {
			        // Draw position at random without replacement from range
			        int draw = drawRandL + rand.nextInt(drawRandR - drawRandL + 1);
			        while (hash.containsKey(new Integer(draw))) {
			            draw = drawRandL + rand.nextInt(drawRandR - drawRandL + 1);
			        }
			        posArr[d] = draw;
			        hash.put(new Integer(draw), new Integer(1));
			    }
			}

            // Set up a replicate of a previous random drawing
			if (drawRep) {
			    if (!repHash.containsKey(seqName)) {
			        System.out.println("Can't find positions for " + seqName);
			        System.exit(1);
			    }
			    Vector posVec = (Vector) repHash.get(seqName);
			    int vecSize = posVec.size();
			    posArr = new int[vecSize];
			    for (int p = 0; p < vecSize; p++) {
			        posArr[p] = Integer.parseInt((String) posVec.elementAt(p));
			    }
			}

            // If neither drawRand or drawRep was used, this will only scan the positions set up in posArr
            // from the command line arguments which can be one more or positions, based on how the user
            // set their scans up ( p = 0 as default, means the ref point is tss)
			for (int p = 0; p < posArr.length; p++) {
				TagScoresRunner run = new TagScoresRunner(S[i], Scores, nucsDownStream, BG_WIN,
					seqName, pwms, nWinVars, winVarL, winVarR,
					winVarPWMIndex, winVarStrand, rand, drawRand, drawRep,
					drawRandN, drawRandL, drawRandR, repHash,
					posArr, i, p, flankingFeatureWindowLength, numFlankingWindows);

                try {
    				results.add(threadPool.submit(run));
                    futureSeqNames.add(seqName);
                }
                catch (Exception ex) {
                	ex.printStackTrace();
                    System.err.println("Error processing example "+seqName+": "+ex.getMessage()+". Quitting.");
                    threadPool.shutdownNow();
                    while(!threadPool.isTerminated()) {
                        threadPool.awaitTermination(2, TimeUnit.SECONDS);
                    }
                    System.exit(1);
                }

                /** 
                 * @author mitra
                 * This Map file is generated by Jason in order to spit out the genomic coordinates of
                 * ROE windows then we use it for computation of OC-Overlap with respect to these coordinates.
                 * 
                 * Note : This part was computing wrong coordinates for genes on anti-sense strand!! It is been 
                 * corrected here.
                 */
                if (map_Fname != null) {
                    RefSeqData seqInfo = fastaCoords.get(seqName);
                    // Shouldn't be missing from the hashmap that we set up above...error introduced in the code most likely 
                    
                    if (seqInfo == null) {
                        System.out.println(seqName);
                        return;
                    }
                    int tssLoc = S[i].length - nucsDownStream; // for + strand 
                    int winL = 0;
                    int relLocL = 0;
                    int relLocR = 0;

                    for(int wv = 0; wv < nWinVars; wv++) {
                        int winVarPWMIdx = winVarPWMIndex[wv];
                        String pwmIdx = pwms.labels[winVarPWMIdx];
                        int w = pwms.pwms[winVarPWMIdx].length;
                        
                		if (wv % 7 == 0) {
                			winL = winVarL[wv] + posArr[p] + tssLoc; 
                			if (winL < 0) { winL = 0; }
                			relLocL = winVarL[wv];
                		} else if (wv % 7 == 6) {
                			int winR = winVarR[wv] + posArr[p] + tssLoc;
                			if (winR > S[i].length - w) { winR = S[i].length - w; }
                			int genomic_start = seqInfo.start + winL;
                            int genomic_end = seqInfo.start + winR;
                            
                            relLocR = winVarR[wv];

                            if (seqInfo.strand.equalsIgnoreCase("-")) {
                            	int upstreamLen = S[i].length - nucsDownStream - 1;
                            	tssLoc = S[i].length - upstreamLen;
                            	
                            	winL = (-1 * relLocR) + posArr[p] + tssLoc; 
                            	winR = (-1 * relLocL) + posArr[p] + tssLoc; 
                            	if (winL < 0) { winL = 0; }
                            	if (winR > S[i].length - w) { winR = S[i].length - w; }
                            	genomic_start = seqInfo.start + winL;
                                genomic_end = seqInfo.start + winR;
                            }
                            mapOutFile.println(seqName + "_" + posArr[p] + "\t" + winVarNames[wv] + "\t" + relLocL + "\t" + relLocR + "\t" + 
                                    tssLoc + "\t" + winL + "\t" + winR + "\t" + seqInfo.id + "\t" + genomic_start + "\t" + genomic_end);
                            relLocL = relLocR = winL = winR = 0;                        }
                    }
                }
			}
        }
        if (map_Fname != null) mapOutFile.close();
      
        // Values calculated and needed for outputing results in Taj's binary format
        long numberNonzeros = 0; // How many features generated non-zero values
        int numFeatures = -1; // Number of columns printed out - queried from the results array length found below...
        int totalLabelLength = 0; // Byte length of the final header string
        int totalNumRows = results.size();  // Number of sequences that features were generated for

        // Specify charset formatting for binary output - can change it up later if you want here - don't know how it
        // will affect downstream applications though - check this out to be sure before you change it.
        Charset charset = Charset.forName("US-ASCII"); 

        String[] seqContentNames = {"GCcontent", "CAcontent", "GAcontent"};
        if (outputBin) {
            // Calculate header byte length for binary output
            for (int wv = 0; wv < nWinVars; wv++) {
                if (flankingFeatureWindowLength == 0) {
                    if (!byNT) {
                        String colName = winVarNames[wv];
                        byte[] columnName = Arrays.copyOf(colName.getBytes(charset), Math.min(255, colName.length()));
                        totalLabelLength = columnName.length + 1; // +1 for the byte written to record string length in final output
                    } else {
                        for (int j=winVarL[wv]; j <= winVarR[wv]; j++) {
                            String colName = winVarNames[wv] + "_" + j;
                            byte[] columnName = Arrays.copyOf(colName.getBytes(charset), Math.min(255, colName.length()));
                            totalLabelLength = columnName.length + 1; // +1 for the byte written to record string length in final output
                        }
                    }
                }
                else {
                    for (int i = -numFlankingWindows; i <= numFlankingWindows; i++) {
                        String colName = winVarNames[wv] + "_" + i;
                        byte[] columnName = Arrays.copyOf(colName.getBytes(charset), Math.min(255, colName.length()));
                        totalLabelLength = columnName.length + 1; // +1 for the byte written to record string length in final output
                    }
                }
            }

            // If we are reporting features by individual nucleotides - assume scanning only occurs within ROEs and no sequence features
            // (i.e. GC/CA/GA content) are generated
            if (!byNT) {
                if (flankingFeatureWindowLength == 0) {
                    for (int j = 0; j < seqContentNames.length; j++) {
                        String colName = seqContentNames[j];
                        byte[] columnName = Arrays.copyOf(colName.getBytes(charset), Math.min(255, colName.length()));
                        totalLabelLength = columnName.length + 1; // +1 for the byte written to record string length in final output
        
                    }
                }
                else {
                    for (int i = -numFlankingWindows; i <= numFlankingWindows; i++) {
                        for (int j = 0; j < seqContentNames.length; j++) {
                            String colName = seqContentNames[j] + "_" + i;
                            byte[] columnName = Arrays.copyOf(colName.getBytes(charset), Math.min(255, colName.length()));
                            totalLabelLength = columnName.length + 1; // +1 for the byte written to record string length in final output
                        }
                    }
                }
            }
        }

        // Get all the results!
        ArrayList <TagScoreResult> finalResults = new ArrayList <TagScoreResult>();
        int currentSeq = 0; // If error is thrown, keep track of sequence we were on so error prints out
                            // the actual sequence being worked on
        for (int j = 0; j < results.size(); j++) {
            Future<TagScoreResult> result = results.get(j);

            TagScoreResult res = null;
            try {
                res = result.get(); // get() blocks until the result is available

                finalResults.add(res);

                numberNonzeros += res.numberNonzeros;

                // Set the numFeatures equal to number of values obtained from results - all scans should have the
                // *same* number of features generated, so only need to set this once
                if (numFeatures == -1) numFeatures = res.results.length;

                currentSeq += 1;
            }
            catch (Exception ex) {
            	ex.printStackTrace();
                String seqName = futureSeqNames.get(currentSeq);
                System.err.println("Error processing example "+seqName+": "+ex.getMessage()+". Quitting.");
                threadPool.shutdownNow();
                while(!threadPool.isTerminated()) {
                    threadPool.awaitTermination(2, TimeUnit.SECONDS);
                }
                System.exit(1);
            }
        }

        PrintWriter outFileVars = null;
        DataOutputStream outFileVarsBin = null;

        // Open file for printing variables
        if (!outputBin) {
            outFileVars = new PrintWriter(new FileWriter(out_Fname));
        } else {
            outFileVarsBin = new DataOutputStream(new BufferedOutputStream(new FileOutputStream(out_Fname)));
        }


        // Print binary header size values (e.g., num. of columns, rows, etc.)
        if (outputBin) {
            outFileVarsBin.writeByte('!');
            outFileVarsBin.writeByte('!');
            outFileVarsBin.writeInt(totalNumRows);
            outFileVarsBin.writeInt(numFeatures);
            outFileVarsBin.writeLong(numberNonzeros);
            outFileVarsBin.writeInt(totalLabelLength);
        }

        // Print header
        for (int wv = 0; wv < nWinVars; wv++) {
            if (flankingFeatureWindowLength == 0) {
                if (!byNT) {
                    if (!outputBin) {
                        outFileVars.print("\t" + winVarNames[wv]);
                    } else {
                        writeStringToFile(outFileVarsBin, winVarNames[wv], charset);
                    }
                } else {
                    if (!outputBin) {
                        for (int j=winVarL[wv]; j <= winVarR[wv]; j++) {
                            outFileVars.print("\t" + winVarNames[wv] + "_"+j);
                        }
                    } else {
                        for (int j=winVarL[wv]; j <= winVarR[wv]; j++) {
                            writeStringToFile(outFileVarsBin, winVarNames[wv] + "_" + j, charset);
                        }
                    }
                }
            }
            else {
                for (int i = -numFlankingWindows; i <= numFlankingWindows; i++) {
                    if (!outputBin) {
                        outFileVars.print("\t" + winVarNames[wv] + "_" + i);
                    } else {
                        writeStringToFile(outFileVarsBin, winVarNames[wv] + "_" + i, charset);
                    }
                }
            }
        }

        if (flankingFeatureWindowLength == 0) {
            if (!byNT) { // If we are reporting features on a by-nucleotide basis, do not genearate sequence content headers / values
                if (!outputBin) {
                    outFileVars.print("\t" + "GCcontent\tCAcontent\tGAcontent");
                } else {
                    for (int j = 0; j < seqContentNames.length; j++) {
                        String colName = seqContentNames[j];
                        writeStringToFile(outFileVarsBin, colName, charset);
                    }
                }
            }
        }
        else {
            for (int i = -numFlankingWindows; i <= numFlankingWindows; i++) {
                if (!outputBin) {
                    outFileVars.print("\tGCcontent_"+i);
                    outFileVars.print("\tCAcontent_"+i);
                    outFileVars.print("\tGAcontent_"+i);
                } else {
                    for (int j = 0; j < seqContentNames.length; j++) {
                        String colName = seqContentNames[j] + "_" + i;
                        writeStringToFile(outFileVarsBin, colName, charset);
                    }
                }
            }
        }

        // COMBINATIONS
        /*
        // Print combination labels
        for (int wv1 = 0; wv1 < nWinVars; wv1++) {
            for (int wv2 = wv1; wv2 < nWinVars; wv2++) {
                outFileVars.print("\t" + winVarNames[wv1] + "_AND_" + winVarNames[wv2]);
            }

            // and GC Content
            outFileVars.print("\tGCcontent_AND_"+winVarNames[wv1]);
        }
        */

        if (!outputBin) outFileVars.println();

        // Print results
        for (int i = 0; i < finalResults.size(); i++) {
            TagScoreResult res = finalResults.get(i);
            if (!outputBin) {
                outFileVars.print(res.sampleName);
            } else {
                writeStringToFile(outFileVarsBin, res.sampleName, charset);
            }
            for (int wv = 0; wv < res.results.length; wv++) {
                if (!outputBin) {
                    outFileVars.print("\t" + res.results[wv]);
                } else {
                    outFileVarsBin.writeDouble(res.results[wv]);
                }
            }
            if (!outputBin) {
                outFileVars.println();
            } else {
                outFileVarsBin.writeByte('!');
                outFileVarsBin.writeByte('!');
                outFileVarsBin.writeByte('!');
                outFileVarsBin.writeByte('!');
            }
        }

        if (!outputBin) {
            outFileVars.flush();
            outFileVars.close();
        } else {
            outFileVarsBin.flush();
            outFileVarsBin.close();
        }

        threadPool.shutdownNow();
        while(!threadPool.isTerminated()) {
            threadPool.awaitTermination(2, TimeUnit.SECONDS);
        }
        System.out.println("Complete");
    }

    private static int writeStringToFile(DataOutputStream out, String str, Charset charset) throws java.io.IOException {
        out.writeByte(Math.min(255, str.length()));

        byte[] columnName = Arrays.copyOf(str.getBytes(charset), Math.min(255, str.length()));
        out.write(columnName);

        return columnName.length + 1; // +1 for the length byte
    }

    private static void parseArgs(String[] args) {
        // Setup CLI to simplify GenFeatures interface;
        Options options = buildOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmdLine;

        try {
            cmdLine = parser.parse(options, args);

            // Remove arguments that were parsed
            args = new String[cmdLine.getArgList().size()];
            args = cmdLine.getArgList().toArray(args);

            if (cmdLine.hasOption("help")) {
                help = true;
            } else {
                if (args.length < 5) {
                    throw new ParseException("Must supply a FWD & REV ROE tables, FASTA file, PWM file and output file name.");
                } else {
                    factorInfo_FWD_Fname = args[0];
                    factorInfo_REV_Fname = args[1];
                    promoterSeqs_Fname = args[2];
                    pwms_Fname = args[3];
                    out_Fname = args[4];
                }

                if (cmdLine.hasOption("nprocs")) {
                   nproc = Integer.parseInt(cmdLine.getOptionValue("nprocs"));
                   if (nproc > Runtime.getRuntime().availableProcessors()) {
                       System.err.println("Warning: your system does not have " + nproc + " CPUs: setting nprocs to " + Runtime.getRuntime().availableProcessors() + ".");
                  
                       nproc = Runtime.getRuntime().availableProcessors();
                   }
                   if (nproc < 1) {
                       System.err.println("Warning: nprocs must be a value >= 1. Setting nprocs = 1");
                       nproc = 1;
                   }
                }

                if (cmdLine.hasOption("BGWIN")) {
                    BG_WIN = Integer.parseInt(cmdLine.getOptionValue("BGWIN"));
                }

                if (cmdLine.hasOption("minScores")) {
                    scoreCutoffs_Fname = cmdLine.getOptionValue("minScores");
                }

                if (cmdLine.hasOption("mapFile")) {
                    map_Fname = cmdLine.getOptionValue("mapFile");
                }

                if (cmdLine.hasOption("nucsAfterTSS")) {
                    nucsDownStream = Integer.parseInt(cmdLine.getOptionValue("nucsAfterTSS"));
                }

                if (cmdLine.hasOption("pseudoCounts")) {
                    PseudoCountsVal = Double.parseDouble(cmdLine.getOptionValue("pseudoCounts"));
                }

                if (cmdLine.hasOption("seqName")) {
                    String seqName = cmdLine.getOptionValue("seqName");
                    if (seqName.charAt(0) != '-') {
                        processPeak = seqName;
                        System.out.println("Processing Peak "+processPeak);
                    }
                }

                if (cmdLine.hasOption("pos")) {
                    posArrString = cmdLine.getOptionValue("pos");
                }

                if (cmdLine.hasOption("flankingWindows")) {
                    String[] flankingWindowParts = cmdLine.getOptionValue("flankingWindows").split(" ");

                    if (flankingWindowParts.length != 2) {
                        throw new ParseException("Incorrect number of flanking window arguments: \"\"|\"FLANK_LENGTH NUM_WIN\")");
                    }
        
                    flankingFeatureWindowLength = Integer.parseInt(flankingWindowParts[0]);
                    numFlankingWindows = Integer.parseInt(flankingWindowParts[1]);
        
                    if (flankingFeatureWindowLength > 0) {
                        if (flankingFeatureWindowLength % numFlankingWindows != 0) {
                            throw new ParseException("Flanking window length not divisible by number of windows.");
                        }
                        else if (numFlankingWindows > flankingFeatureWindowLength) {
                            throw new ParseException("More flanking windows than the flanking region length.");
                        }
                    }
                }

                if (cmdLine.hasOption("bin")) {
                    outputBin = true;
                }

                if (cmdLine.hasOption("byNT")) {
                    byNT = true;
                    flankingFeatureWindowLength = 0;
                    numFlankingWindows = 0;
                }
            }
        }
        catch (ParseException e) {
            System.out.println("Parse Error: " + e.getMessage());
            System.out.println();
            help = true;
        }
        catch (NumberFormatException e) {
            System.out.println("Unexpected formatting problem: " + e.getMessage());
            e.printStackTrace();
            help = true;
        }

        if (help) {
            OptionComparator opt_order = new OptionComparator();
            String cmdLineSyntax = "java -jar GenFeatures tfbs_llscan.jar <FWD_ROEs.table> <REV_ROEs.table> <sequence.fasta> <PWM_FILE> <Outfile.rdat>";
            String header = "\nJava tool for generating features within regions of enrichment for loglikelihood " + 
                            "scans of FASTA sequences. Note: All matrices in <PWM_FILE> are assumed to be frequency matrices - not probabilities\n\nOPTIONS\n\n";
            String footer = "";
            HelpFormatter formatter = new HelpFormatter();
            formatter.setOptionComparator(opt_order);
            formatter.setWidth(120);
            formatter.printHelp(cmdLineSyntax, header, options, footer, true);
            System.exit(0);
        }
    }

    private static class OptionComparator <T extends Option> implements Comparator <T> {
        private ArrayList <String> OPTS_ORDER = new ArrayList <String> ();

        public OptionComparator () {
            OPTS_ORDER.add("help");
            OPTS_ORDER.add("pos");
            OPTS_ORDER.add("nprocs");
            OPTS_ORDER.add("BGWIN");
            OPTS_ORDER.add("minScores");
            OPTS_ORDER.add("mapFile");
            OPTS_ORDER.add("pseudoCounts");
            OPTS_ORDER.add("nucsAfterTSS");
            OPTS_ORDER.add("seqName");
            OPTS_ORDER.add("flankingWindows");
            OPTS_ORDER.add("bin");
            OPTS_ORDER.add("byNT");
        }

        public int compare(T o1, T o2) {
            return OPTS_ORDER.indexOf(o1.getLongOpt()) - OPTS_ORDER.indexOf(o2.getLongOpt());
        }
    }

    private static Options buildOptions() {
        Options options = new Options();

        // help
        options.addOption(Option.builder("h")
            .longOpt("help")
            .hasArg(false)
            .desc("Print this helpful help message")
            .required(false)
            .build());

        // the number of processors to use
        options.addOption(Option.builder("n")
            .longOpt("nprocs")
            .hasArg(true)
            .argName("INT")
            .desc("Number of processors to use (default: 1)")
            .required(false)
            .type(Integer.class)
            .build());

        // string to specify where to center scans for generating features
        options.addOption(Option.builder()
            .longOpt("pos")
            .hasArg(true)
            .argName("Position String")
            .desc("Position String: String describing where to scan from. One of:\n" +
                  "\"N\": A number indicating the location to sample from (relative to TSS of sequence)\n\n" +
                  "\"Rand [ndraws] [draw_start] [draw_stop]\": Choose [ndraws] random locations between [draw_start] and [draw_stop])\n\n" +
			      "\"Range [i] [range_start] [range_stop]\": Scan every [i] nucleotides between [range_start] and [range_stop])\n" + 
			      "\"Replicate [filename]\": Scan every nucleotide defined in [filename])\n" + 
                  "(default: 0 - i.e., only generate features from TSS of sequence)")
            .required(false)
            .build());

        // sequence to limit scans to
        options.addOption(Option.builder()
            .longOpt("seqName")
            .hasArg(true)
            .argName("SEQNAME")
            .desc("Limit ROE scans to sequence <SEQNAME>.")
            .required(false)
            .build());

        // the offset to adjust the TSS position in the input sequence
        options.addOption(Option.builder("N")
            .longOpt("nucsAfterTSS")
            .hasArg(true)
            .argName("INT")
            .desc("Offset from end point of input sequence for TSS (default: midpoint of input sequence)")
            .required(false)
            .type(Integer.class)
            .build());

        // the pseudo-count value to adjsut PWM frequencies by
        options.addOption(Option.builder()
            .longOpt("pseudoCounts")
            .hasArg(true)
            .argName("COUNT")
            .desc("Total pseudo-counts to add to PWM matrix prior to normalization (default: 0.25)")
            .required(false)
            .type(Double.class)
            .build());
            
        // length of background window used for local background sequence calculations
        options.addOption(Option.builder("B")
            .longOpt("BGWIN")
            .hasArg(true)
            .argName("INT")
            .desc("Background Window size used for calculating local sequence background distribution (default: 250)")
            .required(false)
            .build());

        // score threshold file
        options.addOption(Option.builder()
            .longOpt("minScores")
            .hasArg(true)
            .argName("FILENAME")
            .desc("file to read mininum threshold scores for scans to be added (default: all thresholds set to 0)")
            .required(false)
            .build());

        // score threshold file
        options.addOption(Option.builder()
            .longOpt("mapFile")
            .hasArg(true)
            .argName("FILENAME")
            .desc("file to output how scanned features are mapped to genomic coordinates.  (Note: requires that FASTA headers of your input sequences are formatted as follows \"><FASTAID>_<REFSEQ>_<START>\" " +
                  "where  <REFSEQ> is the reference sequence this FASTA sequence comes from, and <START> is the integer location where the sequence begins (1-indexed)")
            .required(false)
            .build());

        // flanking windows string for generating features from flanking sequence
        options.addOption(Option.builder()
            .longOpt("flankingWindows")
            .hasArg(true)
            .argName("WindowLength NumWindows")
            .desc("If specified, generates additional features for each TSS centered upstream and downstream. " +
                  "Generates NumWindows additional features on each side of TSS, going up to WindowLength away from site." +
                  "(default: no flanking windows are used)")
            .required(false)
            .build());

        // Output data in binary format
        options.addOption(Option.builder("b")
            .longOpt("bin")
            .hasArg(false)
            .desc("Output data in binary format.")
            .required(false)
            .build());

        // Output scores on a per-nucleotide basis for each ROE scanned
        options.addOption(Option.builder("y")
            .longOpt("byNT")
            .hasArg(false)
            .desc("Output scores on a per-nucleotide basis for each ROE scanned.  Turns off reporting of flanking features and sequence content features (i.e. GC/GA/CA content)")
            .required(false)
            .build());

        return options;
    }

    private static class TagScoreResult {
    	double[] results;
        long numberNonzeros = 0;
    	String sampleName;
    	
		public TagScoreResult(double[] results, String sampleName, long numberNonzeros) {
			super();
			this.results = results;
			this.sampleName = sampleName;
            this.numberNonzeros = numberNonzeros;
		}
    	
    }
    
    private static class TagScoresRunner implements Callable<TagScoreResult>{
    	private String sampleName;
    	
    	char[] S;
    	DoubleColReturn Scores;
		int nucsAfterTSS;
    	int BG_WIN;
    	String seqName;
    	PWMReturn pwms;
		int nWinVars;
    	int[] winVarL;
    	int[] winVarR;
    	int[] winVarPWMIndex;
		String[] winVarStrand;
    	Random rand;
    	boolean drawRand;
		boolean drawRep;
    	int drawRandN;
    	int drawRandL;
    	int drawRandR;
		Hashtable repHash;
    	int[] posArr;
    	int i;
    	int p;
        int flankingFeatureWindowLength;
        int numFlankingWindows;
        long numberNonzeros = 0;

		public TagScoresRunner(char[] s, DoubleColReturn scores,
				int nucsAfterTSS, int bG_WIN, String seqName, PWMReturn pwms,
				int nWinVars, int[] winVarL, int[] winVarR,
				int[] winVarPWMIndex, String[] winVarStrand, Random rand,
				boolean drawRand, boolean drawRep, int drawRandN,
				int drawRandL, int drawRandR, Hashtable repHash, int[] posArr,
				int i, int p, int flankingFeatureWindowLength, int numFlankingWindows) {
			super();
			S = s;
			Scores = scores;
			this.nucsAfterTSS = nucsAfterTSS;
			BG_WIN = bG_WIN;
			this.seqName = seqName;
			this.pwms = pwms;
			this.nWinVars = nWinVars;
			this.winVarL = winVarL;
			this.winVarR = winVarR;
			this.winVarPWMIndex = winVarPWMIndex;
			this.winVarStrand = winVarStrand;
			this.rand = rand;
			this.drawRand = drawRand;
			this.drawRep = drawRep;
			this.drawRandN = drawRandN;
			this.drawRandL = drawRandL;
			this.drawRandR = drawRandR;
			this.repHash = repHash;
			this.posArr = posArr;
			this.i = i;
			this.p = p;
            this.flankingFeatureWindowLength = flankingFeatureWindowLength;
            this.numFlankingWindows = numFlankingWindows;
		}

		public TagScoreResult call() throws BadCharException {
			// Compute and store local background values for this sequence
			double[][] B = Utils.getWholeSeqLocalBackground(S, BG_WIN);
			double[][][] B_M1 = Utils.getWholeSeqLocalM1Background(S, BG_WIN);
	
			int pos = posArr[p];
	
            int numFeatures = nWinVars + 3; // number of windows, plus GC, CA, GA contents
            if (this.numFlankingWindows > 0) {
                numFeatures = (nWinVars + 3)*(2*numFlankingWindows+1); // +1 for center window
            }

			// Calculate TF window score variables
			//double[] featureVals = new double[nWinVars + 3];

			double[] featureVals = null;
            ArrayList <Double> featureValsArray = new ArrayList<Double>();
            if (!byNT) featureVals = new double[numFeatures];

			for (int wv = 0; wv < nWinVars; wv++) {
			    int L = winVarL[wv];
			    int R = winVarR[wv];
			    int pwmIdx = winVarPWMIndex[wv];
			    String strand = winVarStrand[wv];
	
			    double[][] ThetaT;
			    // SHAWN
			   	//System.out.print("Label: " + pwms.labels[pwmIdx] + " " + seqName + " ");
			    if (strand.equals("FWD")) {
			        ThetaT = pwms.pwms[pwmIdx];
			    } else { // REV
			        ThetaT = pwms.revcmp_pwms[pwmIdx];
			    }
	
			        // SHAWN for debugging
					/*for(int rowindex = 0; rowindex < ThetaT.length; rowindex++) {
						double[] row = ThetaT[rowindex]; 
						for(int colindex = 0; colindex < row.length; colindex++) {
							double el = row[colindex];
							el = ((int)(el*100.0))/100.0;
							System.out.print(el + " ");
						}
						System.out.println();
					}*/

                for (int flankWin = -numFlankingWindows; flankWin <= numFlankingWindows; flankWin++) {
                    int featureIndex = wv*(2*numFlankingWindows+1) + (flankWin+numFlankingWindows);

                    int flankingPos = 0;
                    if (numFlankingWindows > 0) {
                        flankingPos = flankWin*(flankingFeatureWindowLength/numFlankingWindows);
                    }

                    if (!byNT) {
                        featureVals[featureIndex] = cumScoreInWindow(pos + L + flankingPos, pos + R + flankingPos, S, ThetaT, strand,
                                B, B_M1, Scores.values[pwmIdx], nucsAfterTSS);

                        if (featureVals[featureIndex] != 0) numberNonzeros++;
                    } else {
                        double[] ntScores = scoresInWindow(pos + L + flankingPos, pos + R + flankingPos, S, ThetaT, strand,
                                B, B_M1, Scores.values[pwmIdx], nucsAfterTSS);

                        for (int i = 0; i < ntScores.length; i++) {
                            featureValsArray.add(ntScores[i]);
                            if (ntScores[i] != 0) numberNonzeros++;
                        }
                    }
                }
			}
	
            if (!byNT) {
    			// Calculate GC Content
    			int GC_WIN = 100;
    
                for (int flankWin = -numFlankingWindows; flankWin <= numFlankingWindows; flankWin++) {
                    int featureIndex = nWinVars*(2*numFlankingWindows+1) + 3*(flankWin+numFlankingWindows);
    
                    int flankingPos = pos;
    
                    if (numFlankingWindows > 0) {
                        flankingPos += flankWin*(flankingFeatureWindowLength/numFlankingWindows);
                    }
                    
                    double gccontent = getGCContent(flankingPos - GC_WIN, flankingPos + GC_WIN, S, nucsAfterTSS);
                    featureVals[featureIndex] = gccontent; // set last element to contain GC content
                    if (featureVals[featureIndex] != 0) numberNonzeros++;
    
                    double CAcontent = getSeqContent(flankingPos - GC_WIN, flankingPos + GC_WIN, S, nucsAfterTSS, new char[]{'C', 'A'});
                    featureVals[featureIndex+1] = CAcontent; 
                    if (featureVals[featureIndex+1] != 0) numberNonzeros++;
    
                    double GAcontent = getSeqContent(flankingPos - GC_WIN, flankingPos + GC_WIN, S, nucsAfterTSS, new char[]{'G', 'A'});
                    featureVals[featureIndex+2] = GAcontent; 
                    if (featureVals[featureIndex+2] != 0) numberNonzeros++;
                }
            } else {
                featureVals = new double[featureValsArray.size()];
                for (int i = 0; i < featureValsArray.size(); i++) {
                    featureVals[i] = featureValsArray.get(i);
                }
            }
	
			sampleName = seqName + "_" + pos;
			    
			/*
			// Print var line
			String seqLab = seqName + "_" + pos;
			outFileVars.print(seqLab);
			for (int wv = 0; wv < nWinVars; wv++) {
			    outFileVars.print("\t" + featureVals[wv]);
			}
			outFileVars.print("\t" + gccontent);
			*/
			
			// COMBINATIONS
			/*
			// Print combinations
			for (int wv1 = 0; wv1 < nWinVars; wv1++) {
			    for (int wv2 = wv1; wv2 < nWinVars; wv2++) {
			        outFileVars.print("\t" + (winVarVals[wv1] * winVarVals[wv2]));
			    }
	
			    // and GC Content
			    outFileVars.print("\t" + winVarVals[wv1] * gccontent);
			}
			*/
	
			/*
			outFileVars.println();
			outFileVars.flush();
			*/
	
			return new TagScoreResult(featureVals, sampleName, numberNonzeros);
		}
    }

    public static double[] scoresInWindow(int L, int R, char[] S, double[][] ThetaT, String strand, double[][] B, double[][][] B_M1, double Score, int beyondTSS) throws BadCharException {
        int w = ThetaT.length;
        int chrInd, chrInd_prev;
        double theta_freq, bg_freq;
        double numLogSum, denLogSum;

        boolean hitMaskedChar;
        boolean USE_MARKOV_1_BG = false;

        if (B_M1 != null) {
            USE_MARKOV_1_BG = true;
        }

        // Convert window positions to array indices
        int tssLoc = S.length - beyondTSS;
        int winL = L + tssLoc;
        int winR = R + tssLoc;

        if (winL < 0) { winL = 0; }
        if (winR > S.length - w) { winR = S.length - w; }

        double scores[] = new double[winR - winL + 1]; // scores for nt in the window

        // Compute scores for each possible motif site in window
        for (int j = winL; j <= winR; j++) {
            numLogSum = 0.0;
            denLogSum = 0.0;
            // Take product of freq/bg probabilities over motif width
            hitMaskedChar = false;
            for (int mp = 0; mp < w; mp++) {
                chrInd = Setup.charInd(S[j+mp]);
                if (chrInd < 0) { // encountered a masked character
                    hitMaskedChar = true;
                    break;
                }
                theta_freq = ThetaT[mp][chrInd];
                if (USE_MARKOV_1_BG) {
                    if (mp == 0) {
                        bg_freq = B[j][chrInd];
                    } else {
                        chrInd_prev = Setup.charInd(S[j + mp - 1]);
                        bg_freq = B_M1[j][chrInd_prev][chrInd];
                    }
                }
                else {
                    bg_freq = B[j][chrInd];
                }
                numLogSum += Math.log(theta_freq);
                denLogSum += Math.log(bg_freq);
            }

            if (!hitMaskedChar) {
                double score = numLogSum - denLogSum;
                if (score > Score) {
                    scores[j - winL] = score; 
                }
                if (Double.isInfinite(score)) {
                    System.out.println(
                            "WARNING: Encountered an infinite score at position " + j);
                }
            }
        }

        return scores;
    }

    public static double cumScoreInWindow(int L, int R, char[] S, double[][] ThetaT, String strand, double[][] B, double[][][] B_M1, double Score, int beyondTSS) throws BadCharException {
        int w = ThetaT.length;
        int chrInd, chrInd_prev;
        double theta_freq, bg_freq, score;
        double numLogSum, denLogSum;

        double cumScore = 0.0;

        boolean hitMaskedChar;
        boolean USE_MARKOV_1_BG = false;

        if (B_M1 != null) {
            USE_MARKOV_1_BG = true;
        }

        // Convert window positions to array indices
        int tssLoc = S.length - beyondTSS;
        int winL = L + tssLoc;
        int winR = R + tssLoc;

        if (winL < 0) { winL = 0; }
        if (winR > S.length - w) { winR = S.length - w; }

        // Compute scores for each possible motif site in window
        for (int j = winL; j <= winR; j++) {
            numLogSum = 0.0;
            denLogSum = 0.0;
            // Take product of freq/bg probabilities over motif width
            hitMaskedChar = false;
            for (int mp = 0; mp < w; mp++) {
                chrInd = Setup.charInd(S[j+mp]);
                if (chrInd < 0) { // encountered a masked character
                    hitMaskedChar = true;
                    break;
                }
                theta_freq = ThetaT[mp][chrInd];
                if (USE_MARKOV_1_BG) {
                    if (mp == 0) {
                        bg_freq = B[j][chrInd];
                    } else {
                        chrInd_prev = Setup.charInd(S[j + mp - 1]);
                        bg_freq = B_M1[j][chrInd_prev][chrInd];
                    }
                }
                else {
                    bg_freq = B[j][chrInd];
                }
                numLogSum += Math.log(theta_freq);
                denLogSum += Math.log(bg_freq);
            }
            if (!hitMaskedChar) {
                score = numLogSum - denLogSum;
                if (score > Score) {
                    // Record cumulative score
                    cumScore += score;
                }
                if (Double.isInfinite(score)) {
                    System.out.println(
                            "WARNING: Encountered an infinite score at position " + j);
                }
            }
        }
		// SHAWN
		//System.out.println("Strand: " + strand +  " L: " + L + " R: " + R + " beyondTSS: " + beyondTSS + " score: " + cumScore);
        return cumScore;

    }

    public static double getSeqContent(int L, int R, char[] S, int beyondTSS, char[] bases) {

        // Convert window positions to array indices
        int tssLoc = S.length - beyondTSS;
        int winL = L + tssLoc;
        int winR = R + tssLoc;

        if (winL < 0) { winL = 0; }
        if (winR > S.length) { winR = S.length; }

        // Compute gcContent in window
        double seqContent = Utils.seqPercent(String.valueOf(S), winL, winR, bases);

        return seqContent;

    }

    public static double getGCContent(int L, int R, char[] S, int beyondTSS) {
        // Convert window positions to array indices
        int tssLoc = S.length - beyondTSS;
        int winL = L + tssLoc;
        int winR = R + tssLoc;

        if (winL < 0) { winL = 0; }
        if (winR > S.length) { winR = S.length; }

        // Compute gcContent in window
        double gcContent = Utils.gcPercent(String.valueOf(S), winL, winR);

        return gcContent;

    }

}
