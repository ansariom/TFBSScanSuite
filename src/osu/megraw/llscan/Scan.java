package osu.megraw.llscan;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Hashtable;
import java.util.List;
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
 * Copyright (C) 2010  Molly Megraw
 *
 * This program is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *
 * Contact: molly.megraw@duke.edu
 *
 */
public class Scan {
    // Command Line Argument Variables

    // Set variable which represents the number of nucleotides after the
    // TSS represented in the FASTA file of sequences.  For example,
    // if the TSS is at the "far right" of each sequence, then set
    // nucsAfterTSS = 0.  If the TSS is 500 nucleotides from the "far right"
    // of each sequence, then set nucsAfterTSS = 500.
    private static int nucsAfterTSS = 0;
    private static int nprocs = 1;

    private static char[][] S;
    private static PWMReturn pwms;
    private static DoubleColReturn minScores;
    private static Background BG;
    private static double[][] bg;
    private static double[][][] bg2;

    private static double PseudoCountsVal = 0.25;
    private static String strand = "BOTH";
    private static String [] strands = {"FWD", "REV"};
    private static String minScores_Fname = null;
    private static String pwms_Fname = null;
    private static String seqs_Fname = null;
    private static String out_Fbase = null;
    private static String bg_Fname = null;
    private static String bg2_Fname = null;
    private static boolean help = false;
    private static boolean usebg2 = false;

    public Scan() {
    }

    public static void main(String[] args) throws java.io.IOException, java.lang.InterruptedException {
        Scan scansites = new Scan();

        parseArgs(args);

        // Read promoter sequences file
        FastaReturn Seqs = Load.loadFastaFile(new File(seqs_Fname));
        S = new char[Seqs.lines.length][];
        for (int i = 0; i < Seqs.lines.length; i++) {
            S[i] = Seqs.lines[i].toCharArray();
        }

        // Get labels for seqs
        String[] seqLabels = new String[Seqs.headers.length];
        for (int i = 0; i < seqLabels.length; i++) {
            String header = Seqs.headers[i].substring(0, Seqs.headers[i].length());
            String[] header_parts = header.split("\\s+");
            seqLabels[i] = header_parts[0];
        }

        // Load background model(s)
        if (bg_Fname != null) {
            bg = new double[1][];
            bg[0] = Load.loadVector(new File(bg_Fname));
        } else {
            bg = new double[1][];
            bg[0] = Background.getBackground_M0(S);
        }
        if (usebg2) {
            bg2 = new double[1][][];
            if (bg2_Fname != null) {
                bg2[0] = Load.loadMatrix(new File(bg2_Fname));
            } else {
                bg2[0] = Background.getBackground_M1(S);
            }
            System.out.println("Using 1st Order Markov Background");
        }

        // false lets the scanner know I'm not using local sequence background distributions
        BG = new Background("", bg, bg2, false);

        // Read PWM file
        //pwms = Load.loadPWMFile(pwms_Fname, PseudoCountsVal);
        pwms = Load.loadPWMFileSimpleHeader(pwms_Fname, PseudoCountsVal);

        // Read score threshold file corresponding to PWM file
        if (minScores_Fname != null) {
            minScores = Load.loadDoubleColumnFile(new File(minScores_Fname));
        } else {
            // No threshold was given - assume that threshold will be zero
            minScores = new DoubleColReturn(pwms.labels);
        }

        // Setup thread management for parallel processing
        ThreadPoolExecutor threadPool = new ThreadPoolExecutor(nprocs, nprocs, 4, TimeUnit.SECONDS, new LinkedBlockingQueue());
        List<Future<ScanResult>> results = new ArrayList<Future<ScanResult>>();
        List<ScanResult> finalResults = new ArrayList<ScanResult>();

        int totalStrands = strand.equals("BOTH")? 2 : 1;

        // Do all the scans!!!!
        for (int i = 0; i < totalStrands; i++) {
            // Only reset strand to ensure both are used if user explictly requested both strands be scanned
            if (totalStrands == 2) strand = strands[i];

            for (int nmat = 0; nmat < pwms.pwms.length; nmat++) {
                for (int nseq = 0; nseq < S.length; nseq++) {
                    //ScanRunner run = new ScanRunner(S[nseq], pwms.pwms[nmat], strand, bg, bg2, minScores.values[nmat], nucsAfterTSS, 
                    //                                pwms.labels[nmat], seqLabels[nseq]);
                    ScanRunner run = new ScanRunner(S[nseq], pwms.pwms[nmat], strand, BG, minScores.values[nmat], nucsAfterTSS, 
                                                    pwms.labels[nmat], seqLabels[nseq]);

                    try {
                        results.add(threadPool.submit(run));
                    } catch (Exception e) {
                        System.err.println("Error processing sequencing " +
                                            seqLabels[nseq] + 
                                           " with PWM " +
                                           pwms.labels[nmat] + 
                                           " on strand " + 
                                           strand + 
                                           ": " +
                                          e.getMessage() + 
                                          ". Quitting.");
                        threadPool.shutdownNow();
                        while(!threadPool.isTerminated()) {
                            threadPool.awaitTermination(2, TimeUnit.SECONDS);
                        }
                        System.exit(1);
                    }
                }
            }
        }

        // Harcode sorting by PWM in simple scanner, but could add command line arguments in
        // future versions if you wanted to add support for different output formats
        String sortBy = "PWM";

        // Hashtable to store header results, i.e. the number of sequences with results for a PWM if sorting by PWM,
        // or the number of PWMS with results for a sequence if sorting by sequence
        Hashtable <String, MutableInteger> headerResults = new Hashtable <String, MutableInteger>(); 

        // Obtain all the scanned results
        for (int i = 0; i < results.size(); i++) {
            Future <ScanResult> result = results.get(i);
            ScanResult res = null;
            try {
                res = result.get(); // get() blocks until the result is available
                // Only keep scans that had results
                if (res.hitLocs.length > 0) {
                    finalResults.add(res);

                    // Generate headers by strand and PWM
                    String headerKey = "";

                    if (sortBy.equals("PWM")) {
                        headerKey = res.strand + "_" + res.pwmLabel;
                    // Generate headers by SEQ - assume this is the case if "PWM" isn't found
                    } else {
                        headerKey = res.strand + "_" + res.seqLabel;
                    }

                    // initialize count at 0
                    if (!headerResults.containsKey(headerKey)) {
                        headerResults.put(headerKey, scansites.new MutableInteger(0));
                    }

                    // Keep track of how many Seqs produce results for this PWM
                    headerResults.get(headerKey).increment();
                }
            } catch (Exception e) {
                System.err.println(e.getMessage());
                e.printStackTrace();
                break;
            }
        }

        // Scans didn't produce any results - shutdown threadpool and exit the application
        // (Make sure to let user know nothing was found)
        if (finalResults.size() == 0) {
            // Clean up our thread mess
            threadPool.shutdownNow();
            while(!threadPool.isTerminated()) {
                threadPool.awaitTermination(2, TimeUnit.SECONDS);
            }
            System.out.println("No scans resulted in identified TFBSs.");
            return;
        }

        // Keep track of current PWM/Seq header we are printing as we loop through results
        String headerKey = "";
        String header = "";

        // Group results by strand and pwm/seq since the final array is a list of scans done separately for each
        // individual strand, pwm, and sequence
        // Results are always sorted by strand first (FWD then REV), but for this output, 
        // ensure that I am printing out results sorted by PWM if the strands are equal (SEQ sorts by sequence scanned)
        ScanResultComparator scanResultComparator = new ScanResultComparator(sortBy);
        Collections.sort(finalResults, scanResultComparator);

        // Print out all the scan results!!!
        strand = finalResults.get(0).strand;
        PrintWriter outFileHitsDetail = new PrintWriter(new FileWriter(out_Fbase + "." + strand + ".locations"));
        PrintWriter outFileScoresDetail = new PrintWriter(new FileWriter(out_Fbase + "." + strand + ".scores"));
        for (int i = 0; i < finalResults.size(); i++) {
            ScanResult result = finalResults.get(i);

            // We've passed the FWD strand results - start outputting to REV strand results file
            if (!result.strand.equals(strand)) {
                outFileHitsDetail.close();
                outFileScoresDetail.close();
                strand = result.strand;
                outFileHitsDetail = new PrintWriter(new FileWriter(out_Fbase + "." + strand + ".locations"));
                outFileScoresDetail = new PrintWriter(new FileWriter(out_Fbase + "." + strand + ".scores"));
            }

            // Print out headers when we encounter a new header - 'start' header is empty string
            // strand + "_" + pwmLabel / seqLabel is the format of the header keys since we have
            // separate files/headers for FWD and REV strands
            if ( (!headerKey.equals(strand + "_" + finalResults.get(i).pwmLabel) && sortBy.equals("PWM")) ||
                 (!headerKey.equals(strand + "_" + finalResults.get(i).seqLabel) && sortBy.equals("SEQ")) ){
                if (sortBy.equals("PWM")) {
                    headerKey = strand + "_" + finalResults.get(i).pwmLabel;
                    header = finalResults.get(i).pwmLabel;
                } else {
                    headerKey = strand + "_" + finalResults.get(i).seqLabel;
                    header = finalResults.get(i).seqLabel;
                }

                outFileHitsDetail.println(">\t" + header + "\t" + headerResults.get(headerKey).intValue());
                outFileScoresDetail.println(">\t" + header + "\t" + headerResults.get(headerKey).intValue());
            }

            // Print out scann results for this PWM/Seq group
            StringBuilder locResultsData = new StringBuilder();
            StringBuilder scoreResultsData = new StringBuilder();
            if (sortBy.equals("PWM")) {
                locResultsData.append(result.seqLabel);
                scoreResultsData.append(result.seqLabel);
            } else {
                locResultsData.append(result.pwmLabel);
                scoreResultsData.append(result.pwmLabel);
            }
            for (int j = 0; j < result.hitLocs.length; j++) {
                locResultsData.append("\t" + result.hitLocs[j]);
                scoreResultsData.append("\t" + Print.df3.format(result.hitScores[j]));
            }     
            outFileHitsDetail.println(locResultsData.toString());
            outFileScoresDetail.println(scoreResultsData.toString());
        }
        outFileHitsDetail.close();
        outFileScoresDetail.close();

        // Clean up our thread mess
        threadPool.shutdownNow();
        while(!threadPool.isTerminated()) {
            threadPool.awaitTermination(2, TimeUnit.SECONDS);
        }
        System.out.println("Scans complete!");
    }

    // primary function for verifying command line arguments
    private static void parseArgs(String[] args) {
        // Setup CLI to simplify Scan interface;
        Options options = ScanOptions.buildOptions();
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
	            if (args.length < 3) {
	                throw new ParseException("Must supply PWM file, FASTA file, and output basename.");
	            } else {
	                pwms_Fname = args[0];
	                seqs_Fname = args[1];
	                out_Fbase = args[2];
	            }

                if (cmdLine.hasOption("nprocs")) {
                   nprocs = Integer.parseInt(cmdLine.getOptionValue("nprocs"));
                   if (nprocs > Runtime.getRuntime().availableProcessors()) {
                       System.err.println("Warning: your system does not have " + nprocs + " CPUs: setting nprocs to " + Runtime.getRuntime().availableProcessors() + ".");
                       nprocs = Runtime.getRuntime().availableProcessors();
                   }
                   if (nprocs < 1) {
                       System.err.println("Warning: nprocs must be a value >= 1. Setting nprocs = 1");
                       nprocs = 1;
                   }
                }

	            if (cmdLine.hasOption("bg")) {
	                bg_Fname = cmdLine.getOptionValue("bg");
	            }

	            if (cmdLine.hasOption("bg2")) {
	                bg2_Fname = cmdLine.getOptionValue("bg2");
                    usebg2 = true;
	            } else if (cmdLine.hasOption("usebg2")) {
                    usebg2 = true;
                }

	            if (cmdLine.hasOption("strand")) {
	                strand = cmdLine.getOptionValue("strand");
	                if (!strand.equals("FWD") && !strand.equals("REV") && !strand.equals("BOTH")) {
	                    strand = "BOTH";
	                }
	            }

	            if (cmdLine.hasOption("minScores")) {
	                minScores_Fname = cmdLine.getOptionValue("minScores");
	            }

	            if (cmdLine.hasOption("pseudoCounts")) {
	                PseudoCountsVal = Double.parseDouble(cmdLine.getOptionValue("pseudoCounts"));
	            }

	            if (cmdLine.hasOption("nucsAfterTSS")) {
	                nucsAfterTSS = Integer.parseInt(cmdLine.getOptionValue("nucsAfterTSS"));
	            }
            }
        }
        catch (ParseException e) {
            System.out.println("Parse Error: " + e.getMessage());
            help = true;
        }
        catch (NumberFormatException e) {
            System.out.println("Unexpected formatting problem: " + e.getMessage());
            e.printStackTrace();
            System.exit(1);
        }

        if (help) {
            ScanOptions so = new ScanOptions();
            Comparator <Option> opt_order = so.new OptionComparator <Option> ();
            String cmdLineSyntax = "java -jar tfbs_llscan.jar Scan <PWM FILE> <FASTA FILE> <FILEBASE>";
            String header = "\nJava tool for generating loglikelihood scans of FASTA sequences. Note: All matrices in <PWM_FILE> are assumed to be frequency matrices - not probabilities\n\nOPTIONS\n\n";
            String footer = "";
            HelpFormatter formatter = new HelpFormatter();
            formatter.setOptionComparator(opt_order);
            formatter.setWidth(80);
            formatter.printHelp(cmdLineSyntax, header, options, footer, true);
            System.exit(0);
        }
    }

    // Simple helper class to keep track of header counts inside of a Hashtable used in main function
    public class MutableInteger {
        private int value;

        public MutableInteger() {
            this.value = 0;
        }

        public MutableInteger(int value) {
            this.value = value;
        }

        public void set(int value) {
            this.value = value;
        }

        public void increment() {
            this.value++;
        }

        public void decrement() {
            this.value--;
        }

        public int intValue() {
            return this.value;
        }
    }
}
