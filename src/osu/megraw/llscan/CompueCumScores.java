package osu.megraw.llscan;

import java.io.File;
import java.io.FileWriter;
import java.io.PrintWriter;
import java.math.BigDecimal;
import java.math.MathContext;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.Collections;
import java.util.Comparator;
import java.util.HashSet;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.List;
import java.util.Set;
import java.util.concurrent.Callable;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
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
 * <p>Copyright: Copyright (c) 2005</p>
 *
 * <p>Company: </p>
 *
 * @author Molly Megraw
 * @version 1.0
 */
public class CompueCumScores {

    public static char[][] S;
    public static DoubleColReturn Scores;

    // Optional command line arguments / variables
    public static boolean help = false;
    public static int nucsAfterTSS = -1;
    public static int BG_WIN = 250;
    public static String[] strands = {"FWD", "REV"};
    public static String strand = "BOTH";
    public static String scoreCutoffs_Fname = null;
    public static int threadPoolSize = 1;
    public static int seqLength = -1;
    public static int maxUpstreamDist = -1;
    public static int maxDownstreamDist = -1;
    public static double pseudoCountsVal = 0.25;
    public static String seqName = null; // Limit sequences scanned to one sequence
    public static String plotDir = null; // Directory to store plot files made for ROEs

    // Required command line arguments
    public static String pwms_Fname;
    public static String promoterSeqs_Fname;
    public static String out_Fname;

    public CompueCumScores() {
    }

    public static void main(String[] args) throws java.io.IOException, BadCharException, InterruptedException, ExecutionException {
        CompueCumScores ss = new CompueCumScores();

        // Slurp in all the command line arguments
        parseArgs(args);

        ExecutorService pool = Executors.newFixedThreadPool(threadPoolSize);

        // Setup output directories for plotting if user specified to make plots
//        if (plotDir != null) setupPlotDir();

        // Read promoter sequences file
        FastaReturn Seqs;
        Seqs = Load.loadFastaFile(new File(promoterSeqs_Fname));
        char[][] allS = new char[Seqs.lines.length][];
        for (int i = 0; i < Seqs.lines.length; i++) {
            allS[i] = Seqs.lines[i].toCharArray();
        }
        // Get labels for seqs
        String[] allSeqLabels = new String[Seqs.headers.length];
        for (int i = 0; i < allSeqLabels.length; i++) {
            String header = Seqs.headers[i].substring(0, Seqs.headers[i].length());
            String[] header_parts = header.split("\\s+");
            allSeqLabels[i] = header_parts[0];
        }


        String[] seqLabels;
        if (seqName == null || seqName.equals("-")) {
            S = allS;
            seqLabels = allSeqLabels;
        }

        // Code that handles keeping an individual sequence from the promoter sequence file
        else {
            int index=0;
            for (index = 0; index < Seqs.headers.length; index++) {
                if (seqName.equals(allSeqLabels[index])) {
                    break;
                }
            }

            if (index == Seqs.headers.length) {
                System.err.println("Could not find sequence "+seqName);
                return;
            }
            else {
                S = new char[1][];
                seqLabels = new String[1];

                S[0] = allS[index];
                seqLabels[0] = allSeqLabels[index];

                allS = null;
                allSeqLabels = null;
            }
        }

        // Ensure all sequences are the same length to avoid weird errors
        for (int i = 0; i < S.length; i++) {
            if (seqLength == -1) seqLength = S[i].length;

            if (S[i].length != seqLength) {
                System.err.println("Error: input sequences *x* be the same length.");
                System.exit(0);
            }
        }

        // no TSS offset given - autoset it to midpoint of input sequence length
        if (nucsAfterTSS == -1) {
            // For odd length sequences, half the length of the sequence rounded down will be midpoint since
            // sequences are 0-indexed (i.e. a sequence 10001 nt long needs nucsAfterTSS=5000 to set the TSS
            // at the midpoint, i.e. character 5001 in a 1-index based setting)
            nucsAfterTSS = seqLength / 2;  // integer arithmetic - automatically floors
        }

        // Variable included in case future developers want to make this change dynamically
        // (i.e. make this a variable set at run-time at the command line)
        int minFlanking = 1000;

        // Ensure that I use atleast 1kb of the flanking sequence to calculate background distributions
        maxUpstreamDist = (seqLength - (nucsAfterTSS + 1) - minFlanking) * -1;
        maxDownstreamDist = (nucsAfterTSS + 1) - minFlanking;

        if (maxUpstreamDist >= 0 || maxDownstreamDist <= 0) {
            System.err.println("Error: input sequences must have 1KB of flanking sequence to calculate the background distribution.");
            System.exit(0);
        }

        // Generate local sequence background distributions for all our promoter sequences
        // NOTE: Orignal ROE Finder code *DID NOT HANDLE BACKGROUND MODELS OF 0th ORDER*
        //       If you want to do that - you will need to changes this section by setting B_M1
        //       to null - this should instantiate the Background object to assume a 0th order
        //       when ScanRunner uses it.

        System.out.println("Generating local background sequence distributions");

        // Generate local sequence background distributions for all our promoter sequences
        Hashtable <String, Background> bgModels = new Hashtable <String, Background>();
        if (BG_WIN > 0) {
            for (int i = 0; i < seqLabels.length; i++) {
                int current = i + 1;
                System.err.print("\rGetting BG for " + seqLabels[i] + ": " + current + " / " + seqLabels.length); 
                double[][] B = Utils.getWholeSeqLocalBackground(S[i], BG_WIN);
                double[][][] B_M1 = Utils.getWholeSeqLocalM1Background(S[i], BG_WIN);
                bgModels.put(seqLabels[i], new Background(seqLabels[i], B, B_M1));
                if (i < seqLabels.length - 1) {
                    System.err.print("\r                                                               ");
                } else {
                    System.err.println();
                }
            }
        } else {
            bgModels.put("", new Background()); // Store equal background model under an empty string
        }

        // Read PWM file and process each entry
        PWMReturn pwms = Load.loadPWMFileSimpleHeader(pwms_Fname, pseudoCountsVal);

        // Read score threshold file corresponding to PWM file
        if (scoreCutoffs_Fname != null) {
            Scores = Load.loadDoubleColumnFile(new File(scoreCutoffs_Fname));
        } else {
            // No threshold was given - assume that threshold will be zero
            Scores = new DoubleColReturn(pwms.labels);
        }

        // Setup / run all the scans!
        ThreadPoolExecutor threadPool = new ThreadPoolExecutor(threadPoolSize, threadPoolSize, 4, TimeUnit.SECONDS, new LinkedBlockingQueue());
        List<Future<ScanResult>> futResults = new ArrayList<Future<ScanResult>>();

        List<ScanResult> results = new ArrayList<ScanResult>();
        int totalStrands = (strand.equals("BOTH"))? 2 : 1;
        for (int strandNum = 0; strandNum < totalStrands; strandNum++) {
            strand = (totalStrands == 2)? strands[strandNum] : strand;
            for (int nmat = 0; nmat < pwms.pwms.length; nmat++) {
                for (int nseq = 0; nseq < S.length; nseq++) {
                    Background BG = bgModels.get(seqLabels[nseq]);
                    ScanRunner run = new ScanRunner(S[nseq], pwms.pwms[nmat], strand, BG, Scores.values[nmat], nucsAfterTSS, pwms.labels[nmat], seqLabels[nseq]);
                    ScanResult result = run.call();
                    if (result.hitLocs.length > 0) results.add(result);
//                    futResults.add(threadPool.submit(run));
                }
            }
        }

//        // Get all the scan results!
//        List<ScanResult> results = new ArrayList<ScanResult>();
//        for (int i = 0; i < futResults.size(); i++) {
//            ScanResult result = null;
//            try {
//                result = futResults.get(i).get(); // get() blocks until the result is available
//
//                // Only keep scans that had results
//                if (result.hitLocs.length > 0) results.add(result);
//            } catch (Exception e) {
//                System.err.println(e.getMessage());
//                e.printStackTrace();
//                break;
//            }
//        }

        // Open file for printing observations
        PrintWriter outFileLocs = null;
        PrintWriter outFileScores = null;
        PrintWriter outFileCumScores = null;
//        PrintWriter outFileTables = null;

//        List<Future<String>> submittedPlots = new ArrayList<Future<String>>();
        
        for (int strandNum = 0; strandNum < totalStrands; strandNum++) {
            strand = (totalStrands == 2)? strands[strandNum] : strand;
//            System.out.println(strand);

            // We're printing out next strand results - close prev strand files and open new ones
            if (strandNum == 1) {
                outFileLocs.close();
                outFileCumScores.close();
//                outFileTables.close();
            }

            outFileLocs = new PrintWriter(new FileWriter(out_Fname + "." + strand + ".locs"));
            outFileCumScores = new PrintWriter(new FileWriter(out_Fname + "." + strand + ".cumscores"));
//            outFileTables = new PrintWriter(new FileWriter(out_Fname + "." + strand + ".table"));

            // Calculate the cumulative score for all sequences scanned by each PWM for the given strand
            Hashtable <String, Hashtable <Integer, Double>> pwmResults = Utils.getCumScore(results, strand);
    
            // Print ROE table header
//            outFileTables.println("MaxPeakLoc\tHalfWidth\tLeft\tRight");
            for (int nmat = 0; nmat < pwms.pwms.length; nmat++) {
                Hashtable <Integer, Double> cshash = null;
    
                if (pwmResults.containsKey(pwms.labels[nmat])) cshash = pwmResults.get(pwms.labels[nmat]);
                
                if (cshash != null) {
                    ArrayList<Integer> keys = new ArrayList<Integer>(cshash.keySet());
                    Collections.sort(keys);
                    ArrayList <Double> values = new ArrayList <Double>();
    
	                // Print locations and cumulative scores
	                outFileLocs.print(pwms.labels[nmat]);
	                outFileCumScores.print(pwms.labels[nmat]);
	                
	                for (int i = 0; i < keys.size(); i++) {
                        Integer loc =  keys.get(i);
                        Double cumscore = cshash.get(loc);
                        values.add(cumscore);
                        outFileLocs.print("\t" + loc );
                        outFileCumScores.print("\t" + Print.df2.format(cumscore.doubleValue()));
    
                    }
//	                if (plotDir != null) {
//	                	ROETable tableFinder = new ROETable(pwms.labels[nmat], strand, keys, values);
//	                	submittedPlots.add(threadPool.submit(tableFinder));
//	                	System.out.println("plot submitted: " + pwms.labels[nmat]);
//	                }
	                outFileLocs.print("\n");
	                outFileCumScores.print("\n");
                }
            }
//            List<String> callBackList = new ArrayList<String>();
//            for (int i = 0; i < submittedPlots.size(); i++) {
//                String result = null;
//                try {
//                    result = submittedPlots.get(i).get(); // get() blocks until the result is available
//                    callBackList.add(result);
//                    System.out.println("returned : " + result);
//                } catch (Exception e) {
//                    System.err.println(e.getMessage());
//                    e.printStackTrace();
//                    break;
//                }
//            }
            // Put pwm labels in order in a file to produce ROEs woth respect to that order
//            String pwmLabelFile = plotDir + "/pwm_labels.txt"; 
//            PrintWriter pwmLabelWriter = new PrintWriter(new FileWriter(pwmLabelFile));
//            pwmLabelWriter.println("pwm");
//            for (String pwmLabel : pwms.labels) {
//            	pwmLabelWriter.println(pwmLabel);
//			}
//            pwmLabelWriter.close();
//           
//            //print out final ROE table
//            printROETable(out_Fname + "." + strand + ".table", strand, pwmLabelFile);
        }

//        outFileTables.close();
        outFileLocs.close();
        outFileCumScores.close();

        System.out.println("Processed all PWMs, shutting down now...");

        pool.shutdown();

        try {
            while(!pool.isTerminated()) {
                pool.awaitTermination(2, TimeUnit.SECONDS);
                pool.shutdownNow();
            }
        }
        catch (InterruptedException ex) {
            System.out.println("Got interrupted: "+ex);
            pool.shutdownNow();
        }
    }
    
    // this method concatinate all the tables that we got for each factor and write out into final ROE table
    // this method should be replaced with another one because the order of factors(PWMs) should be the same as pwm file
    /**
     * @deprecated
     * @param outFileTable
     * @param strand
     */
    private static void printROETable_notSorted(String outFileTable, String strand) {
    	try {
    		String rFile = plotDir + "/" + strand + ".mergeAllTables.R";
            PrintWriter rFH = new PrintWriter(new FileWriter(rFile));
            rFH.println("p <- paste(\"\\\\D*\\\\.\",\"" + strand + "\", \"\\\\.table\", sep=\"\")");
            rFH.println("datadir <- \"" + plotDir + "\"");
            rFH.println("flist <- list.files(datadir, pattern = p)");
            rFH.println("all <- data.frame(MaxPeakLoc=c(), HalfWidth=c(), Left=c(), Right=c())");
            rFH.println("for (f in flist) {");
            rFH.println("rfile <- paste()");
            rFH.println("outfile <- paste(datadir, \"/\", f, sep=\"\")");
            rFH.println("tbl <- read.table(outfile, row.names = 1, header=T, sep = \"\t\")");
            rFH.println("all <- rbind(all, tbl)");
            rFH.println("write.table(all, file=\"" + outFileTable + "\" , row.names=T, col.names=T, quote=F, sep=\"\t\")");
            rFH.println("}");
            
            rFH.close();

//            SysCom cmd = Utils.runSystemCommand("R --quiet --slave -f " + rFile);
            SysCom cmd = Utils.runSystemCommand("R --slave -f " + rFile);

            File rFileTmp = new File(rFile);
            rFileTmp.delete();
            
    	} catch (Exception e) {
    		e.printStackTrace();
		}
    	
	}
    
    private static void printROETable(String outFileTable, String strand, String pwmLabelFile) {
    	try{
    		String rFile = plotDir + "/" + strand + ".mergeAllTables.R";
            PrintWriter rFH = new PrintWriter(new FileWriter(rFile));
        	rFH.println("datadir <- \"" + plotDir + "\"");
        	rFH.println("pwms <- read.table(\"" + pwmLabelFile + "\", header = T)");
        	rFH.println("strand <- \"" + strand + "\"");
        	rFH.println("all <- data.frame(MaxPeakLoc=c(), HalfWidth=c(), Left=c(), Right=c())");
        	rFH.println("for (pwmName in pwms$pwm) {");
        	rFH.println("roe_file <- paste(datadir, \"/\", pwmName, \".\", strand, \".table\", sep = \"\")");
            rFH.println("tbl <- read.table(roe_file, row.names = 1, header=T, sep = \"\\t\")");
            rFH.println("all <- rbind(all, tbl)");
            rFH.println("}");
            rFH.println("write.table(all, file=\"" + outFileTable + "\" , row.names=T, col.names=T, quote=F, sep=\"\\t\")");
            rFH.close();

            SysCom cmd = Utils.runSystemCommand("R --slave -f " + rFile);

            File rFileTmp = new File(rFile);
//            rFileTmp.delete();
    	}catch (Exception e) {
    		e.printStackTrace();
		}
    }

    public static class ROETable implements Callable<String> {
        private String pwmLabel;
        private String strand;
        private ArrayList <Integer> hitLocs;
        private ArrayList <Double> hitScores;

        public ROETable(String pwmLabel, String strand, ArrayList<Integer> hitLocs, ArrayList <Double> hitScores) {
            this.pwmLabel = pwmLabel;
            this.strand = strand;
            this.hitLocs = hitLocs;
            this.hitScores = hitScores;
        }

        public String call() {
            try {
            	// Generate distribution files to use in R for finding ROE Peaks
            	// These dist files can be used in PLOT section too
            	String rFile = plotDir + "/" + pwmLabel + "." + strand + ".computeTable.tmp.R";
            	String outFileName = plotDir + "/" + pwmLabel + "." + strand + ".table";
                //File f = new File(plotDir + "/" + rFile);
                PrintWriter rFH = new PrintWriter(new FileWriter(rFile));

                String strandDir = (strand.equals("FWD"))? "fwd" : "rev";
                String tableFileDir = plotDir + "/" + strandDir;

                // print out calculated variables for scans and regions of enrichment (ROE) to R file
                
                String distFile = plotDir + "/" + pwmLabel + "." + strand + ".dist";
//                System.out.println(distFile);
                PrintWriter distWriter = new PrintWriter(new FileWriter(distFile));
                for (int i = 0; i < hitLocs.size(); i++) {
                	distWriter.println(hitLocs.get(i) + "\t" + hitScores.get(i));
                }
                distWriter.close();
                
                /**
                 * Mitra: This part comes from R code compute_table.R
                 */
                rFH.println("argmax <- function(x, y, w=1, ...) {\n" +
                		"\t require(zoo)\n" +
                		"\t n <- length(y)\n" +
                		"\t y.smooth <- loess(y ~ x, ...)$fitted\n" +
                		"\t y.max <- rollapply(zoo(y.smooth), 2*w+1, max, align=\"center\")\n" +
                		"\t delta <- y.max - y.smooth[-c(1:w, n+1-1:w)]\n" + 
                		"\t i.max <- which(delta <= 0) + w\n" + 
                		"\t list(x=x[i.max], i=i.max, y.hat=y.smooth)\n" +
                		" }");
                
                rFH.println("dist_tbl <- read.table(\""+ distFile + "\", header=F)");
                rFH.println("locs <- dist_tbl[,1] \n" + 
                		    "scores <- dist_tbl[,2] \n" + 
                		    "varuse <- \"\" \n" + 
                		    "maxdist <- length(locs) / 2 \n");
                
                rFH.println("outfname <- \"" + outFileName + "\"");
                rFH.println("outTable <- matrix(nrow=1, ncol=4)");

                String plotFileDir = plotDir + "/" + strandDir;
                rFH.println("## Use the smoothing method to find peaks \n" + 
                			"w <- 1000 \n" + 
                			"span <- 0.05 \n" + 
                			"peaks <- argmax(locs, scores, w=w, span=span)");
                rFH.println("smooth_plotfile <- paste(\"" + plotFileDir + "/\", \"" + pwmLabel + "\", \".smoothed\", \".jpg\", sep=\"\")");
                rFH.println("jpeg(file=smooth_plotfile) \n" + 
                			"plot(locs, scores, cex=0.75, col=\"Gray\", main=paste(\"w = \", w, \", span = \", span, sep=\"\")) \n " + 
                			"lines(locs, peaks$y.hat,  lwd=2) \n" + 
                			"scores.min <- min(scores) \n" + 
                			"sapply(peaks$i, function(i) lines(c(locs[i],locs[i]), c(scores.min, peaks$y.hat[i]),col=\"Red\", lty=2)) \n" + 
                			"points(locs[peaks$i], peaks$y.hat[peaks$i], col=\"Red\", pch=19, cex=1.25) \n" + 
                			"dev.off()\n");
                rFH.println("## process found peaks");
                rFH.println("smooth_scores <- peaks$y.hat");
                
                rFH.println("if (length(peaks$i) == 0) { \n" + 
                			"\t print(\"no peaks were found!\") \n" +
                			"\t varuse = \"NIX\" \n" + 
            				"} else if (length(peaks$i) > 1){ ");
                rFH.println("\t print(\"more than one peak detected! Get the closest one to max value\") \n " +
                			"\t max_idx <- which(peaks$y.hat == max(smooth_scores)) \n" + 
                			"\t mode_idxes <- peaks$i \n" + 
                			"\t distance_to_max <- 10000 \n" + 
                			"\t closest <- \"\"");
                rFH.println("\t for (idx in mode_idxes) { \n" + 
                			"\t\t distance = abs(idx - max_idx) \n" + 
                			"\t\t if (distance < distance_to_max) { \n" + 
                				"\t\t\t distance_to_max <- distance \n" + 
                				"\t\t\t closest <- idx \n" + 
            				"\t\t } \n" + 
            				"\t} \n" +
            				"\t peak_mod_index <- closest \n" + 
            				"\t peak_mode_loc <- locs[peak_mod_index] \n" +
            				"} else { ");
                rFH.println("\t peak_mod_index <- peaks$i[1] \n" + 
                			"\t peak_mode_loc <- locs[peak_mod_index] \n" + 
            				"}\n");
                
//                rFH.println("maxbgval <- mean(density(smooth_scores)$x)\n");
                
                rFH.println("#### Use old method to find region but using smoothed values");
                rFH.println("if ((peak_mode_loc < -maxdist) || (peak_mode_loc > maxdist) || varuse == \"NIX\") { \n " +
                			"\t varuse <- \"NIX\" \n" + 
            				"} else { \n" + 
                			"\t maxbgval <- mean(density(smooth_scores[1:peak_mod_index])$x) \n" +
            				"\t buffer <- 2 \n" + 
            				"\t # Find Left \n " + 
            				"\t leftind <- peak_mod_index \n" +
            				"\t nBelow <- 0 \n " +
            				"\t while (nBelow < buffer) { \n" + 
            				"\t\t if (leftind > length(scores) | leftind < 1) { \n" + 
            				"\t\t\t varuse <- \"NIX\" \n" + 
            				"\t\t\t break \n " + 
            				"\t\t } \n" + 
            				"\t\t if (smooth_scores[leftind] < maxbgval) { \n" + 
            				"\t\t\t nBelow <- nBelow + 1 \n" + 
            				"\t\t } \n" + 
            				"\t\t leftind <- leftind - 1 \n" + 
            				"\t } " ); 
                
                rFH.println("\t # Find Right");
                rFH.println("\t rightind <- peak_mod_index");
                rFH.println("\t maxbgval <- mean(density(smooth_scores[peak_mod_index:length(smooth_scores)])$x)");
                rFH.println("\t nBelow <- 0");
                rFH.println("\t while (nBelow < buffer) {");
                rFH.println("\t\t if (rightind > length(scores)) {");
                rFH.println("\t\t\t varuse <- \"NIX\"");
                rFH.println("\t\t break");
                rFH.println("\t }");

                rFH.println("\t\t if (smooth_scores[rightind] < maxbgval) { \n" + 
                			"\t\t\t nBelow <- nBelow + 1 \n" + 
            				"\t\t} \n" + 
            				"\t\t rightind <- rightind + 1 \n" + 
            				"\t} \n" +
            				"\t if ((leftind == peak_mod_index) || (rightind == peak_mod_index)) { \n" + 
            				"\t\t varuse <- \"NIX\" \n" + 
            				"\t } \n" +
            				"} \n");

                rFH.println("");
                
                //------------------
                rFH.println("if (varuse == \"NIX\") {");
                rFH.println("\t outTable[1,1] <- NA");
                rFH.println("\t outTable[1,2] <- NA");
                rFH.println("\t outTable[1,3] <- NA");
                rFH.println("\t outTable[1,4] <- NA");
                rFH.println(" } else {");
                rFH.println("\t left <- locs[leftind]");
                rFH.println("\t right <- locs[rightind]");
                rFH.println("\t halfWidth <- (right - left)/2.0;");
                rFH.println("\t if (halfWidth > 250) {");
                rFH.println("\t\t halfWidth <- 250 \n" +
                			"\t\t left <- peak_mode_loc - halfWidth \n " + 
                			"\t\t right <- peak_mode_loc + halfWidth \n" + 
                			"\t } \n" +
                			"\t outTable[1,1] <- format(peak_mode_loc, digits=2) \n" + 
                			"\t outTable[1,2] <- format(halfWidth, digits=2) \n" + 
                			"\t outTable[1,3] <- format(left, digits=2) \n" + 
                			"\t outTable[1,4] <- format(right, digits=2) \n" +
                			"} \n");
                rFH.println("");
                
                rFH.println("# Print table of values");
                rFH.println("dimnames(outTable)[[2]] <- c(\"MaxPeakLoc\", \"HalfWidth\", \"Left\", \"Right\")");
                rFH.println("dimnames(outTable)[[1]] <- \"" + pwmLabel + "\"");
                rFH.println("write.table(outTable, file=outfname, row.names=T, col.names=T, quote=F, sep=\"\t\")");
                
                rFH.println("# Plot ROEs");
                // Setup Title info
                rFH.println("if (varuse != \"NIX\") {");
                rFH.println("subT <- paste(\"Interval = [\", toString(format(left, digits=2)), \",\", toString(format(right, digits=2)), \"]\")");
                rFH.println("mainT <- paste(\"" + pwmLabel + "\", subT, sep=\"\\n\")");

                // Set the file name / type for the plot output
                
                rFH.println("plot_filename <- paste(\"" + plotFileDir + "/\", \"" + pwmLabel + "\", \".pk\", \".jpg\", sep=\"\")");
                rFH.println("jpeg(file=plot_filename)");

                // Plot out the data
                rFH.println("par(mfrow=c(2,1))");
                rFH.println("plot(locs, scores, main=mainT, xlab=\"locs\", ylab=\"scores\")");
                rFH.println("points(locs[locs >= left & locs <= right], scores[locs >= left & locs <= right], col='red', pch=1)");
                rFH.println("plot(locs[(locs >= left - 20) & (locs <= right + 20)], scores[(locs >= left - 20) & (locs <= right + 20)], main=mainT, xlab=\"locs\", ylab=\"scores\")");
                rFH.println("points(locs[locs >= left & locs <= right], scores[locs >= left & locs <= right], col='red', pch=1)");
                rFH.println("points(peak_mode_loc, maxbgval, pch=20, col='green')");
                rFH.println("arrows(left, maxbgval, right, maxbgval, col='green', length=0.1, code=3)");
                rFH.println("results <- dev.off()");
                rFH.println("}");

                rFH.close();

                // Could add error checking if you want, but this command is pretty safe
                SysCom cmd = Utils.runSystemCommand("R --slave -f " + rFile);

                File rFileTmp = new File(rFile);
//                rFileTmp.delete();
            }
            catch (java.io.IOException e) {
                e.printStackTrace();
            }
            return pwmLabel + "," + strand + " - done";
        }
        
    }
    

    public static void setupPlotDir() {
        try {
            File f = new File(plotDir);

            // Make sure directory exists
            if (f.exists() && !f.isDirectory()) {
                System.err.println("Error: file '" + plotDir + "' exists and isn't a directory.  Change directory name or remove existsing file.");
                System.exit(0);
            } else if (!f.exists()) {
                f.mkdir();
            }

            if (strand.equals("BOTH")) {
                f = new File(plotDir + "/fwd");
                // Make sure directory exists
                if (f.exists() && !f.isDirectory()) {
                    System.err.println("Error: file '" + plotDir + "/fwd' exists and isn't a directory.  Change directory name or remove existsing file.");
                    System.exit(0);
                } else if (!f.exists()) {
                    f.mkdir();
                }

                f = new File(plotDir + "/rev");
                // Make sure directory exists
                if (f.exists() && !f.isDirectory()) {
                    System.err.println("Error: file '" + plotDir + "/rev' exists and isn't a directory.  Change directory name or remove existsing file.");
                    System.exit(0);
                } else if (!f.exists()) {
                    f.mkdir();
                }
            } else {
                String strandDir = (strand.equals("FWD"))? "fwd" : "rev";
                f = new File(plotDir + "/" + strandDir);
                // Make sure directory exists
                if (f.exists() && !f.isDirectory()) {
                    System.err.println("Error: file '" + plotDir + "/" + strandDir + "' exists and isn't a directory.  Change directory name or remove existsing file.");
                    System.exit(0);
                } else if (!f.exists()) {
                    f.mkdir();
                }
            }
        }
        catch (Exception e) {
            e.printStackTrace();
            System.exit(0);
        }
    }

    private static void parseArgs(String[] args) {
        // Setup CLI to simplify Scan interface;
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
                if (args.length < 3) {
                    throw new ParseException("Must supply a PWM File, FASTA file, and output basename.");
                } else {
                    pwms_Fname = args[0];
                    promoterSeqs_Fname = args[1];
                    out_Fname = args[2];
                }

                if (cmdLine.hasOption("nprocs")) {
                   threadPoolSize = Integer.parseInt(cmdLine.getOptionValue("nprocs"));
                   if (threadPoolSize > Runtime.getRuntime().availableProcessors()) {
                       System.err.println("Warning: your system does not have " + threadPoolSize + " CPUs: setting nprocs to " + Runtime.getRuntime().availableProcessors() + ".");
                  
                       threadPoolSize = Runtime.getRuntime().availableProcessors();
                   }
                   if (threadPoolSize < 1) {
                       System.err.println("Warning: nprocs must be a value >= 1. Setting nprocs = 1");
                       threadPoolSize = 1;
                   }
                }

                if (cmdLine.hasOption("BGWIN")) {
                    BG_WIN = Integer.parseInt(cmdLine.getOptionValue("BGWIN"));
                }

                if (cmdLine.hasOption("minScores")) {
                    scoreCutoffs_Fname = cmdLine.getOptionValue("minScores");
                }
                
                if (cmdLine.hasOption("strand")) {
                    strand = cmdLine.getOptionValue("strand");
                    if (!strand.equals("FWD") && !strand.equals("REV") && !strand.equals("BOTH")) {
                        strand = "BOTH";
                    }
                }

                if (cmdLine.hasOption("nucsAfterTSS")) {
                    nucsAfterTSS = Integer.parseInt(cmdLine.getOptionValue("nucsAfterTSS"));
                }

                if (cmdLine.hasOption("pseudoCounts")) {
                    pseudoCountsVal = Double.parseDouble(cmdLine.getOptionValue("pseudoCounts"));
                }

                if (cmdLine.hasOption("seqName")) {
                    seqName = cmdLine.getOptionValue("seqName");
                }

                if (cmdLine.hasOption("plotDir")) {
                    plotDir = cmdLine.getOptionValue("plotDir");
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
            help = true;
        }

        if (help) {
            OptionComparator opt_order = new OptionComparator();
            String cmdLineSyntax = "java -jar tfbs_llscan.jar ROEFinder <PWM FILE> <FASTA FILE> <FILEBASE>";
            String header = "\nJava tool for generating regions of enrichment for loglikelihood scans of FASTA sequences. Note: All matrices in <PWM_FILE> are assumed to be frequency matrices - not probabilities\n\nOPTIONS\n\n";
            String footer = "";
            HelpFormatter formatter = new HelpFormatter();
            formatter.setOptionComparator(opt_order);
            formatter.setWidth(80);
            formatter.printHelp(cmdLineSyntax, header, options, footer, true);
            System.exit(0);
        }
    }

    private static class OptionComparator <T extends Option> implements Comparator <T> {
        private ArrayList <String> OPTS_ORDER = new ArrayList <String> ();

        public OptionComparator () {
            OPTS_ORDER.add("help");
            OPTS_ORDER.add("nprocs");
            OPTS_ORDER.add("BGWIN");
            OPTS_ORDER.add("strand");
            OPTS_ORDER.add("minScores");
            OPTS_ORDER.add("pseudoCounts");
            OPTS_ORDER.add("nucsAfterTSS");
            OPTS_ORDER.add("plotDir");
            OPTS_ORDER.add("seqName");
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

        // sequence to limit scans to
        options.addOption(Option.builder()
            .longOpt("seqName")
            .hasArg(true)
            .argName("SEQNAME")
            .desc("Limit ROE scans to sequence <SEQNAME>.")
            .required(false)
            .build());

        // the offset to adjust the TSS position in the input sequence
        options.addOption(Option.builder()
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

        // the offset used to adjust loglikelihood site positions relative to TSS
        options.addOption(Option.builder("n")
            .longOpt("nprocs")
            .hasArg(true)
            .argName("INT")
            .desc("Number of processors to use (default: 1)")
            .required(false)
            .type(Integer.class)
            .build());
            
        // length of background window used for local background sequence calculations
        options.addOption(Option.builder("B")
            .longOpt("BGWIN")
            .hasArg(true)
            .argName("INT")
            .desc("Background Window size used for calculating local sequence background, set this to '0' if you want to use a neutral distribution (i.e. each base is equally likely) (default: 250)")
            .required(false)
            .build());

        // strand to scan
        options.addOption(Option.builder("s")
            .longOpt("strand")
            .argName("FWD|REV|BOTH")
            .hasArg(true)
            .desc("the strand of the input sequence to be scanned (default: BOTH)")
            .required(false)
            .build());

        // score threshold file
        options.addOption(Option.builder()
            .longOpt("minScores")
            .hasArg(true)
            .argName("FILENAME")
            .desc("file to read mininum threshold scores for scans to be add (default: all thresholds set to 0)")
            .required(false)
            .build());

        // score threshold file
        options.addOption(Option.builder()
            .longOpt("plotDir")
            .hasArg(true)
            .argName("DIRNAME")
            .desc("directory to store ROE plots for each PWM (default: do not generate plots)")
            .required(false)
            .build());

        return options;
    }
}
