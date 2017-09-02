package osu.megraw.llscan;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;
import java.util.concurrent.Future;
import java.util.concurrent.LinkedBlockingQueue;
import java.util.concurrent.ThreadPoolExecutor;
import java.util.concurrent.TimeUnit;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


/**
 * This class generates features based on tiled array windowing around the TSS
 * @author mitra
 *
 */
public class GenFeaturesTiledWindows {
	String[] strands = {"FWD", "REV"};
	public char[][] sequenceCharArr; // holds input sequences from input FASTA file
	public DoubleColReturn scoreCutOffs;
	public int nucsDownStream = -1;
	public String scoreCutoffs_Fname = null;
    public int BG_WIN = 250;
    public int nproc = 1;
    public double PseudoCountsVal = 0.01;
    public String map_Fname = null;
    public String pwms_Fname;
    public String out_Fname;
    
    public boolean USE_MARKOV_1_BG = true;
	private boolean help = false;
	private Integer tilingNucsUpsream;
	private Integer tilingNucsDownstream;
	private Integer winsWidth;
	private String inputSeqName;
	String[] seqLabels;
	private int seqLength = -1;
	HashMap<String, LefRightPair<Integer, Integer>> windowsHash = new HashMap<>();
	
	HashMap <String, RefSeqData> fastaCoords = new HashMap <String, RefSeqData>();
	PWMReturn pwms = null;

    public static double[] bgdefault = {0.25, 0.25, 0.25, 0.25};
    public static double[][] bgM1default = {{0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25}};

    
	public static void main(String[] args) {
		try {
			GenFeaturesTiledWindows genFeaturesTiledWindows = new GenFeaturesTiledWindows();
			genFeaturesTiledWindows.parseArgs(args);
			genFeaturesTiledWindows.readPromoterSeq();
			genFeaturesTiledWindows.readPWMsInfo();
			genFeaturesTiledWindows.createWins();
			genFeaturesTiledWindows.getWinScores();
			System.out.println("DONE!");
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	// This is main method to compute log-likelihood scores within each window for each sequence in parallel
	private void getWinScores() throws Exception {
//		ThreadPoolExecutor threadPool = new ThreadPoolExecutor(nproc, nproc, 4, TimeUnit.SECONDS, new LinkedBlockingQueue());
//        List<Future<LoglikScoreResult>> results = new ArrayList<Future<LoglikScoreResult>>();
//        ArrayList<String> futureSeqNames = new ArrayList<String>();
        
     // file to output mapped genomic coordinates to
        PrintWriter mapOutFile = null;
        if (map_Fname != null) {
            mapOutFile = new PrintWriter(new FileWriter(map_Fname));
            // header for mapping file
            mapOutFile.println("seq_id\tfeature_id\trelative_start\trelative_end\ttss_start\tabsolute_start\tabsolute_end\tref_seq\tgenomic_start\tgenomic_end");
        }
        
        List<LoglikScoreResult> final_results = new ArrayList<>();
        // Iterate over sequences 
        for (int i = 0; i < sequenceCharArr.length; i++) {
        	String seqName = seqLabels[i];
        	ComputeLoglikScoreThread computeLoglikScoreThread = new ComputeLoglikScoreThread(sequenceCharArr[i], scoreCutOffs, nucsDownStream, BG_WIN, seqName, pwms, windowsHash, winsWidth);
        	LoglikScoreResult loglikScoreResult = computeLoglikScoreThread.call();
        	final_results.add(loglikScoreResult);
//        	futureSeqNames.add(seqName);
//        	try {
//				results.add(threadPool.submit(computeLoglikScoreThread));
//                futureSeqNames.add(seqName);
//            }
//            catch (Exception ex) {
//                System.err.println("Error processing example "+seqName+": "+ex.getMessage()+". Quitting.");
//                threadPool.shutdownNow();
//                while(!threadPool.isTerminated()) {
//                    threadPool.awaitTermination(2, TimeUnit.SECONDS);
//                }
//                System.exit(1);
//            }
        	
            if (map_Fname != null) {
                RefSeqData seqInfo = fastaCoords.get(seqName);
                if (seqInfo == null) {
                    System.out.println(seqName);
                    return;
                }

                for (String strand : strands) {
                	for (String pwmID : pwms.labels) {
                		int pwmIdx = pwms.labelIndex.get(pwmID);
                		int w = pwms.pwms[pwmIdx].length;
                		for (int winNo = 1; winNo <= windowsHash.size(); winNo++) {
                			int left = windowsHash.get(String.valueOf(winNo)).getL();
                			int right = windowsHash.get(String.valueOf(winNo)).getR();        	        	
        					int tssLoc = sequenceCharArr[i].length - nucsDownStream;
        					int winL = left + tssLoc;
                            int winR = right + tssLoc;
                            if (winL < 0) { winL = 0; }
                            if (winR > sequenceCharArr[i].length - w) { winR = sequenceCharArr[i].length - w; }
                            
                            if (seqInfo.strand.equalsIgnoreCase("-")) {
                            	int upstreamLen = sequenceCharArr[i].length - nucsDownStream - 1;
                            	tssLoc = sequenceCharArr[i].length - upstreamLen;
                            	
                            	winL = (-1 * windowsHash.get(String.valueOf(winNo)).getR()) + tssLoc; 
                            	winR = (-1 * windowsHash.get(String.valueOf(winNo)).getL()) + tssLoc;
                            	if (winL < 0) { winL = 0; }
                            	if (winR > sequenceCharArr[i].length - w) { winR = sequenceCharArr[i].length - w; }
                            }
                            int genomic_start = seqInfo.start + winL;
                            int genomic_end = seqInfo.start + winR;
                            
                            String featureId = pwmID + "_" + strand + "_" + winNo + "_tile" + winsWidth;

                            mapOutFile.println(seqName + "_0" + "\t" + featureId + "\t" + left + "\t" + right + "\t" + 
                                               tssLoc + "\t" + winL + "\t" + winR + "\t" + seqInfo.id + "\t" + genomic_start + "\t" + genomic_end);
                            
        				}
        	        }
                } 
            }
        }
        if (map_Fname != null) mapOutFile.close();
        
        process_results(final_results);
//        process_results(results, threadPool, futureSeqNames);
        
//        threadPool.shutdownNow();
//        while(!threadPool.isTerminated()) {
//            threadPool.awaitTermination(2, TimeUnit.SECONDS);
//        }
        System.out.println("Complete");

    }
      


	private void process_results(List<LoglikScoreResult> finalResults) throws Exception {
        int numFeatures = -1; // Number of columns printed out - queried from the results array length found below...

//		ArrayList <LoglikScoreResult> finalResults = new ArrayList <LoglikScoreResult>();
		String[] seqContentNames = {"GCcontent", "CAcontent", "GAcontent"};
	    int currentSeq = 0; // If error is thrown, keep track of sequence we were on so error prints out
	       						// the actual sequence being worked on
    	// Write out features
        PrintWriter outFileVars = new PrintWriter(new FileWriter(out_Fname));
        
        // Print header
        for (String strand : strands) {
        	for (String pwmID : pwms.labels) {
    	        for (int winNo = 1; winNo <= windowsHash.size(); winNo++) {
    	        	String featureId = pwmID + "_" + strand + "_" + winNo + "_tile" + winsWidth;
    	        	outFileVars.print("\t" + featureId);
    	        }
        	}	        	
        }
//        outFileVars.print("\t" + "GCcontent\tCAcontent\tGAcontent");
        outFileVars.println();
        
        // Print results
        for (int i = 0; i < finalResults.size(); i++) {
            LoglikScoreResult res = finalResults.get(i);
            outFileVars.print(res.sampleName);
	        for (String strand : strands) {
	        	for (String pwmID : pwms.labels) {
	    	        for (int winNo = 1; winNo <= windowsHash.size(); winNo++) {
	    	        	String featureId = pwmID + "_" + strand + "_" + winNo + "_tile" + winsWidth;
	    	        	outFileVars.print("\t" + res.featureHash.get(featureId));
	    	        }
	        	}	        	
	        }
            outFileVars.println();
        }
        outFileVars.flush();
        outFileVars.close();

	}
	private void process_results(List<Future<LoglikScoreResult>> results, ThreadPoolExecutor threadPool, 
			ArrayList<String> futureSeqNames) throws Exception {
		int numFeatures = -1; // Number of columns printed out - queried from the results array length found below...
		
		ArrayList <LoglikScoreResult> finalResults = new ArrayList <LoglikScoreResult>();
		String[] seqContentNames = {"GCcontent", "CAcontent", "GAcontent"};
		int currentSeq = 0; // If error is thrown, keep track of sequence we were on so error prints out
		// the actual sequence being worked on
		for (int j = 0; j < results.size(); j++) {
			Future<LoglikScoreResult> result = results.get(j);
			
			LoglikScoreResult res = null;
			try {
				res = result.get(); // get() blocks until the result is available
				
				finalResults.add(res);
				
//	    		numberNonzeros += res.numberNonzeros;
				
				// Set the numFeatures equal to number of values obtained from results - all scans should have the
				// *same* number of features generated, so only need to set this once
				if (numFeatures == -1) numFeatures = res.featureHash.size();
				
				currentSeq += 1;
			} catch (Exception ex) {
				String seqName = futureSeqNames.get(currentSeq);
				System.err.println("Error processing example "+seqName+": "+ex.getMessage()+". Quitting.");
				threadPool.shutdownNow();
				while(!threadPool.isTerminated()) {
					threadPool.awaitTermination(2, TimeUnit.SECONDS);
				}
				System.exit(1);
			}
		}
		
		// Write out features
		PrintWriter outFileVars = new PrintWriter(new FileWriter(out_Fname));
		
		// Print header
		for (String strand : strands) {
			for (String pwmID : pwms.labels) {
				for (int winNo = 1; winNo <= windowsHash.size(); winNo++) {
					String featureId = pwmID + "_" + strand + "_" + winNo + "_tile" + winsWidth;
					outFileVars.print("\t" + featureId);
				}
			}	        	
		}
//        outFileVars.print("\t" + "GCcontent\tCAcontent\tGAcontent");
		outFileVars.println();
		
		// Print results
		for (int i = 0; i < finalResults.size(); i++) {
			LoglikScoreResult res = finalResults.get(i);
			outFileVars.print(res.sampleName);
			for (String strand : strands) {
				for (String pwmID : pwms.labels) {
					for (int winNo = 1; winNo <= windowsHash.size(); winNo++) {
						String featureId = pwmID + "_" + strand + "_" + winNo + "_tile" + winsWidth;
						outFileVars.print("\t" + res.featureHash.get(featureId));
					}
				}	        	
			}
			outFileVars.println();
		}
		outFileVars.flush();
		outFileVars.close();
		
	}

	private void createWins() {
		// Compute the windows coordinates relative to TSS mode and store them in windowHash
		int winCount = (tilingNucsUpsream + tilingNucsDownstream) / winsWidth;
		
		int currentWinL = -1 * tilingNucsUpsream;
		for (int winNo = 1; winNo <= winCount; winNo++) {
			int left = currentWinL;
			int right = currentWinL + winsWidth;
			
			LefRightPair<Integer, Integer> winCoordinate = new LefRightPair<Integer, Integer>(left, right);
			windowsHash.put(String.valueOf(winNo), winCoordinate);
			
			currentWinL = right;
		}
	}


	private void readPWMsInfo() throws Exception {
        // Read PWM file, store rev comp PWMs for scanning opposite strand
        pwms = Load.loadPWMFileSimpleHeader(pwms_Fname, PseudoCountsVal);
        pwms.ComputeRevCmpPWMs();

        // Read score threshold file corresponding to PWM file
        if (scoreCutoffs_Fname != null) {
            scoreCutOffs = Load.loadDoubleColumnFile(new File(scoreCutoffs_Fname));
        } else {
            // No threshold was given - assume that threshold will be zero
            scoreCutOffs = new DoubleColReturn(pwms.labels);
        }		
	}


	/** Read promoter sequences and load required variables */
	private void readPromoterSeq() throws IOException {
		FastaReturn Seqs = Load.loadFastaFile(new File(inputSeqName));
		sequenceCharArr = new char[Seqs.lines.length][];
        for (int i = 0; i < Seqs.lines.length; i++) {
            sequenceCharArr[i] = Seqs.lines[i].toCharArray();

            if (seqLength == -1) seqLength = sequenceCharArr[i].length;

            if (seqLength != sequenceCharArr[i].length) {
                System.err.println("Error: input FASTA sequences of different length found with (" + Seqs.headers[i] + ")");
                return;
            }
        }
        
        // no TSS offset given - autoset it to midpoint of input sequence length
        if (nucsDownStream == -1) {
            nucsDownStream = seqLength / 2;  // integer arithmetic - automatically floors
        }
        
        // Store FASTA seq's coordinates in hashmap
        seqLabels = new String[Seqs.headers.length];
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
                    String ref = genomic_parts[1]; // Chromosome
                    int start = Integer.parseInt(genomic_parts[2]); //promoter start location
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
	}
	
	private void parseArgs(String[] args) {
		Options options = buildOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmdLine;
        
        try {
            cmdLine = parser.parse(options, args);

            // Remove arguments that were parsed
            args = new String[cmdLine.getArgList().size()];
            args = cmdLine.getArgList().toArray(args);

            if (cmdLine.hasOption("help")) {
                this.help  = true;
            } else {
                if (args.length < 3) {
                    throw new ParseException("Please supply nucsUpstream, nucsDownStream, window width, FASTA file, PWM file and output file name.");
                } else {
                    this.inputSeqName = args[0];
                    pwms_Fname = args[1];
                    out_Fname = args[2];
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
                
                if (cmdLine.hasOption("nucsUp")) {
                	tilingNucsUpsream = Integer.parseInt(cmdLine.getOptionValue("nucsUp"));
                }

                if (cmdLine.hasOption("nucsDown")) {
                	tilingNucsDownstream = Integer.parseInt(cmdLine.getOptionValue("nucsDown"));
                }
                
                if (cmdLine.hasOption("winWidth")) {
                	winsWidth = Integer.parseInt(cmdLine.getOptionValue("winWidth"));
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
                
//                if (cmdLine.hasOption("pos")) {
//                    this.posArrString = cmdLine.getOptionValue("pos");
//                }

//                if (cmdLine.hasOption("flankingWindows")) {
//                    String[] flankingWindowParts = cmdLine.getOptionValue("flankingWindows").split(" ");
//
//                    if (flankingWindowParts.length != 2) {
//                        throw new ParseException("Incorrect number of flanking window arguments: \"\"|\"FLANK_LENGTH NUM_WIN\")");
//                    }
//        
//                    flankingFeatureWindowLength = Integer.parseInt(flankingWindowParts[0]);
//                    numFlankingWindows = Integer.parseInt(flankingWindowParts[1]);
//        
//                    if (flankingFeatureWindowLength > 0) {
//                        if (flankingFeatureWindowLength % numFlankingWindows != 0) {
//                            throw new ParseException("Flanking window length not divisible by number of windows.");
//                        }
//                        else if (numFlankingWindows > flankingFeatureWindowLength) {
//                            throw new ParseException("More flanking windows than the flanking region length.");
//                        }
//                    }
//                }

//                if (cmdLine.hasOption("bin")) {
//                    this.outputBin = true;
//                }
//
//                if (cmdLine.hasOption("byNT")) {
//                    this.byNT = true;
//                    this.flankingFeatureWindowLength = 0;
//                    this.numFlankingWindows = 0;
//                }
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
            OptionComparator<Option> opt_order = new OptionComparator<>();
            String cmdLineSyntax = "java -jar GenFeaturesTiledWindows tfbs_llscan.jar <upstream nucs> <downstream nucs> <window width> <sequence.fasta> <PWM_FILE> <Outfile.rdat>";
            String header = "\nJava tool for generating features within non-overlapping windows of given width from -nucsUpstream to +nucsDownstream of TSS  for loglikelihood " + 
                            "scans of FASTA sequences.\n\nOPTIONS\n\n";
            String footer = "";
            HelpFormatter formatter = new HelpFormatter();
            formatter.setOptionComparator(opt_order);
            formatter.setWidth(120);
            formatter.printHelp(cmdLineSyntax, header, options, footer, true);
            System.exit(0);
        }

	}

	public GenFeaturesTiledWindows() {
	}
	
    private Options buildOptions() {
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
//        options.addOption(Option.builder()
//            .longOpt("pos")
//            .hasArg(true)
//            .argName("Position String")
//            .desc("Position String: String describing where to scan from. One of:\n" +
//                  "\"N\": A number indicating the location to sample from (relative to TSS of sequence)\n\n" +
//                  "\"Rand [ndraws] [draw_start] [draw_stop]\": Choose [ndraws] random locations between [draw_start] and [draw_stop])\n\n" +
//			      "\"Range [i] [range_start] [range_stop]\": Scan every [i] nucleotides between [range_start] and [range_stop])\n" + 
//			      "\"Replicate [filename]\": Scan every nucleotide defined in [filename])\n" + 
//                  "(default: 0 - i.e., only generate features from TSS of sequence)")
//            .required(false)
//            .build());

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
            .desc("Total pseudo-counts to add to PWM matrix prior to normalization (default: 0.01)")
            .required(false)
            .type(Double.class)
            .build());
        
        // How many upstream nts take for windowing promoter?
        options.addOption(Option.builder()
        		.longOpt("nucsUp")
        		.hasArg(true)
        		.argName("INT")
        		.desc("length of upstream region of TSS for dividing into non-overlapping windows (for example 1000 nucleotide upstream of TSS)")
        		.required(true)
        		.type(Double.class)
        		.build());
        
        // How many downstream nts take for windowing promoter?
        options.addOption(Option.builder()
        		.longOpt("nucsDown")
        		.hasArg(true)
        		.argName("INT")
        		.desc("length of downstream region of TSS for dividing into non-overlapping windows (for example 100 nucleotide downstream of TSS)")
        		.required(true)
        		.type(Double.class)
        		.build());
        
        // Width of of each window
        options.addOption(Option.builder()
        		.longOpt("winWidth")
        		.hasArg(true)
        		.argName("INT")
        		.desc("width of non-overlapping windows")
        		.required(true)
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

//        // flanking windows string for generating features from flanking sequence
//        options.addOption(Option.builder()
//            .longOpt("flankingWindows")
//            .hasArg(true)
//            .argName("WindowLength NumWindows")
//            .desc("If specified, generates additional features for each TSS centered upstream and downstream. " +
//                  "Generates NumWindows additional features on each side of TSS, going up to WindowLength away from site." +
//                  "(default: no flanking windows are used)")
//            .required(false)
//            .build());
//
//        // Output data in binary format
//        options.addOption(Option.builder("b")
//            .longOpt("bin")
//            .hasArg(false)
//            .desc("Output data in binary format.")
//            .required(false)
//            .build());
//
//        // Output scores on a per-nucleotide basis for each ROE scanned
//        options.addOption(Option.builder("y")
//            .longOpt("byNT")
//            .hasArg(false)
//            .desc("Output scores on a per-nucleotide basis for each ROE scanned.  Turns off reporting of flanking features and sequence content features (i.e. GC/GA/CA content)")
//            .required(false)
//            .build());

        return options;
    }

    private static class OptionComparator <T extends Option> implements Comparator <T> {
        private ArrayList <String> OPTS_ORDER = new ArrayList <String> ();

        public OptionComparator () {
            OPTS_ORDER.add("help");
            OPTS_ORDER.add("pos");
            OPTS_ORDER.add("nucsUp");
            OPTS_ORDER.add("nucsDown");
            OPTS_ORDER.add("winWidth");
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
}
