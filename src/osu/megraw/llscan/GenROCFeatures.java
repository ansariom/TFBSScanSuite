package osu.megraw.llscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;


/**
 * Start with open regions in each promoter and compute cumulitive log-lik scores for each region-promoter
 * Now the features are going to be dependent on open regions
 * @author mitra
 *
 */
public class GenROCFeatures {
	String[] strands = {"FWD", "REV"};
	public char[][] sequenceCharArr; // holds input sequences from input FASTA file
	public DoubleColReturn scoreCutOffs;
	public int nucsAfterTSS = -1;
	public String scoreCutoffs_Fname = null;
    public int BG_WIN = 250;
    public int nproc = 1;
    public double PseudoCountsVal = 0.001;
    public String map_Fname = null;
    public String pwms_Fname;
    public String out_Fname;
	private Integer tilingNucsUpsream;
	private Integer tilingNucsDownstream;
	private Integer winsWidth;

    
    public boolean USE_MARKOV_1_BG = true;
	private boolean help = false;
	private String inputSeqName;
	String[] seqLabels;
	private int seqLength = -1;
	HashMap<String, LefRightPair<Integer, Integer>> windowsHash = new HashMap<>();
	
	HashMap <String, RefSeqData> fastaCoords = new HashMap <String, RefSeqData>();
	PWMReturn pwms = null;
	private String rootOcRegionsFile;
	private String leafOcRegionsFile;
	private HashMap<String, List<Coordinate>> leafOCHash = new HashMap<>();
	private HashMap<String, List<Coordinate>> rootOCHash = new HashMap<>();

    public static double[] bgdefault = {0.25, 0.25, 0.25, 0.25};
    public static double[][] bgM1default = {{0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25}};

	public static void main(String[] args) {
		try {
			GenROCFeatures genROCFeatures = new GenROCFeatures();
			genROCFeatures.parseArgs(args);
			genROCFeatures.readPromoterSeq();
			genROCFeatures.readPWMsInfo();
			genROCFeatures.readOCRegions();
			genROCFeatures.createWins();
			genROCFeatures.getOCScores();
			System.out.println("DONE!");
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	
	private void readOCRegions() {
		// Read root and lead oc regions separately
		try {
			BufferedReader leafFile = new BufferedReader(new FileReader(leafOcRegionsFile));
			BufferedReader rootFile = new BufferedReader(new FileReader(rootOcRegionsFile));
			
			String line = null;
			while ((line = leafFile.readLine()) != null) {
				line.trim();
				String[] parts = line.split("\t");
				String tssName = parts[0];
				String chrom = parts[2];
				int relLeft = Integer.valueOf(parts[5]);
				int relRight = Integer.valueOf(parts[6]);
				
				Coordinate coordinate = new Coordinate(chrom, relLeft, relRight);
				List<Coordinate> coordinateList = leafOCHash.get(tssName);
				if (coordinateList == null) 
					coordinateList = new ArrayList<>();
				coordinateList.add(coordinate);
				leafOCHash.put(tssName, coordinateList);
			}
			
			leafFile.close();
			
			line = null;
			while ((line = rootFile.readLine()) != null) {
				line.trim();
				String[] parts = line.split("\t");
				String tssName = parts[0];
				String chrom = parts[2];
				int relLeft = Integer.valueOf(parts[5]);
				int relRight = Integer.valueOf(parts[6]);
				
				Coordinate coordinate = new Coordinate(chrom, relLeft, relRight);
				List<Coordinate> coordinateList = rootOCHash.get(tssName);
				if (coordinateList == null) 
					coordinateList = new ArrayList<>();
				coordinateList.add(coordinate);
				rootOCHash.put(tssName, coordinateList);
			}
			
			rootFile.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
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

	private void getOCScores() throws Exception {
		List<LoglikScoreResult> final_results = new ArrayList<>();
        // Iterate over sequences 
        for (int i = 0; i < sequenceCharArr.length; i++) {
        	String seqName = seqLabels[i];
        	List<Coordinate> leafCoords = leafOCHash.get(seqName);
        	List<Coordinate> rootCoords = rootOCHash.get(seqName);
        	ComputeLoglikScoreROC computeLoglikScoreThread = new ComputeLoglikScoreROC(sequenceCharArr[i], scoreCutOffs, nucsAfterTSS, BG_WIN, seqName, pwms, windowsHash, winsWidth, leafCoords, rootCoords, tilingNucsUpsream);
        	LoglikScoreResult loglikScoreResult = computeLoglikScoreThread.call();
        	final_results.add(loglikScoreResult);
        }
        process_results(final_results);
	}
	
	private void process_results(List<LoglikScoreResult> finalResults) throws Exception {
    	// Write out features
        PrintWriter outFileVars = new PrintWriter(new FileWriter(out_Fname));
        
        // Print header
        for (String strand : strands) {
        	for (String pwmID : pwms.labels) {
    	        for (int winNo = 1; winNo <= windowsHash.size(); winNo++) {
    	        	String featureId = pwmID + "_" + strand + "_" + winNo + "_ROC_LEAF_" + winsWidth;
    	        	outFileVars.print("\t" + featureId);
    	        	featureId = pwmID + "_" + strand + "_" + winNo + "_ROC_ROOT_" + winsWidth;
    	        	outFileVars.print("\t" + featureId);
    	        }
        	}	        	
        }
        outFileVars.println();
        
        // Print results
        for (int i = 0; i < finalResults.size(); i++) {
            LoglikScoreResult res = finalResults.get(i);
            outFileVars.print(res.sampleName);
	        for (String strand : strands) {
	        	for (String pwmID : pwms.labels) {
	    	        for (int winNo = 1; winNo <= windowsHash.size(); winNo++) {
	    	        	String featureId = pwmID + "_" + strand + "_" + winNo + "_ROC_LEAF_" + winsWidth;
//	    	        	System.out.println(featureId);
	    	        	outFileVars.print("\t" + res.featureHash.get(featureId));
	    	        	featureId = pwmID + "_" + strand + "_" + winNo + "_ROC_ROOT_" + winsWidth;
	    	        	outFileVars.print("\t" + res.featureHash.get(featureId));
	    	        }
	        	}	        	
	        }
            outFileVars.println();
        }
        outFileVars.flush();
        outFileVars.close();

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
        if (nucsAfterTSS == -1) {
            nucsAfterTSS = seqLength / 2;  // integer arithmetic - automatically floors
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
                    String ref = genomic_parts[genomic_parts.length - 2];
                    int start = Integer.parseInt(genomic_parts[genomic_parts.length - 1]);
                    RefSeqData refData = new RefSeqData(ref, start);
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
                
                if (cmdLine.hasOption("minScores")) {
                    scoreCutoffs_Fname = cmdLine.getOptionValue("minScores");
                }

                if (cmdLine.hasOption("leafOC")) {
                    leafOcRegionsFile = cmdLine.getOptionValue("leafOC");
                }
                
                if (cmdLine.hasOption("rootOC")) {
                    rootOcRegionsFile = cmdLine.getOptionValue("rootOC");
                }

                if (cmdLine.hasOption("nucsAfterTSS")) {
                    nucsAfterTSS = Integer.parseInt(cmdLine.getOptionValue("nucsAfterTSS"));
                }

                if (cmdLine.hasOption("pseudoCounts")) {
                    PseudoCountsVal = Double.parseDouble(cmdLine.getOptionValue("pseudoCounts"));
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
            String cmdLineSyntax = "java -jar GenROCFeaturesTile tfbs_llscan.jar <upstream nucs> <downstream nucs> <window width> <sequence.fasta> <PWM_FILE> <Outfile.rdat>";
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
	        	        
	        // Width of of each window
	        options.addOption(Option.builder()
	        		.longOpt("leafOC")
	        		.hasArg(true)
	        		.argName("String")
	        		.desc("OC regions for leaf")
	        		.required(true)
	        		.type(Double.class)
	        		.build());
	        
	        options.addOption(Option.builder()
	        		.longOpt("rootOC")
	        		.hasArg(true)
	        		.argName("String")
	        		.desc("OC regions for root")
	        		.required(true)
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

	        return options;
	    }


	    private static class OptionComparator <T extends Option> implements Comparator <T> {
	        private ArrayList <String> OPTS_ORDER = new ArrayList <String> ();

	        public OptionComparator () {
	            OPTS_ORDER.add("help");
	            OPTS_ORDER.add("nprocs");
	            OPTS_ORDER.add("BGWIN");
	            OPTS_ORDER.add("minScores");
	            OPTS_ORDER.add("leafOC");
	            OPTS_ORDER.add("rootOC");
	            OPTS_ORDER.add("nucsUp");
	            OPTS_ORDER.add("nucsDown");
	            OPTS_ORDER.add("winWidth");
	            OPTS_ORDER.add("pseudoCounts");
	            OPTS_ORDER.add("nucsAfterTSS");
	            OPTS_ORDER.add("seqName");
	        }

	        public int compare(T o1, T o2) {
	            return OPTS_ORDER.indexOf(o1.getLongOpt()) - OPTS_ORDER.indexOf(o2.getLongOpt());
	        }
	    }

}
