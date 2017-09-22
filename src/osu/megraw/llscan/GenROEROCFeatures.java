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
public class GenROEROCFeatures {
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

    
    public boolean USE_MARKOV_1_BG = true;
	private boolean help = false;
	private String inputSeqName;
	String[] seqLabels;
	private int seqLength = -1;
	HashMap<String, LefRightPair<Integer, Integer>> windowsHash = new HashMap<>();
	
	HashMap <String, RefSeqData> fastaCoords = new HashMap <String, RefSeqData>();
	PWMReturn pwms = null;
	private String rootROCDir;
	private String leafROCDir;
	private HashMap<String, List<Coordinate>> leafOCHash = new HashMap<>();
	private HashMap<String, List<Coordinate>> rootOCHash = new HashMap<>();
	private String roeFWDTable;
	private String roeREVTable;
	private String[] winVarNames;
	private int[] winVarL;
	private int[] winVarR;
	private int[] winVarPWMIndex;
	private String[] winVarStrand;
	private int nWinVars;

    public static double[] bgdefault = {0.25, 0.25, 0.25, 0.25};
    public static double[][] bgM1default = {{0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25},
                                           {0.25, 0.25, 0.25, 0.25}};

	public static void main(String[] args) {
		try {
			GenROEROCFeatures genROCFeatures = new GenROEROCFeatures();
			genROCFeatures.parseArgs(args);
			genROCFeatures.readPromoterSeq();
			genROCFeatures.readPWMsInfo();
//			genROCFeatures.readOCRegions();
			genROCFeatures.createWins();
			genROCFeatures.getOCScores();
			System.out.println("DONE!");
		} catch (Exception e) {
			e.printStackTrace();
		}

	}
	
	private void readOCRegions(String tssName) {
		// Read root and lead oc regions separately
		try {
			String leafInfile = leafROCDir + "/" + tssName + ".overlap";
			String rootInfile = rootROCDir + "/" + tssName + ".overlap";
			
			BufferedReader leafFile = new BufferedReader(new FileReader(leafInfile));
			BufferedReader rootFile = new BufferedReader(new FileReader(rootInfile));
			leafOCHash = new HashMap<String, List<Coordinate>>();
			String line = null;
			while ((line = leafFile.readLine()) != null) {
				line.trim();
				String[] parts = line.split("\t");
				String winFacInfo = parts[0].split("%")[1];
				String chrom = parts[2];
				int relLeft = Integer.valueOf(parts[5]);
				int relRight = Integer.valueOf(parts[6]);
				
				Coordinate coordinate = new Coordinate(chrom, relLeft, relRight);
				
				List<Coordinate> coordinateList = leafOCHash.get(winFacInfo);
				if (coordinateList == null) {
					coordinateList = new ArrayList<>();
				}
				coordinateList.add(coordinate);
				leafOCHash.put(winFacInfo, coordinateList);
			}
			
			leafFile.close();
			
			line = null;
			rootOCHash = new HashMap<String, List<Coordinate>>();
			while ((line = rootFile.readLine()) != null) {
				line.trim();
				String[] parts = line.split("\t");
				String winFacInfo = parts[0].split("%")[1];
				String chrom = parts[2];
				int relLeft = Integer.valueOf(parts[5]);
				int relRight = Integer.valueOf(parts[6]);
				
				Coordinate coordinate = new Coordinate(chrom, relLeft, relRight);
				
				List<Coordinate> coordinateList = rootOCHash.get(winFacInfo);
				if (coordinateList == null) {
					coordinateList = new ArrayList<>();
				}
				coordinateList.add(coordinate);
				rootOCHash.put(winFacInfo, coordinateList);
			}
			
			rootFile.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
		
	}

	private void createWins() throws IOException {
        DoubleMatReturn facInfo_FWD = Load.loadDoubleMatrixFile(new File(roeFWDTable), true);
        DoubleMatReturn facInfo_REV = Load.loadDoubleMatrixFile(new File(roeREVTable), true);

        // Process factorInfo file into variable window definitions
        int NPI = 5; // Number of windows used to cover the peak region
        int N = NPI + 2; // total number of windows

        DoubleMatReturn[] facInfoArr = {facInfo_FWD, facInfo_REV};
        String[] strandArr = {"FWD", "REV"};

        nWinVars = 0;
        for (int fi = 0; fi < facInfoArr.length; fi++) {
            DoubleMatReturn facInfo = facInfoArr[fi];
            nWinVars += N*(facInfo.rowLabels.length - facInfo.nflags); // N * number of factors with numeric parameters
        }
        winVarNames = new String[nWinVars];
        winVarL = new int[nWinVars];
        winVarR = new int[nWinVars];
        winVarPWMIndex = new int[nWinVars];
        winVarStrand = new String[nWinVars];
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

	}

	private void getOCScores() throws Exception {
		List<LoglikScoreResult> final_results = new ArrayList<>();
        // Iterate over sequences 
        for (int i = 0; i < sequenceCharArr.length; i++) {
        	String seqName = seqLabels[i];
        	
        	// pull related open regions by seq name from input directory
        	readOCRegions(seqName+ "_0"); 
        	//
        	ComputeLoglikScoreROEROC computeLoglikScoreROEROC = new ComputeLoglikScoreROEROC(sequenceCharArr[i], scoreCutOffs, nucsAfterTSS, BG_WIN, seqName, pwms, winVarNames, 
        			winVarL,  winVarR,  winVarPWMIndex,  winVarStrand, leafOCHash, rootOCHash, nWinVars);
        	LoglikScoreResult loglikScoreResult = computeLoglikScoreROEROC.call();
        	final_results.add(loglikScoreResult);
        }
        process_results(final_results);
	}
	
	private void process_results(List<LoglikScoreResult> finalResults) throws Exception {
    	// Write out features
        PrintWriter outFileVars = new PrintWriter(new FileWriter(out_Fname));
        
        // Print header
		for (int wv = 0; wv < nWinVars; wv++) {
        	String featureId = winVarNames[wv] + "_LEAF";
        	outFileVars.print("\t" + featureId);
        	featureId = winVarNames[wv] + "_ROOT";
        	outFileVars.print("\t" + featureId);
    	        	
        }
        outFileVars.println();
        
        // Print results
        for (int i = 0; i < finalResults.size(); i++) {
            LoglikScoreResult res = finalResults.get(i);
            outFileVars.print(res.sampleName);
    		for (int wv = 0; wv < nWinVars; wv++) {
            	String featureId = winVarNames[wv] + "_LEAF";
            	outFileVars.print("\t" + res.featureHash.get(featureId));
            	featureId = winVarNames[wv] + "_ROOT";
            	outFileVars.print("\t" + res.featureHash.get(featureId));
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
                	roeFWDTable = args[0];
                    roeREVTable = args[1];
                    inputSeqName = args[2];
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

                if (cmdLine.hasOption("leafOC")) {
                    leafROCDir = cmdLine.getOptionValue("leafOC");
                }
                
                if (cmdLine.hasOption("rootOC")) {
                    rootROCDir = cmdLine.getOptionValue("rootOC");
                }

                if (cmdLine.hasOption("nucsAfterTSS")) {
                    nucsAfterTSS = Integer.parseInt(cmdLine.getOptionValue("nucsAfterTSS"));
                }

                if (cmdLine.hasOption("pseudoCounts")) {
                    PseudoCountsVal = Double.parseDouble(cmdLine.getOptionValue("pseudoCounts"));
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
