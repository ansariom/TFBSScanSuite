package osu.megraw.llscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.List;
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
public class Thresholds {
    public Thresholds() {
    }

    public static void main(String[] args) throws IOException , InterruptedException {
        // Required parameters
        String matrix_pwm = "";
        String background_fasta = "";
        double critical_value = 0.0;

        // Optional parameters
        boolean help = false;                 // -h
        int nprocs = 1;                       // -n
        String decision_type = "FP";          // -d 
        boolean use_prob = false;             // -P
        double pwm_pseudo_counts = 0.25;      // -C
        double pseudo_counts = 0.0;           // -c
        int histo_bins = 500;                 // -b
        int histo_samples = 50000;            // -s
        int order = 1;                        // -o
        boolean use_midpoint = false;         // -M
        boolean output_histograms = false;    // -H
        String histogram_output_file; 
        String background_output_file = "";   // -B
        String threshold_output_file = "";    // -T

        // Setup CLI to simplify Threshold interface;
        Options options = buildOptions();
        CommandLineParser parser = new DefaultParser();
        CommandLine cmdLine;

        try {
            cmdLine = parser.parse(options, args);

            // Remove arguments that were parsed
            args = new String[cmdLine.getArgList().size()];
            args = cmdLine.getArgList().toArray(args);

            if (cmdLine.hasOption("h")) {
                help = true;
            } else {
                // Handle required arguments
                if (args.length < 3) {
                    throw new ParseException("Must supply a matrix file, background fasta file, and a critical value.");
                } else {
                    matrix_pwm = args[0];
                    background_fasta = args[1];
                    critical_value = Double.parseDouble(args[2]);
                }

                // -n, -d, -P, -C, -c, -b, -s, -o, -M, -H, -B, -T
                if (cmdLine.hasOption("n")) {
                    nprocs = Integer.parseInt(cmdLine.getOptionValue("n"));
                    if (nprocs > Runtime.getRuntime().availableProcessors()) {
                        System.err.println("Warning: your system does not have " + nprocs + " CPUs: setting nprocs to " + Runtime.getRuntime().availableProcessors() + ".");
                        nprocs = Runtime.getRuntime().availableProcessors();
                    }
                    if (nprocs < 1) {
                        System.err.println("Warning: nprocs must be a value >= 1. Setting nprocs = 1");
                        nprocs = 1;
                    }
                }

                if (cmdLine.hasOption("d")) {
                    decision_type = cmdLine.getOptionValue("d");
                    if (!decision_type.equals("FP") && !decision_type.equals("FN") && !decision_type.equals("MAX")) {
                        throw new ParseException("-d must have value FP, FN, or MAX.");
                    }
                }
                if (cmdLine.hasOption("P")) {
                    use_prob = true;
                }
                if (cmdLine.hasOption("C")) {
                    pwm_pseudo_counts = Double.parseDouble(cmdLine.getOptionValue("C"));
                }
                if (cmdLine.hasOption("c")) {
                    pseudo_counts = Double.parseDouble(cmdLine.getOptionValue("c"));
                }
                if (cmdLine.hasOption("b")) {
                    histo_bins = Integer.parseInt(cmdLine.getOptionValue("b"));
                }
                if (cmdLine.hasOption("s")) {
                    histo_samples = Integer.parseInt(cmdLine.getOptionValue("s"));
                }
                if (cmdLine.hasOption("o")) {
                    order = Integer.parseInt(cmdLine.getOptionValue("o"));
                }
                if (cmdLine.hasOption("M")) {
                    use_midpoint = true;
                }
                if (cmdLine.hasOption("H")) {
                    output_histograms = true;
                }
                if (cmdLine.hasOption("B")) {
                    background_output_file = cmdLine.getOptionValue("B");
                }
                if (cmdLine.hasOption("T")) {
                    threshold_output_file = cmdLine.getOptionValue("T");
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
            Comparator <Option> opt_order = new Thresholds.OptionComparator <Option> ();
            String cmdLineSyntax = "java -jar tfbs_scan.jar Thresholds <matrix.pwm> <background.fasta> <critical-value> where <critical-value> is the FP rate";
            String header = "";
            String footer = "";
            HelpFormatter formatter = new HelpFormatter();
            formatter.setOptionComparator(opt_order);
            formatter.setWidth(80);
            formatter.printHelp(cmdLineSyntax, header, options, footer, true);
            System.exit(0);
        }

        // The strings needed to build the get-WMM-cutoff command to run
        String wmmBin = "get-WMM-cutoff"; // executable - assumes 'get-WMM-cutoff' is in your PATH
        ArrayList<String> wmmOptions = buildWMMOptions(decision_type, use_prob, pseudo_counts, histo_bins, histo_samples, order, use_midpoint, background_output_file); // Options user has supplied

        // Read PWM file
        PWMReturn pwms;
        //pwms = Load.loadPWMFile(matrix_pwm, pwm_pseudo_counts);
        pwms = Load.loadPWMFileSimpleHeader(matrix_pwm, pwm_pseudo_counts);

        // Setup thread management for parallel processing
        ThreadPoolExecutor threadPool = new ThreadPoolExecutor(nprocs, nprocs, 4, TimeUnit.SECONDS, new LinkedBlockingQueue());
        List<Future<String>> results = new ArrayList<Future<String>>();

        // Go through each PWM and calculate a threshold
        for (int i = 0; i < pwms.pwms.length; i++) {
            // Build the command-line argument needed to run 'get-WMM-cutoff'
            List <String> wmmArgs = new ArrayList <String> ();
            wmmArgs.add(wmmBin);
            wmmArgs.addAll(wmmOptions);

            // Specify the histogram bin for this file
            if (output_histograms) {
                wmmArgs.add("-h");
                wmmArgs.add(pwms.labels[i] + ".hbin");
            }

            wmmArgs.add(pwms.labels[i] + ".wmm"); // Name of the wmm file - made and destroyed in ThresholdRunner
            wmmArgs.add(background_fasta);
            wmmArgs.add(Double.toString(critical_value));

            ThresholdRunner run = new ThresholdRunner(wmmArgs, pwms.labels[i], pwms.pwms[i]);

            try {
                results.add(threadPool.submit(run));
            } catch (Exception e) {
                System.err.println("Error processing pwm " +
                                    pwms.labels[i] +
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

        // Print results!
        PrintWriter writer = new PrintWriter(System.out);
        if (!threshold_output_file.equals("")) {
            writer.close();
            writer = new PrintWriter(threshold_output_file);
        }
        for (int i = 0; i < results.size(); i++) {
            Future <String> result = results.get(i);
            String res = null;
            try {
                res = result.get(); // blocks
            } catch (Exception e) {
                System.err.println("Error processing pwm: " + e.getMessage() + ". Quitting.");
                break;
            }
            writer.println(res);
        }
        writer.close();

        // Clean up our thread mess
        threadPool.shutdownNow();
        while(!threadPool.isTerminated()) {
            threadPool.awaitTermination(2, TimeUnit.SECONDS);
        }
    }

    private static class ThresholdRunner implements Callable <String> {
        List <String> command;
        String pwm_id;
        double [][] pwm;

        public ThresholdRunner (List <String> command, String pwm_id, double[][] pwm) {
            this.command = command;
            this.pwm_id = pwm_id;
            this.pwm = pwm;
        }

        public String call() throws IOException, InterruptedException {
            // File to hold the PWM transformed in the WMM format needed for Bill Majoro's get-WMM-cutoff tool
            String matrix_wmm = pwm_id + ".wmm";
    
            // Convert our PWM frequency matrix file to a WMM file on the fly
            writePWMtoWMMFile(matrix_wmm, pwm);
    
            // Kick off process
            ProcessBuilder wmmProcessObj = new ProcessBuilder(command);
            Process wmmProcess = wmmProcessObj.start();
            wmmProcess.waitFor();
    
            // Read results from process
            BufferedReader wmmResult = new BufferedReader(new InputStreamReader(wmmProcess.getInputStream()));
            String result = wmmResult.readLine();
            wmmResult.close();

            // Clean up intermediate file data
            File wmmFile = new File(matrix_wmm);
            wmmFile.delete();

            String result_line = pwm_id + "\t" + result;
    
            return result_line;
        }
    }

    private static ArrayList<String> buildWMMOptions(String decision_type, boolean use_prob, double pseudo_counts, int histo_bins, int histo_samples, int order, boolean use_midpoint, String background_output_file){
        ArrayList <String> wmmOptions = new ArrayList<String> ();

        wmmOptions.add("-d");
        wmmOptions.add(decision_type);

        if (use_prob) {
            wmmOptions.add("-P");
        }

        wmmOptions.add("-c");
        wmmOptions.add(Double.toString(pseudo_counts));

        wmmOptions.add("-b");
        wmmOptions.add(Integer.toString(histo_bins));

        wmmOptions.add("-s");
        wmmOptions.add(Integer.toString(histo_samples));

        wmmOptions.add("-o");
        wmmOptions.add(Integer.toString(order));

        if (!use_midpoint) {
            wmmOptions.add("-M");
        }

        if (!background_output_file.equals("")) {
            wmmOptions.add("-B");
            wmmOptions.add(background_output_file);
        }

        return wmmOptions;
    }

    private static void writePWMtoWMMFile(String outfile, double [][]pwm) throws IOException {
        PrintWriter writer = new PrintWriter(outfile);

        // Write out the header of the WMM file
        writer.println("WMM");
        writer.println("PROMOTER");
        // 0 = cutoff (doesn't appear to be used in get-WMM-cutoff)
        // pwm.length = the length of the pwm sequence
        // 5 = number of columns represented as A,C,G,N, and T (in that order) - N isn't found in the PWM matrix format, only four columns exist
        writer.println("0 " + pwm.length + " 5");

        // pwm[0].length = consensus window length which is equal to the pwm sequence length
        // 0 and 0 refer to the consensus offset (which we aren't setting, hence needs to be 0),
        // and the consensus sequence length, respectively - neither are needed so just set to 0
        // + refers to the strand, and all threshold calculations are made with respect to the FWD strand
        writer.println(pwm.length + " 0 0 +");

        for (int j = 0; j < pwm.length; j++) {
            String logPs = "";
            for (int i = 0; i < 5; i++) {
                // Write out 'N' value as '-inf' since we don't need it
                if (i == 3) {
                    logPs += "\t-inf";
                } else {
                    int nt_index = i;

                    // Adjust nt_index to account for the '-inf' value needed for WMM matrix which isn't present in PWM matrix
                    if (nt_index == 4) {
                        nt_index = 3;
                    }
                    // Only prepend a tab to values betyond the first value
                    if (nt_index != 0) {
                        logPs += "\t";
                    }
                    logPs += Double.toString(Math.log(pwm[j][nt_index]));
                }
            }
            writer.println(logPs);
        }
        writer.close();
    }

    private static class OptionComparator <T extends Option> implements Comparator<T> {
        private ArrayList <String> OPTS_ORDER = new ArrayList <String> ();

        public OptionComparator () {
            OPTS_ORDER.add("help");
            OPTS_ORDER.add("n");
            OPTS_ORDER.add("d");
            OPTS_ORDER.add("P");
            OPTS_ORDER.add("C");
            OPTS_ORDER.add("c");
            OPTS_ORDER.add("b");
            OPTS_ORDER.add("s");
            OPTS_ORDER.add("o");
            OPTS_ORDER.add("M");
            OPTS_ORDER.add("T");
            OPTS_ORDER.add("H");
            OPTS_ORDER.add("B");
        }

        public int compare(T o1, T o2) {
            return OPTS_ORDER.indexOf(o1.getOpt()) - OPTS_ORDER.indexOf(o2.getOpt());
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
            .hasArg(true)
            .argName("INT")
            .desc("Number of processors to use (default: 1)")
            .required(false)
            .type(Integer.class)
            .build());
            
        // type of critical value
        options.addOption(Option.builder("d")
            .hasArg(true)
            .desc("interpret <critical-value> as: FP = false positive rate, FN = false negative rate, MAX = max of FP and FN. default: FP")
            .argName("FP|FN|MAX")
            .required(false)
            .build());

        // type of probabilities to use
        options.addOption(Option.builder("P")
            .hasArg(false)
            .desc("use raw probabilities instead of LLR ratios")
            .required(false)
            .build());

        // add probability pseudocounts to your PWM
        options.addOption(Option.builder("C")
            .hasArg(true)
            .argName("NUM")
            .desc("add <NUM> pseudcounts to your PWM frequency scores prior to normalization (default=0.25)")
            .type(Double.class)
            .required(false)
            .build());

        // pseudocounts to add
        options.addOption(Option.builder("c")
            .hasArg(true)
            .argName("NUM")
            .desc("use pseudocount <NUM> for histogram bins (default=0.0)")
            .type(Double.class)
            .required(false)
            .build());

        // number of bins for histograms
        options.addOption(Option.builder("b")
            .hasArg(true)
            .argName("INT")
            .desc("use <INT> bins in the histograms (default=500)")
            .type(Integer.class)
            .required(false)
            .build());

        // number of samples for probability cutoff computations
        options.addOption(Option.builder("s")
            .hasArg(true)
            .argName("INT")
            .desc("use <INT> samples when computing histogram (default=50000)")
            .type(Integer.class)
            .required(false)
            .build());

        // background model order
        options.addOption(Option.builder("o")
            .hasArg(true)
            .argName("INT")
            .desc("order of background model (default=1)")
            .type(Integer.class)
            .required(false)
            .build());

        // midpoint rule
        options.addOption(Option.builder("M")
            .hasArg(false)
            .desc("Apply midpoint rule, even for separable distributions")
            .required(false)
            .build());

        // threshold output file
        options.addOption(Option.builder("T")
            .hasArg(true)
            .argName("outfile")
            .desc("Output threshold values to <outfile> (default: print to STDOUT)")
            .required(false)
            .build());

        // histogram output file
        options.addOption(Option.builder("H")
            .hasArg(false)
            .desc("dump histogram files (outputs histogram to <matrix_id.hbin>)")
            .required(false)
            .build());

        // background output file
        options.addOption(Option.builder("B")
            .hasArg(true)
            .argName("outfile")
            .desc("save background model to <outfile>")
            .required(false)
            .build());

        return options;
    }
}
