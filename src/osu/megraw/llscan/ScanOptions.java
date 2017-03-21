package osu.megraw.llscan;

import java.util.ArrayList;
import java.util.Comparator;

// Use Apache Commons CLI Parser
import org.apache.commons.cli.Option;
import org.apache.commons.cli.Options;

public class ScanOptions {
    public ScanOptions() {
    }

    public static Options buildOptions() {
        Options options = new Options();

        // help
        options.addOption(Option.builder("h")
            .longOpt("help")
            .hasArg(false)
            .desc("Print this helpful help message")
            .required(false)
            .build());

        // nprocs
        options.addOption(Option.builder("n")
            .longOpt("nprocs")
            .hasArg(true)
            .desc("Number of threads to use for loglikelihood scan. default: 1")
            .argName("INT")
            .required(false)
            .build());

        // bg_Filename
        options.addOption(Option.builder()
            .longOpt("bg")
            .hasArg()
            .argName("FILE")
            .desc("Mononucleotide background model file. default: generate from input sequences")
            .required(false)
            .build());

        // usebg2
        options.addOption(Option.builder()
            .longOpt("usebg2")
            .hasArg(false)
            .desc("Use dinucleotide background model. If no model file supplied, automatically generates model from input sequences.")
            .required(false)
            .build());

        // bg2_Filename
        options.addOption(Option.builder()
            .longOpt("bg2")
            .hasArg()
            .argName("FILE")
            .desc("Dinucleotide backgroud model file. deafult: generate from input sequences (automatically sets --usebg2 flag)")
            .required(false)
            .build());

        // strand
        options.addOption(Option.builder()
            .longOpt("strand")
            .hasArg()
            .argName("FWD|REV|BOTH")
            .desc("Strand of sequence to scan. default: BOTH")
            .required(false)
            .build());

        // minScores_Fname
        options.addOption(Option.builder()
            .longOpt("minScores")
            .hasArg()
            .argName("FILE")
            .desc("Score threshold file. default: 0.0 for all PWMS")
            .required(false)
            .build());

        // PseudoCountsVal
        options.addOption(Option.builder()
            .longOpt("pseudoCounts")
            .hasArg()
            .argName("NUM")
            .desc("Total pseudo-counts to add to frequency matrix prior to normalization. default: 0.25")
            .type(Double.class)
            .required(false)
            .build());

        // nucsAfterTSS
        options.addOption(Option.builder()
            .longOpt("nucsAfterTSS")
            .hasArg()
            .argName("INT")
            .desc("Number of nucleotides from the right of sequences defined as the TSS. default: 0")
            .type(Integer.class)
            .required(false)
            .build());

        return options;
    }

    public class OptionComparator <T extends Option> implements Comparator<T> {
        private ArrayList <String> OPTS_ORDER = new ArrayList <String> ();

        public OptionComparator () {
            OPTS_ORDER.add("help");
            OPTS_ORDER.add("nprocs");
            OPTS_ORDER.add("bg");
            OPTS_ORDER.add("usebg2");
            OPTS_ORDER.add("bg2");
            OPTS_ORDER.add("strand");
            OPTS_ORDER.add("minScores");
            OPTS_ORDER.add("pseudoCounts");
            OPTS_ORDER.add("nucsAfterTSS");
        }

        public int compare(T o1, T o2) {
            return OPTS_ORDER.indexOf(o1.getLongOpt()) - OPTS_ORDER.indexOf(o2.getLongOpt());
        }
    }
}
