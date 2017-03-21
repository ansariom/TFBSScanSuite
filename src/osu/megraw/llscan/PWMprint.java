package osu.megraw.llscan;

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
public class PWMprint {
    public PWMprint() {
    }

    public static void main(String[] args) throws java.io.IOException {
        PWMprint pwmp = new PWMprint();

        PWMReturn pwms;

        String usage = "usage: java -jar tfbs_scan.jar PWMprint <pseudo-counts> <PWM FILE> <PWM OUTFILE>\n\n" +
                       "The PWMprint utility generates a normalized position weight matrix file by reading in\n" +
                       "the <PWM FILE>, adding <pseudo-counts> (which is a <DOUBLE>, i.e a number allowing decimal\n" +
                       "places, e.g. 3.14) to each value in the matrix, and generating a final normalized matrix\n" +
                       "in <PWM OUTFILE> outputing probabilities for each nucleotide at each position. This utility\n" +
                       "current only accepts the general PWM Matrix format supported by the TFBS scanning suite.\n\n" +
                       "The format for each matrix in the <PWM FILE> is as follows (minus the lines with '----'):\n" +
                       "------------------------------\n" +
                       "> pwm_id\n" +
                       "=\\t<DOUBLE>\\t<DOUBLE>\\t\\t<DOUBLE>\\t<DOUBLE>\n" +
                       "...\n" +
                       "...\n" +
                       "\n" +
                       "------------------------------\n\n" +
                       "Take note of the final line - each matrix *must* be followed by an empty line, *including*\n" +
                       "the last matrix in the <PWM FILE>!\n";

        if (args.length < 3) {
            System.out.println(usage);
            return;
        }

        double PseudoCountsVal = Double.parseDouble(args[0]);
        String pwms_Fname = args[1];
        String out_Fname = args[2];

        // Read PWM file
        pwms = Load.loadPWMFileSimpleHeader(pwms_Fname, PseudoCountsVal);

        // Print PWM file
        Print.printPWMsToFile(pwms, out_Fname, Print.df4);
    }
}
