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
public class MakePWMs {
    public MakePWMs() {
    }

    public static void main(String[] args) throws java.io.IOException {
        MakePWMs makePWMs = new MakePWMs();

        String usage = "usage: java -jar tfbs_scan.jar MakePWMs <SEQUENCE FILE> <OUTFILE>\n\n" +
                       "The MakePWMs utility generates a position weight matrix file and outputs\n" +
                       "a frequency matrix to <OUTFILE> for each matrix within <SEQUENCE FILE>.\n\n" +
                       "The format for the <SEQUENCE FILE> is as follows for each matrix (minus the lines with '----'):\n" +
                       "------------------------------\n" +
                       ">pwm_id\n" +
                       "SEQUENCE1\n" +
                       "SEQUENCE2\n" +
                       "SEQUENCE3\n" +
                       "...\n" +
                       "...\n" +
                       "\n" +
                       "------------------------------\n\n" +			
                       "Take note of the empty line at the end.  Each list of sequences for each pwm *must*\n" +
                       "end with an empty line *including* the final matrix sequence list! Also, if a matrix\n" +
                       "has sequences of variable length, ensure that the longest sequence is listed first!\n" +
                       "The first sequences sets the length of the PWM generated.\n";

        if (args.length < 2) {
            System.out.println(usage);
            return;
        }

        // Define input files
        String inFile = args[0];

        // Define output files
        String outFile = args[1];

        // Read sequences
        FreqReturn fr = Load.loadPWMSeqFile(inFile);

        // Print out frequency file
        Print.printFreqsToFile(fr, outFile);
    }
}
