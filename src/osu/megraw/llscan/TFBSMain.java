package osu.megraw.llscan;

import java.util.HashSet;

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
public class TFBSMain {
    public static HashSet <String> commands;

    public static void main(String[] args) throws java.io.IOException, java.lang.InterruptedException, osu.megraw.llscan.BadCharException {
        String[] command_list = { "Background", 
                                  "GenFeatures", 
                                  "GenFeaturesTiledWins",
                                  "MakePWMs", 
                                  "PWMprint", 
                                  "ROEFinder", 
                                  "Scan", 
                                  "Thresholds",
                                  "GenFeaturesByNT"};

        HashSet <String> commands = new HashSet <String>();
        for (int i = 0; i < command_list.length; i++) commands.add(command_list[i]);

        String usage = "usage: java -jar tfbs_scan.jar <command>\n\n" +
                       "This is the main interface to the TFBS-scanning utilties package.\n" +
                       "This package inlcudes the following commands:\n\n" +
                       "Scan:        The main scanning utility which scans TFBS sites.\n\n" +
                       "ROEFinder:   Generates regions of enrichment for a list of TFBSs.\n\n" +
                       "GenFeatures: Generates TFBS features within regions of enrichment and\n" +
                       "             additional sequences features (GC/GA/CA content).\n\n" +
                       "GenFeaturesTiledWins: Generates TFBS features within tiled windows \n" +
                       "             additional sequences features (GC/GA/CA content).\n\n" +
                       "Additional utility commands exist to aid researchers in general scaning:\n\n" +
                       "Background:  This utiltity generates 0th and 1st order background Markov\n" +
                       "             Models for a FASTA formatted sequence file.\n\n" +
                       "Thresholds:  This utility is an interface to the get-WMM-cutoff binary\n" +
                       "             which generates loglikelihood score cutoffs at a specified\n" +
                       "             critical value (FP, FN, and MAX - the maximum of FP & FN).\n\n" +
                       "PWMprint:    This utility generates normalized PWMs and can add 'pseudo counts'\n" +
                       "             to matrices prior to normalization to avoid problems with logs of '0'.\n\n" +
                       "MakePWMs:    This utility generates PWMs from a list of sequences.\n\n" +
                       "GenFeaturesByNT: Generates loglik scores for each individual nucleotide within ROE windows\n\n";

        if (args.length < 1) {
            System.out.print(usage);
            return;
        }

        String command = args[0];

        if (commands.contains(command)) {
            // Strip off first command line argument, and pass on the rest of the arguments to the main java program
            String[] new_args = new String[args.length - 1];
            for (int i = 1; i < args.length; i++) {
                new_args[i - 1] = args[i];
            }
            // Primary scanning programs
            if (command.equals("Scan")) {
                Scan.main(new_args);
            } else if (command.equals("ROEFinder")) {
                ROEFinder.main(new_args);
            } else if (command.equals("GenFeatures")) {
                GenFeatures.main(new_args);
            } else if (command.equals("GenFeaturesTiledWins")) {
                GenFeaturesTiledWindows.main(new_args);
            // Additional helper utilities
            } else if (command.equals("Background")) {
                Background.main(new_args);
            } else if (command.equals("Thresholds")) {
                Thresholds.main(new_args);
            } else if (command.equals("PWMprint")) {
                PWMprint.main(new_args);
            } else if (command.equals("MakePWMs")) {
                MakePWMs.main(new_args);
            } else if (command.equals("GenFeaturesByNT")) {
                GenFeaturesByNT.main(new_args);
            }
        } else {
            System.out.print(usage);
            return;
        }
    }
}
