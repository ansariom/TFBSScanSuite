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
public class SeqArrayReturn {

    public String[] labels;
    public String[] headers;
    public char[][][] seqs;

    public SeqArrayReturn() {
    }

    public SeqArrayReturn(String[] in_labels, String[] in_headers, char[][][] in_seqs) {
        labels = in_labels;
        headers = in_headers;
        seqs = in_seqs;
    }
}
