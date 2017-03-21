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
public class DoubleColReturn {

    public String[] labels;
    public double[] values;

    public DoubleColReturn() {
    }

    public DoubleColReturn(String[] in_labels, double[] in_values) {
        labels = in_labels;
        values = in_values;
    }

    // Helper constructor to create default array of zeroes
    public DoubleColReturn(String[] in_labels) {
        this.labels = in_labels;
        this.values = new double[this.labels.length];
        for (int i = 0; i < this.values.length; i++) {
            this.values[i] = 0.0;
        }
    }
}
