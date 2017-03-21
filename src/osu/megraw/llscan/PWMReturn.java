package osu.megraw.llscan;

import java.util.HashMap;

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
public class PWMReturn {

    public HashMap<String,Integer> labelIndex = new HashMap<String,Integer>();

    public String[] labels;
    public String[] headers;
    public double[][][] pwms;
    public double[][][] revcmp_pwms;

    public PWMReturn() {
    }

    public PWMReturn(String[] in_labels, String[] in_headers, double[][][] in_pwms) {
        labels = in_labels;
        headers = in_headers;
        pwms = in_pwms;

        for (int i = 0; i < in_labels.length; i++) {
            labelIndex.put(in_labels[i], i);
        }
    }

    public void ComputeRevCmpPWMs() {
        revcmp_pwms = new double[pwms.length][][];
        for (int nmat = 0; nmat < pwms.length; nmat++) {
            revcmp_pwms[nmat] = Utils.PWMRevCmp(pwms[nmat]);
        }
    }    
}
