package osu.megraw.llscan;

import java.util.Vector;

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
public class PWMHitRec {
    public double proportionHit;
    public int[][] hitLocs;
    public Vector locVec;
    public double[][] hitScores;
    public Vector scoreVec;

    public PWMHitRec() {
    }

    public PWMHitRec(double in_proportionHit, int[][] in_hitLocs, Vector in_locVec, double[][] in_hitScores, Vector in_scoreVec) {
        proportionHit = in_proportionHit;
        hitLocs = in_hitLocs;
        locVec = in_locVec;
        hitScores = in_hitScores;
        scoreVec = in_scoreVec;
    }
}
