package osu.megraw.llscan;

import java.util.Comparator;

public class ScanResultComparator implements Comparator <ScanResult> {
    boolean byPWM = true;

    public ScanResultComparator () {
    }

    public ScanResultComparator (String sortType) {
        if (sortType.equals("PWM")) {
            byPWM = true;
        } else if (sortType.equals("SEQ")) {
            byPWM = false;
        }
    }

    public int compare(ScanResult o1, ScanResult o2) {
        int order = 0;

        if (o1.strand.equals(o2.strand)) {
            // sort results by the PWM used in scanning
            if (byPWM) {
                order = o1.pwmLabel.compareTo(o2.pwmLabel);
            // sort results by the sequence scanned
            } else {
                order = o1.seqLabel.compareTo(o2.seqLabel);
            }
        } else {
            // Ensure FWD strand is always listed first
            if (o1.strand.equals("FWD")) {
                order = -1;
            } else {
                order = 1;
            }
        }

        return order;
    }
}
