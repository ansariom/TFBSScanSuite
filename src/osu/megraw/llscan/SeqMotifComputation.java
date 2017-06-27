package osu.megraw.llscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileWriter;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

public class SeqMotifComputation {
	
	private static String sitesTableFile;
	private static String pwmFile;
	private static String histDir;
	private static String outFile;
	private static double[] bgFrequenciesMarkov0 = {0.3486, 0.1599,	0.1498,	0.3417 };
	
	public static void main(String[] args) {
		if (args.length < 4) {
			System.out.println("SeqMotifComputation <sitesTableFile> <PWM_FILE> <HIST_DIR> <OUTFILE>");
			System.exit(1);
		}
		sitesTableFile = args[0];
		pwmFile = args[1];
		histDir = args[2];
		outFile = args[3];
		
		SeqMotifComputation seqMotifComputation = new SeqMotifComputation();
		seqMotifComputation.processSiteTable(sitesTableFile);	
		
	}
	
	private void processSiteTable(String siteTable) {
		HashMap<String, double[][]> pwmMatrixHash = new HashMap<String, double[][]>();
		
		try {
			FileWriter writer = new FileWriter(new File(outFile));
			InputStream inputStream = new FileInputStream(new File(siteTable));
			BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(inputStream));
			
			String header = bufferedReader.readLine();
			writer.write(header + "\tmotif_seq_revcom\tflanking_seq_revcom\tFP_rate_M0\tFN_rate_M0\tloglik_score_M0\n");

			String line = null;
			while ((line = bufferedReader.readLine()) != null) {
				String[] parts = line.split("\t");
				String motifSeq = parts[25];
				String flankingSeq = parts[26];
				String winStrand = parts[14];
				String pwmName = parts[13];
				
				double[][] pwm = pwmMatrixHash.get(pwmName);
				if (pwm == null) {
					pwm = Utils.getPWM(pwmFile, pwmName);
					pwmMatrixHash.put(pwmName, pwm);
				}
				
				if (winStrand.equalsIgnoreCase("REV"))
					motifSeq = Utils.getRevComplement(motifSeq);
				double logLikeScoreM0 = Utils.computeLogLikScoreForStringUsingM0(motifSeq, motifSeq.length(), pwm, bgFrequenciesMarkov0);
				double fpRateM0 = computeFP(logLikeScoreM0, loadBGHistogram(histDir + "/" + pwmName + ".hbin"));
				double fnRateM0 = computeFP(logLikeScoreM0, loadFGHistogram(histDir + "/" + pwmName + ".hbin"));
				writer.write(line.substring(0, line.length() - 1) + "\t" +  
						Utils.getRevComplement(motifSeq) + "\t" + Utils.getRevComplement(flankingSeq) + "\t" +
						fpRateM0 + "\t" + fnRateM0 + "\t" + logLikeScoreM0 + "\n");
				
				
			}
			writer.flush();
			writer.close();
			inputStream.close();
			bufferedReader.close();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public static double computeFP(double loglikScore, List<String> bins) {
		double fp = 0.0;
		BigDecimal auP = new BigDecimal(0);
		
		for (int i = 0; i < bins.size(); i++) {
			BigDecimal lb = new BigDecimal(bins.get(i).split("_")[0]);
			if (new BigDecimal(loglikScore).compareTo(lb) == 1) {
				BigDecimal binHeight = new BigDecimal(bins.get(i).split("_")[1]);
				auP = auP.add(binHeight);
			} else
				break;
		}
		
		fp = new BigDecimal(1).subtract(auP).doubleValue();
		return fp;
	}

	/**
	 * Background histograms start from first line up to (nbin) (depends on number of bins)
	 * @param histFile
	 * @return
	 */
	public static List<String> loadBGHistogram(String histFile) {
		List<String> histList = new ArrayList<String>();
		
		try {
			BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(histFile))));
			String line = bufferedReader.readLine();
			
			while (line != null && !line.equalsIgnoreCase("")) {
				BigDecimal lb = new BigDecimal(line.split("\t")[0].split(",")[1]);
				BigDecimal ub = new BigDecimal(line.split("\t")[0].split(",")[2]);
				histList.add(lb.toString() + "_" + line.split("\t")[1]);
				line = bufferedReader.readLine();
			}
			bufferedReader.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return histList;
	}

	/**
	 * Foreground histograms start from line (nbin)+ (depends on number of bins)
	 * @param histFile
	 * @return
	 */
	public static List<String> loadFGHistogram(String histFile) {
		List<String> histList = new ArrayList<String>();
		
		try {
			BufferedReader bufferedReader = new BufferedReader(new InputStreamReader(new FileInputStream(new File(histFile))));
			String line = bufferedReader.readLine();
			
			while (line != null && !line.equalsIgnoreCase("")) {
				BigDecimal lb = new BigDecimal(line.split("\t")[0].split(",")[1]);
				BigDecimal ub = new BigDecimal(line.split("\t")[0].split(",")[2]);
				histList.add(lb.toString() + "_" + line.split("\t")[1]);
				line = bufferedReader.readLine();
			}
			bufferedReader.close();
		} catch (Exception e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		return histList;
	}

}
