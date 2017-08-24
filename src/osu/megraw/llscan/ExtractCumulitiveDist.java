package osu.megraw.llscan;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

public class ExtractCumulitiveDist {
	private static String pwmFile;
	private static String scoreDir;
	private static String outDir;

	private static final String FWD_SCORE = "FWD.cumscores";
	private static final String REV_SCORE = "REV.cumscores";
	private static final String FWD_LOCS = "FWD.locs";
	private static final String REV_LOCS = "REV.locs";

	public static void main(String[] args) {
		// pwmFile = args[0];
		// scoreDir = args[1];
		// outDir = args[2];
		pwmFile = "/Users/mitra/Downloads/mitra";
		scoreDir = "/Users/mitra/Downloads/mitra";
		outDir = "";
		loadScoresLocs();
	}

	private static void loadScoresLocs() {
		File folder = new File(scoreDir);

		File[] fwdScoreFiles = folder.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(FWD_SCORE);
			}
		});

		File[] fwdLocFiles = folder.listFiles(new FilenameFilter() {

			@Override
			public boolean accept(File dir, String name) {
				return name.endsWith(FWD_LOCS);
			}
		});

		// HAsh <pwm_name, locs_score_hash>
		HashMap<String, String[]> scoreHash = new HashMap<>();
		HashMap<String, String[]> locsHash = new HashMap<>();
		
		for (File fwdScoreFile : fwdScoreFiles) {
			String scoreFilePath = fwdScoreFile.getAbsolutePath();
			System.out.println(scoreFilePath);
			String locsFilePath = scoreFilePath.substring(0, scoreFilePath.indexOf(FWD_SCORE)) + FWD_LOCS;
			System.out.println(locsFilePath);
			try {
				String line = null;
				FileReader scoreReader = new FileReader(scoreFilePath);
				FileReader locsReader = new FileReader(locsFilePath);

				BufferedReader scoresBufferedReader = new BufferedReader(scoreReader);

				while ((line = scoresBufferedReader.readLine()) != null) {
					String pwmName = line.split("\t")[0];
					String key = fwdScoreFile.getName() + "=" + pwmName;
					String[] scores = line.substring(pwmName.length() + 1).split("\t");
					scoreHash.put(key, scores);
				}
				scoresBufferedReader.close();
				
				BufferedReader locsBufferedReader = new BufferedReader(locsReader);

				while ((line = locsBufferedReader.readLine()) != null) {
					String pwmName = line.split("\t")[0];
					String key = fwdScoreFile.getName() + "=" + pwmName;
					String[] locs = line.substring(pwmName.length() + 1).split("\t");
					locsHash.put(key, locs);
				}
				locsBufferedReader.close();

			} catch (Exception e) {
				e.printStackTrace();
			}

		}
		
		HashMap<String, Double[]> fwdHash = new HashMap<>();
		Iterator<String> scoreKeys = scoreHash.keySet().iterator();
		while (scoreKeys.hasNext()) {
			String scoreKey = (String) scoreKeys.next();
			String pwmName = scoreKey.split("=")[1];
			String[] scores = scoreHash.get(scoreKey);
			Double[] item = fwdHash.get(pwmName);
			if (item != null) {
				for (Double score : item) {
					
				}
			}
				
		}
	}
}
