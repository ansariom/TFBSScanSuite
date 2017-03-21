# TFBSScanSuite

## Compile
Create exeucutable jar file by running the make command as follows:
```bash
  make all
```

## Usage

java -jar tfbs_scan.jar GenFeaturesTiledWins --nucsUp 500 --nucsDown 100 --winWidth 100 --mapFile sample_run/map.txt -n 5 sample_run/seqs_2000_region.fa sample_run /pwms.mat sample_run /feature_out.txt
