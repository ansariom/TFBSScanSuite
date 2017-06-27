#!/bin/bash

java -jar tfbs_scan.jar GenFeaturesTiledWins --nucsUp 500 --nucsDown 100 --winWidth 100 --mapFile sample_run/map.txt -n 5 sample_run/seqs_2000_region.fa sample_run/pwms.mat sample_run/feature_out.txt

java -jar tfbs_scan.jar GenFeaturesByNT -N 5000 -B 250 sample/peat_roe_small_fwd.tbl sample/peat_roe_small_rev.tbl sample/small_seq.fa sample/pwm_orig.txt out.txt
