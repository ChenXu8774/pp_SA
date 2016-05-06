#!/bin/csh -f


find /gpfs/mnt/gpfs02/phenix/spin/spin2/chenxu/spinAnalyzer/pAu_SA/scan/*/ -name "4*.root" >> tmp.lst

sort  tmp.lst >> diMuon_his.lst

rm -rf tmp.lst
