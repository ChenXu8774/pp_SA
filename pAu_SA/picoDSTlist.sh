#!/bin/sh
rm -f Temp_picoDST.list
for i in `cat goodrun.list`

do
	find /gpfs/mnt/gpfs02/phenix/spin/spin2/chenxu/taxi/Run15pAu200MuonsMUPro104/7822/data/ -name "$i*.root" >> Temp_picoDST.list
done

sort Temp_picoDST.list >> picoDST.list

rm -f Temp_picoDST.list

