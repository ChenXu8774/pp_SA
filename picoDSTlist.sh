#!/bin/sh
rm -f Temp_picoDST.list
for i in `cat goodrun.list`
do
	find /gpfs/mnt/gpfs02/phenix/fvtx/subsys/fvtx/xuanli/Taxi/Run15pp200MuonsMUPro105/8194/ -name "$i.root" >> Temp_picoDST.list
done

sort Temp_picoDST.list >> picoDST.list

rm -f Temp_picoDST.list

