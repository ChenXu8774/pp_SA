#!/bin/csh -f

setenv runlistfile `echo goodrun.list`
setenv runlist `cat $runlistfile`

foreach runnumber ( $runlist )

	setenv fill_n `psql -h phnxdb0.phenix.bnl.gov -Uphnxrc -d daq -c "select * from run where Runnumber= $runnumber and eventsinrun>=-1 and RunType='PHYSICS'" | grep "$runnumber" | awk '{print $60}'`


	echo $runnumber	$fill_n >>runfilltable.list
end


