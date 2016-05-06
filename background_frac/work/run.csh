#!/bin/tcsh -f


#1: binnined fitting; 2: paper plot; 3: test
set run_mode = 1;

if ( $#argv >= 1 ) then
	set run_mode = $1;
endif


if ( $run_mode == 1 ) then
	set out_dir = output;
	mkdir -p $out_dir;
	cd $out_dir;
	rm *;
	touch out.log

	( time root -l -q -b ../Run_gprFit.C\(1,0,0.0,10.,0\) ) >> out.log
	( time root -l -q -b ../Run_gprFit.C\(1,1,0.0,10.,0\) ) >> out.log

	( time root -l -q -b ../Run_gprFit.C\(1,0,0.0,2.0,0\) ) >> out.log
#	( time root -l -q -b ../Run_gprFit.C\(1,0,2.0,4.0,0\) ) >> out.log
	( time root -l -q -b ../Run_gprFit.C\(1,0,2.0,10.,0\) ) >> out.log

	( time root -l -q -b ../Run_gprFit.C\(1,1,0.0,2.0,0\) ) >> out.log
#	( time root -l -q -b ../Run_gprFit.C\(1,1,2.0,4.0,0\) ) >> out.log
	( time root -l -q -b ../Run_gprFit.C\(1,1,2.0,10.,0\) ) >> out.log

	( time root -l -q -b ../Run_gprFit.C\(2,0,1.2,1.8,0\) ) >> out.log
	( time root -l -q -b ../Run_gprFit.C\(2,0,1.8,2.2,0\) ) >> out.log

	( time root -l -q -b ../Run_gprFit.C\(2,1,1.2,1.8,0\) ) >> out.log
	( time root -l -q -b ../Run_gprFit.C\(2,1,1.8,2.2,0\) ) >> out.log

	ln -sf fitparam_fit.dat fitparam.dat
	hadd fit_result.root *_fit.root
endif

if ( $run_mode == 2 ) then
	set out_dir = fitplot;
	mkdir -p $out_dir;
	cd $out_dir;
	rm *;
	touch out.log
	( time root -l -q -b ../Run_gprFit.C\(1,2,0.0,10.,1\) ) >> out.log
	ln -sf fitparam_fit.dat fitparam.dat
	hadd fit_result.root *_fit.root
endif

if ( $run_mode == 3 ) then
	set out_dir = test;
	mkdir -p $out_dir;
	cd $out_dir;
	rm *;
	touch out.log
	( time root -l -q -b ../Run_gprFit.C\(1,0,0.0,2.0,0\) ) >> out.log
	ln -sf fitparam_fit.dat fitparam.dat
	hadd fit_result.root *_fit.root
endif
