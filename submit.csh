#!/bin/csh -f
setenv main_dir `pwd`
setenv i 2001
foreach file_list (`cat picoDST.list`)
  setenv file_l `basename $file_list`
  mkdir -p scan/$i
  chmod 777 scan/$i
  cd scan/$i
  echo $file_list > picoDST.list
  ln -sf $main_dir/condor .
  ln -sf $main_dir/jobscript .
  ln -sf $main_dir/Fun4FVTX_RecoDST_SpinAna.C .
  #./jobscript
  condor_submit condor
  cd $main_dir 
  setenv i `expr 1 + $i`
  end
