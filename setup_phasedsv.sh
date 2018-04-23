#!/usr/bin/env bash
export PS1=""
source /home/cmb-16/mjc/shared/virtualenvs/python2.7/bin/activate
source /home/cmb-16/mjc/mchaisso/projects/phasedsv/pbsamstream/setup_pbsamstream.sh
source /usr/usc/java/1.8.0_45/setup.sh 
source /home/cmb-16/mjc/mchaisso/projects/blasr/setup_blasr.sh
BASEDIR=/home/cmb-16/mjc/mchaisso/projects/phasedsv
export PATH=$PATH:$BASEDIR/vt
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/usr/usc/gnu/gcc/5.3.0/lib64/

