virtualenv environments/python2.7
source environments/python2.7/bin/activate
p=`which python`
$p -m pip install --upgrade pip
$p -m pip install pandas
$p -m pip install numpy
$p -m pip install pysam
$p -m pip install networkx
$p -m pip install intervaltree
$p -m pip install biopython
