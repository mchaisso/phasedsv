virtualenv environments/python2.7 --python=`which python`
source environments/python2.7/bin/activate
p=`which python`
$p -m pip install --no-cache-dir  --upgrade pip
$p -m pip install --no-cache-dir --no-build-isolation pandas 
$p -m pip install --no-cache-dir numpy
$p -m pip install --no-cache-dir pysam
$p -m pip install --no-cache-dir networkx
$p -m pip install --no-cache-dir intervaltree
$p -m pip install --no-cache-dir biopython
