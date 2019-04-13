# Setup conda environment: python2

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

pip install --upgrade pip

conda install -y \
    python==2.7.15 \
    six==1.12.0 \
    scipy==1.1.0 \
    numpy==1.15.4 \
    matplotlib==2.2.3 \
    pandas==0.23.4 \
    pysam==0.15.1 \
    networkx==2.2 \
    intervaltree==2.1.0 \
    biopython==1.72 \
    Cython==0.29.2 \
    ipython==5.8.0 \
    h5py==2.9.0
