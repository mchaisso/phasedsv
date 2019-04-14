# Setup conda environment: pacbio

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

pip install --upgrade pip

conda install -y \
    genomicconsensus==2.3.2 \
    pbbam==0.23.0 \
    blasr==5.3.2 \
    pbmm2==1.0.0 \
    bax2bam==0.0.9
