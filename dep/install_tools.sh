# Setup conda environment: tools

conda config --add channels defaults
conda config --add channels bioconda
conda config --add channels conda-forge

pip install --upgrade pip

conda install -y \
    canu==1.8 \
    bedtools==2.27.1 \
    samtools==1.9 \
    tabix==0.2.6 \
    vt==2015.11.10 \
    trf==4.09 \
    repeatmasker==4.0.8 \
    bioawk==1.0 \
    ucsc-bedToBigBed==377
