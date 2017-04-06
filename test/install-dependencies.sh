#!/bin/bash
set -ex

# samtools 1.4
wget --no-check-certificate https://github.com/samtools/samtools/releases/download/1.4/samtools-1.4.tar.bz2
tar -xvf samtools-1.4.tar.bz2
cd samtools-1.4 && ./configure && make
sudo make PREFIX=/usr install
cd ..

# bcftools 1.4
wget --no-check-certificate https://github.com/samtools/bcftools/releases/download/1.4/bcftools-1.4.tar.bz2
tar -xvf bcftools-1.4.tar.bz2
cd bcftools-1.4 && make
sudo make PREFIX=/usr install
cd ..

# kraken
wget --no-check-certificate https://ccb.jhu.edu/software/kraken/dl/kraken-0.10.5-beta.tgz
tar -xvf kraken-0.10.5-beta.tgz
cd kraken-0.10.5-beta
sudo ./install_kraken.sh /usr/local/bin/
cd ..

# trimmomatic
wget --no-check-certificate http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip
unzip Trimmomatic-0.36.zip
cd Trimmomatic-0.36
printf '#!/bin/bash\njava -jar /home/travis/gopath/src/github.com/will-rowe/gopherSeq/Trimmomatic-0.36/trimmomatic-0.36.jar $@' > trimmomatic
chmod +rx trimmomatic
sudo ln -s $PWD/trimmomatic /usr/local/bin
cd ..

# multiqc
pip install multiqc
