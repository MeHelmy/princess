
#!/bin/sh


# Welcome to Miniconda3 4.6.14
# Creatinf directory for CLAIR
work_dir="$PWD"
clair_dir="$PWD/bin"
modules_dir="$PWD/bin/modules"
mkdir -p $clair_dir $modules_dir

# cd $clair_dir
# git clone git@github.com:HKU-BAL/Clair.git
# cd $work_dir
cd $modules_dir
# download pretrained model (for Illumina)
mkdir illumina && cd illumina
wget http://www.bio8.cs.hku.hk/clair_models/illumina/12345.tar
tar -xf 12345.tar
cd ../

# download pretrained model (for PacBio CCS)
mkdir ccs && cd ccs
wget http://www.bio8.cs.hku.hk/clair_models/pacbio/ccs/15.tar
tar -xf 15.tar
cd ../

# download pretrained model (for PacBio CLR)
mkdir clr && cd clr
wget http://www.bio8.cs.hku.hk/clair_models/pacbio/clr/1234567.tar
tar -xf 1234567.tar
cd ../

# download pretrained model (for ONT)
mkdir ont && cd ont
wget http://www.bio8.cs.hku.hk/clair_models/ont/1234.tar
tar -xf 1234.tar
cd ../
# wget -v -r  --no-parent --reject "index.html*" http://www.bio8.cs.hku.hk/clair/models/
# mv  www.bio8.cs.hku.hk/clair/models/* ./
# rm -r -f www.bio8.cs.hku.hk

cd $clair_dir
wget https://github.com/squeaky-pl/portable-pypy/releases/download/pypy3.6-7.2.0/pypy3.6-7.2.0-linux_x86_64-portable.tar.bz2
mv pypy* pypy.tar.bz2
tar -xvjf pypy.tar.bz2
mv pypy3.* pypy
# tar -xvjf pypy3.6-7.2.0-linux_x86_64-portable.tar.bz2
cd $work_dir
export PATH=$clair_dir/pypy/bin:$PATH


wget https://bootstrap.pypa.io/get-pip.py
python get-pip.py
pypy3 -m ensurepip
#pypy3 -m pip install --no-cache-dir intervaltree blosc
pypy3 -m pip install --no-cache-dir intervaltree
pypy3 -m pip install blosc==1.8.3

echo "Done"
