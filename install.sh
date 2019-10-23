#!/bin/sh


#Welcome to Miniconda3 4.6.14
# Creatinf directory for CLAIR
work_dir="$PWD"
clair_dir="$PWD/bin"
modules_dir="$PWD/bin/modules"
mkdir -p $clair_dir $modules_dir

cd $clair_dir
git clone git@github.com:HKU-BAL/Clair.git
cd $work_dir
cd $modules_dir
wget -v -r  --no-parent --reject "index.html*" http://www.bio8.cs.hku.hk/clair/models/
mv  www.bio8.cs.hku.hk/clair/models/* ./
rm -r -f www.bio8.cs.hku.hk
cd $work_dir


cd $clair_dir
wget https://github.com/squeaky-pl/portable-pypy/releases/download/pypy-7.2.0/pypy-7.2.0-linux_x86_64-portable.tar.bz2
tar -xvjf pypy-7.2.0-linux_x86_64-portable.tar.bz2
cd $work_dir
export PATH=$clair_dir/pypy-7.2.0-linux_x86_64-portable/bin:$PATH
wget https://bootstrap.pypa.io/get-pip.py

pypy get-pip.py
pypy -m pip install --no-cache-dir intervaltree==2.1.0

echo "Done"
