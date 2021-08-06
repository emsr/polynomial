#! /bin/bash

src_dir=`pwd`
echo $src_dir

mkdir -p $HOME/builds/cxx_polynomial/release
cd $HOME/builds/cxx_polynomial/release
cmake -DCMAKE_BUILD_TYPE=release $src_dir
make -j$(nproc)

mkdir -p $HOME/builds/cxx_polynomial/debug
cd $HOME/builds/cxx_polynomial/debug
cmake -DCMAKE_BUILD_TYPE=debug $src_dir
make -j$(nproc)
