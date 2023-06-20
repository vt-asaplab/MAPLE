#!/bin/bash

# Update system 
sudo apt-get update

# Install building tools
sudo apt-get install -y build-essential

# Install python3 
sudo apt-get install -y python3

# Install ZeroMQ
sudo apt-get install -y libzmq3-dev

# Install EMP-Toolkit
# wget https://raw.githubusercontent.com/emp-toolkit/emp-readme/master/scripts/install.py

# python3 install.py --install --tool --ot --agmpc

# Install libgmp libntl libssl, etc.
sudo apt-get install -y autogen automake build-essential cmake git libboost-dev libboost-thread-dev libgmp-dev libntl-dev libsodium-dev libssl-dev libtool m4 texinfo yasm

# Build MAPLE source code
cd BloomFilterKeywordSearch
make clean
make 
cd ../OMAT/Server
make clean 
make
cd ../Client
make clean 
make
