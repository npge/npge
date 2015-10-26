#!/bin/bash

set -x

sudo apt-get update
sudo apt-get --yes install git > npge-1.log 2>&1
git clone https://github.com/npge/npge \
    > npge-2.log 2>&1
cd npge
git pull
./windows/all-steps-impl.sh
