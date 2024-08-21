#!/usr/bin/env bash 

set -e 

pip install --upgrade pip 
./install-smoldyn-mac-silicon.sh  
pip install -r requirements.txt 
pip install -e .