#!/usr/bin/env bash
mkdir -p ~/venvs/
virtualenv -p /usr/bin/python3.6 ~/venvs/tsg_2_2
source ~/venvs/tsg_2_2/bin/activate
pip install --upgrade pip
pip install pip-tools
cd ./app
pip-sync

