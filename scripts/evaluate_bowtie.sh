#!/usr/bin/env bash
cd app
bowtie2 -a --local --mp 2,2 --rdg 10,2 --rfg 10,2 -f -x ./sample_data/reference/reference -U $1 | python ./evaluate.py
