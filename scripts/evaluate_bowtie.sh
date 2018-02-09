#!/usr/bin/env bash
bowtie2 -a --local --mp 2,2 --rdg 10,2 --rfg 10,2 -f -x ./app/sample_data/reference/reference -U $1 | python ./app/evaluate.py
