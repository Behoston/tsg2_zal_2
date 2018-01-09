import re

import click
import math
import requests

GLYCINE_MAX_CHR_8 = (
    'https://www.ncbi.nlm.nih.gov/sviewer/viewer.cgi'
    '?tool=portal'
    '&save=file'
    '&log$=seqview'
    '&db=nuccore'
    '&report=fasta'
    '&id=952545308'
    '&maxplex=1'
)


@click.command()
@click.argument('url', required=False, default=GLYCINE_MAX_CHR_8)
@click.option('--overlap', default=8)
def fetch_fasta(url, overlap):
    source_data = download(url)
    split(source_data, overlap)


def download(url):
    source_data = requests.get(url).text
    with open('sample_data/input.fasta', 'w') as f:
        f.write(source_data)
    return source_data


def split(source_data, overlap):
    base_read_len = 80
    source_data = re.sub('(>.*\n)|(\s)|', '', source_data)
    new_in_actual_read = (base_read_len - overlap)
    with open('sample_data/input.fasta', 'w') as f:
        for i in range(math.ceil(len(source_data) / new_in_actual_read)):
            start = i * new_in_actual_read
            f.write('>\n')
            f.write(source_data[start:start + base_read_len])
            f.write('\n')


if __name__ == '__main__':
    fetch_fasta()
