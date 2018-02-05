import tempfile
from contextlib import contextmanager

import pytest
import sh

from assembly import SAMPLE_FILES, _assembly, DEFAULT_ALGORITHM
from evaluate import evaluate


@contextmanager
def pysam_sucks(bowtie_output: sh.RunningCommand):
    """
       At this moment, pysam is stupid and DOES NOT SUPPORT StringIO,
       so we must dump result to file
    """
    with tempfile.NamedTemporaryFile(suffix='.fasta') as file:
        file.write(bowtie_output.stdout)
        file.seek(0)
        yield file


def test_reference():
    indexed_result = sh.bowtie2(
        "-a",
        "--local",
        "--mp 2,2",
        "--rdg 10,2",
        "--rfg 10,2",
        "-f",
        "-x ./sample_data/reference/reference",
        "-U ./sample_data/reference/reference.fasta",
    )
    with pysam_sucks(indexed_result) as file:
        evaluation_results = evaluate(file)
    assert evaluation_results.overall_score == 1


@pytest.mark.parametrize('file', SAMPLE_FILES)
def test_sample_data(file):
    with tempfile.NamedTemporaryFile(suffix='.fasta') as tmp_file:
        _assembly(
            file,
            tmp_file.name,
            DEFAULT_ALGORITHM,
        )
        indexed_result = sh.bowtie2(
            "-a",
            "--local",
            "--mp 2,2",
            "--rdg 10,2",
            "--rfg 10,2",
            "-f",
            "-x ./sample_data/reference/reference",
            f"-U {tmp_file.name}",
        )
        with pysam_sucks(indexed_result) as file:
            evaluation_results = evaluate(file)
        print(evaluation_results)
        assert evaluation_results.overall_score >= 0.01
