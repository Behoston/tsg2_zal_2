import click

from algorithms import algorithms
from algorithms.error_corrections import CorrectedReads
from io_utils import dump_output
from io_utils import parse_input

DEFAULT_ALGORITHM = 'OLC_NAIVE'


@click.command()
@click.argument('input_file_name')
@click.argument('output_file_name', required=False, default='output.fasta')
@click.option('--algorithm', required=False, default=DEFAULT_ALGORITHM,
              type=click.Choice([key for key in algorithms.keys()]))
@click.option('--no-error_correction', is_flag=True)
def _assembly(input_file_name, output_file_name, algorithm, no_error_correction):
    return assembly(input_file_name, output_file_name, algorithm, error_correction=not no_error_correction)


def assembly(input_file_name, output_file_name, algorithm, error_correction):
    data = parse_input(input_file_name)
    if error_correction:
        data = CorrectedReads(data)
    do_assembly = algorithms[algorithm]
    result = do_assembly(data)
    dump_output(output_file_name, result)


if __name__ == '__main__':
    _assembly()
