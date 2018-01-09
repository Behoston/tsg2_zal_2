import click

from algorithms import algorithms
from io_utils import dump_output
from io_utils import parse_input


@click.command()
@click.argument('input_file_name')
@click.argument('output_file_name', required=False, default='output.fasta')
@click.option('--algorithm', required=False, type=click.Choice([key for key in algorithms.keys()]), default='SCS')
def assembly(input_file_name, output_file_name, algorithm):
    data = parse_input(input_file_name)
    do_assembly = algorithms[algorithm]
    result = do_assembly(data)
    dump_output(output_file_name, result)


if __name__ == '__main__':
    assembly()
