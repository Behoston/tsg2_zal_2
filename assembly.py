import click

from algorithm import do_assembly
from io_utils import dump_output
from io_utils import parse_input


@click.command()
@click.argument('input_file_name')
@click.argument('output_file_name', required=False, default='output.fasta')
def assembly(input_file_name, output_file_name):
    data = parse_input(input_file_name)
    result = do_assembly(data)
    dump_output(output_file_name, result)


if __name__ == '__main__':
    assembly()
