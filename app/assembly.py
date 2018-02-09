import click

from algorithms import algorithms
from io_utils import dump_output
from io_utils import parse_input
from algorithms.error_corrections import CorrectedReads

DEFAULT_ALGORITHM = 'OLC'
SAMPLE_FILES = [
    f'./sample_data/reads_{sample_file}_percent_bad.fasta'
    for sample_file in ['1', '1_to_3', '3_to_5']
]


@click.command()
@click.argument('input_file_name')
@click.argument('output_file_name', required=False, default='output.fasta')
@click.option('--algorithm', required=False, type=click.Choice([key for key in algorithms.keys()]),
              default=DEFAULT_ALGORITHM)
@click.option('--error_correction', is_flag=True)
def assembly(input_file_name, output_file_name, algorithm, error_correction):
    return _assembly(input_file_name, output_file_name, algorithm, error_correction)


def _assembly(input_file_name, output_file_name, algorithm, error_correction):
    data = parse_input(input_file_name)
    if error_correction:
        data = CorrectedReads(data)
    do_assembly = algorithms[algorithm]
    result = do_assembly(data)
    dump_output(output_file_name, result)


@click.command()
@click.option('--error_correction', is_flag=True)
@click.option('--sample', is_flag=True)
@click.argument('input_file_name', required=False)
@click.argument('output_file_name', required=False, default='output.fasta')
@click.option('--algorithm', required=False, type=click.Choice([key for key in algorithms.keys()]),
              default=DEFAULT_ALGORITHM)
@click.pass_context
def assembly_with_sample_option(ctx, sample: bool, **kwargs):
    if not sample:
        assembly()
    else:
        for sample_file in SAMPLE_FILES:
            _assembly(
                input_file_name=sample_file,
                output_file_name=f'./sample_data/output_{sample_file}_percent_bad.fasta',
                algorithm=DEFAULT_ALGORITHM,
                error_correction=True,
            )


if __name__ == '__main__':
    assembly_with_sample_option()
