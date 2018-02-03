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


@click.command()
@click.option('--sample', is_flag=True)
@click.argument('input_file_name', required=False)
@click.argument('output_file_name', required=False, default='output.fasta')
@click.option('--algorithm', required=False, type=click.Choice([key for key in algorithms.keys()]), default='SCS')
@click.pass_context
def assembly_with_sample_option(ctx, sample: bool, **kwargs):
    if not sample:
        assembly()
    else:
        sample_files = ['1', '1_to_3', '3_to_5']
        for sample_file in sample_files:
            ctx.invoke(
                assembly,
                input_file_name=f'./sample_data/reads_{sample_file}_percent_bad.fasta',
                output_file_name=f'./sample_data/output_{sample_file}_percent_bad.fasta',
            )


if __name__ == '__main__':
    assembly_with_sample_option()
