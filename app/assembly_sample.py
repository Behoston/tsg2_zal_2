from assembly import assembly, DEFAULT_ALGORITHM

SAMPLE_FILES = [
    f'./sample_data/reads_{sample_file}_percent_bad.fasta'
    for sample_file in ['1', '1_to_3', '3_to_5']
]

if __name__ == '__main__':
    for sample_file in SAMPLE_FILES:
        output_file_name = sample_file.split('.')[-2].split('/')[-1]
        assembly(
            input_file_name=sample_file,
            output_file_name=f'./sample_data/output_{output_file_name}.fasta',
            algorithm=DEFAULT_ALGORITHM,
            error_correction=True,
        )
