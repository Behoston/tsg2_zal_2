import textwrap


def parse_input(input_file_name):
    dane = []
    with open(input_file_name) as f:
        for line in f:
            line = line.strip()
            if line[0] != '>':
                dane.append(line)
    return dane


def dump_output(output_file_name, data: str or []):
    with open(output_file_name, 'w') as f:
        if isinstance(data, str):
            f.write('>\n')
            for line in textwrap.wrap(data, width=80):
                f.write(line)
                f.write('\n')
        else:
            for contig_number, contig in enumerate(data):
                f.write(f'>contig_{contig_number}\n')
                for line in textwrap.wrap(contig, width=80):
                    f.write(line)
                    f.write('\n')
