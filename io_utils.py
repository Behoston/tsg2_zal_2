def parse_input(input_file_name):
    dane = []
    with open(input_file_name) as f:
        for line in f:
            line=line.strip()
            if line[0]!='>':
                dane.append(line)
    return dane


def dump_output(output_file_name, data):
    with open(output_file_name, 'w') as f:
        f.write(data)
        
