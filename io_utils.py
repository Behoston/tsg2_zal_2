def parse_input(input_file_name):
    with open(input_file_name) as f:
        for line in f:
            pass


def dump_output(output_file_name, data):
    with open(output_file_name, 'w') as f:
        f.write(data)
