from sys import argv

def parse_mapping(mapping_file):
    old_to_new = {}
    with open(mapping_file, 'r') as mapping:
        for line in mapping:
            line = line.strip()
            old, new = line.split('\t')
            old_to_new[old] = new

    return old_to_new


def replace_ids(mapping_file, substrate_sheet, out_file):
    old_to_new = parse_mapping(mapping_file)
    problematic_ids = set()
    with open(substrate_sheet, 'r') as substrates:
        with open(out_file, 'w') as out:
            out.write(substrates.readline())
            for line in substrates:
                line_info = line.split('\t')
                bgc_id = line_info[1]
                if bgc_id.startswith('NEW'):
                    if bgc_id in old_to_new:
                        line_info[1] = old_to_new[bgc_id]
                        out.write('\t'.join(line_info))
                    else:
                        problematic_ids.add(bgc_id)
                        out.write('\t'.join(line_info))
                elif line.strip():
                    out.write(line)

    for id in problematic_ids:
        print(f"Could not find permanent id for {id}.")


if __name__ == "__main__":
    replace_ids(argv[1], argv[2], argv[3])
