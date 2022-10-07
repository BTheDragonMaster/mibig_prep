from add_substrates_to_json import Substrate
from sys import argv
from pikachu.general import read_smiles


def parse_smiles_mapping(mapping_file):
    name_to_smiles = {}
    with open(mapping_file, 'r') as mapping:
        for line in mapping:
            line = line.strip()
            name, smiles = line.split('\t')
            name_to_smiles[name] = smiles

    return name_to_smiles


if __name__ == "__main__":
    substrate_sheet = argv[1]
    name_to_smiles = parse_smiles_mapping(argv[2])
    with open(substrate_sheet, 'r') as substrates:
        substrates.readline()
        missed_names = set()
        for line in substrates:
            substrate = Substrate(line)
            names = substrate.specificities + substrate.minor_specificities
            for name in names:
                if name and name not in name_to_smiles:
                    missed_names.add(name)
                elif name:
                    smiles = name_to_smiles[name]
                    read_smiles(smiles)

        for name in missed_names:
            print(name)




