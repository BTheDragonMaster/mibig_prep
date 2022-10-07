import os
import json
from sys import argv
from add_bioactivity_to_json import add_changelog


def get_structures(smiles_file):
    name_to_smiles = {}
    with open(smiles_file, 'r') as smiles_info:
        smiles_info.readline()
        for line in smiles_info:
            line = line.strip()
            name, smiles = line.split('\t')
            name = name.strip()
            smiles = smiles.strip()
            name_to_smiles[name.lower()] = (smiles, name)

    return name_to_smiles


def fix_substrates(json_dir, smiles_data, json_out):
    name_to_smiles = get_structures(smiles_data)

    if not os.path.exists(json_out):
        os.mkdir(json_out)
    for json_file in os.listdir(json_dir):

        if json_file.endswith('.json'):

            json_path = os.path.join(json_dir, json_file)
            out_path = os.path.join(json_out, json_file)
            bgc = json_file.split('.')[0]

            with open(json_path, 'r') as json_data:
                bgc_data = json.load(json_data)
                if 'nrp' in bgc_data['cluster']:
                    if 'nrps_genes' in bgc_data['cluster']['nrp']:
                        for gene in bgc_data['cluster']['nrp']['nrps_genes']:
                            if 'modules' in gene:
                                for module in gene['modules']:
                                    if 'a_substr_spec' in module:
                                        if 'substrates' in module['a_substr_spec'] and not module['a_substr_spec']['substrates']:
                                            del module['a_substr_spec']['substrates']
                                            print(f"Removed empty substrates from {bgc}")
                                        if not module['a_substr_spec']:
                                            del module['a_substr_spec']
                                            print(f"Removed empty substrate spec from {bgc}")

                                    if 'a_substr_spec' in module:
                                        if 'substrates' in module['a_substr_spec']:
                                            for substrate in module['a_substr_spec']['substrates']:
                                                smiles, name = name_to_smiles[substrate['name'].lower()]
                                                substrate['name'] = name
                                    if 'integrated_monomers' in module:
                                        for integrated in module['integrated_monomers']:
                                            smiles, name = name_to_smiles[integrated['name'].lower()]
                                            integrated['name'] = name

                with open(out_path, 'w') as out:
                    json.dump(bgc_data, out, indent=4, ensure_ascii=False)

if __name__ == "__main__":
    fix_substrates(argv[1], argv[2], argv[3])