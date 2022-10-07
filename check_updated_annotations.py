import json
import os
from add_substrates_to_json import get_substrates
from sys import argv


def check_updated_annotations(json_dir, substrate_sheet):
    bgc_to_substrates = get_substrates(substrate_sheet)
    minimal_entries = []
    for json_file in sorted(list(os.listdir(json_dir))):
        if json_file.endswith('.json'):
            json_path = os.path.join(json_dir, json_file)
            bgc = json_file.split('.')[0]
            if bgc in bgc_to_substrates:
                with open(json_path, 'r') as json_data:
                    bgc_data = json.load(json_data)
                    if 'nrp' in bgc_data['cluster']:
                        if 'nrps_genes' in bgc_data['cluster']['nrp']:
                            for gene in bgc_data['cluster']['nrp']['nrps_genes']:
                                gene_id = gene['gene_id']
                                a_dom_modules = []
                                if 'modules' in gene:
                                    for module in gene['modules']:
                                        if 'a_substr_spec' in module:
                                            a_dom_modules.append(module)
                                if a_dom_modules and gene_id not in bgc_to_substrates[bgc]:
                                    print(f"Gene {gene_id} missing from {bgc}.")
                                elif a_dom_modules and len(a_dom_modules) > len(bgc_to_substrates[bgc][gene_id]):
                                    print(f"More modules in MIBiG than annotated for {gene_id} in {bgc}.")
                                elif a_dom_modules and len(a_dom_modules) < len(bgc_to_substrates[bgc][gene_id]):
                                    print(f"Fewer modules in MIBiG than annotated for {gene_id} in {bgc}.")
                                elif not a_dom_modules and gene_id in bgc_to_substrates[bgc]:
                                    no_substrates_recorded = True
                                    for substrate in bgc_to_substrates[bgc][gene_id]:
                                        if substrate.specificities:
                                            print(substrate.specificities)
                                            no_substrates_recorded = False
                                            break
                                    if not no_substrates_recorded:
                                        print(f"No adenylation domain modules recorded for {gene_id} in {bgc}.")

                                if gene_id in bgc_to_substrates[bgc]:
                                    if len(a_dom_modules) == len(bgc_to_substrates[bgc][gene_id]):
                                        for substrate in bgc_to_substrates[bgc][gene_id]:
                                            if not substrate.specificities:
                                                pass
                                                # print(f"Missing substrate specificity for domain {substrate.domain_nr} in {gene_id} in {bgc}.")
                    else:
                        minimal_entries.append(bgc)
    counter = 0
    for bgc in minimal_entries:
        bgc_has_entry = False
        for protein, substrates in bgc_to_substrates[bgc].items():
            for substrate in substrates:
                if substrate.specificities:
                    bgc_has_entry = True
                    break
            if bgc_has_entry:
                break
        if bgc_has_entry:
            counter += 1
    print(f"Substrates recorded for {counter} clusters.")


if __name__ == "__main__":
    check_updated_annotations(argv[1], argv[2])

