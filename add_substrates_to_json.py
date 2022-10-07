import json
import os
from pprint import pprint
from sys import argv
from shutil import copy

EVIDENCE_CODES = {'structure_inference': "Structure-based inference",
                  "Structure-based inference": "Structure-based inference",
                  'sequence_based_prediction': "Sequence-based prediction",
                  "Sequence-based prediction": "Sequence-based prediction",
                  "homology": "Homology",
                  "atp_ppi_exchange": "ATP-PPi exchange assay",
                  "radio_labelling": "Radio labelling",
                  "feeding_experiments": "Feeding study",
                  "Feeding study": "Feeding study",
                  "knockout": "Knock-out studies",
                  "in_vitro_tailoring_characterisation": "In-vitro experiments",
                  "in_vitro_reconstitution": "In-vitro experiments",
                  "Activity assay": "Activity assay",
                  "mass_spec": "Mass spectrometry",
                  "heterologous_expression": "Heterologous expression",
                  "cloning": "Heterologous expression",
                  "acvs_assay": "ACVS assay",
                  "crystallography": "X-ray crystallography",
                  "steady_state_kinetics": "Steady-state kinetics",
                  "enzyme_coupled_assay": "Enzyme-coupled assay",
                  "hydroxylamine_trapping_assay": "Enzyme-coupled assay",
                  "hydroxamate_formation_assay": "Enzyme-coupled assay",
                  "nmr": "NMR",
                  "hplc": "HPLC",
                  "None": None}

PROTEINOGENIC = ["Alanine",
                 "Arginine",
                 "Glutamate",
                 "Serine",
                 "Tryptophan",
                 "Methionine",
                 "Threonine",
                 "Glycine",
                 "Isoleucine",
                 "Proline",
                 "Glutamine",
                 "Tyrosine",
                 "Phenylalanine",
                 "Cysteine",
                 "Histidine",
                 "Leucine",
                 "Lysine",
                 "Valine",
                 "Asparagine",
                 "Aspartate",
                 "Aspartic acid",
                 "Glutamic acid"]


def get_structures(smiles_file):
    name_to_smiles = {}
    with open(smiles_file, 'r') as smiles_info:
        smiles_info.readline()
        for line in smiles_info:
            line = line.strip()
            name, smiles = line.split('\t')
            name = name.strip()
            smiles = smiles.strip()
            name_to_smiles[name.lower()] = smiles

    return name_to_smiles


def get_substrates(substrates_file):
    bgc_to_substrates = {}
    with open(substrates_file, 'r') as substrates:
        substrates.readline()
        for line in substrates:
            substrate = Substrate(line)
            if substrate.bgc_id not in bgc_to_substrates:
                bgc_to_substrates[substrate.bgc_id] = {}
            if substrate.protein_id not in bgc_to_substrates[substrate.bgc_id]:
                bgc_to_substrates[substrate.bgc_id][substrate.protein_id] = []

            bgc_to_substrates[substrate.bgc_id][substrate.protein_id].append(substrate)

    for bgc_id, protein_to_substrates in bgc_to_substrates.items():
        for protein, substrates in protein_to_substrates.items():
            substrates.sort(key=lambda x: x.domain_nr)

    return bgc_to_substrates


def add_substrate_evidence(substrate, mibig_substrate):
    substrate_evidence = []
    if substrate.evidence_1 and EVIDENCE_CODES[substrate.evidence_1]:
        substrate_evidence.append(EVIDENCE_CODES[substrate.evidence_1])
    if substrate.evidence_2 and EVIDENCE_CODES[substrate.evidence_2]:
        substrate_evidence.append(EVIDENCE_CODES[substrate.evidence_2])
    if substrate.evidence_3 and EVIDENCE_CODES[substrate.evidence_3]:
        substrate_evidence.append(EVIDENCE_CODES[substrate.evidence_3])

    if substrate_evidence:
        mibig_substrate["evidence"] = substrate_evidence


def add_substrate_publications(substrate, mibig_substrate):
    if substrate.pmid or substrate.doi:
        mibig_substrate["publications"] = []
        if substrate.pmid:
            for pmid in substrate.pmid:
                mibig_substrate["publications"].append(f"pubmed:{pmid}")
        elif substrate.doi:
            for doi in substrate.doi:
                mibig_substrate["publications"].append(f"doi:{doi}")


def get_names_and_structures(substrate, name_to_smiles):
    spec_names = []
    structures = []
    for spec in substrate.specificities:
        spec_names.append(spec)
        structures.append(name_to_smiles[spec.lower()])
    for spec in substrate.minor_specificities:
        spec_names.append(spec)
        structures.append(name_to_smiles[spec.lower()])

    return spec_names, structures


def set_specificity(spec_names, structures, mibig_substrate):
    for j, spec in enumerate(spec_names):

        proteinogenic = False
        name = spec.capitalize()
        if name in PROTEINOGENIC:

            if name == 'Aspartate':
                name = 'Aspartic acid'
            elif name == 'Glutamate':
                name = 'Glutamic acid'
            proteinogenic = True
        elif 'wonky' in name:
            name = name.split(' wonky')[0]
        else:
            name = spec

        structure = structures[j]
        mibig_substrate['substrates'].append({"name": name,
                                              "proteinogenic": proteinogenic,
                                              "structure": structure})


def set_integrated_monomer(substrate, name_to_smiles, mibig_module):
    for spec in substrate.module_specificities:
        structure = name_to_smiles[spec.lower()]
        name = spec.capitalize()
        if name in PROTEINOGENIC:

            if name == 'Aspartate':
                name = 'Aspartic acid'
            elif name == 'Glutamate':
                name = 'Glutamic acid'
        elif 'wonky' in name:
            name = name.split(' wonky')[0]
        else:
            name = spec

        mibig_module['integrated_monomers'].append({"name": name,
                                                    "structure": structure})


def bgc_has_specificities(bgc, bgc_to_substrates):

    for gene_id in bgc_to_substrates[bgc]:
        for substrate in bgc_to_substrates[bgc][gene_id]:
            if substrate.specificities:
                return True

    return False


def gene_has_specificities(gene_id, bgc, bgc_to_substrates):
    for substrate in bgc_to_substrates[bgc][gene_id]:
        if substrate.specificities:
            return True

    return False


def add_genes(bgc, bgc_to_substrates, bgc_data, name_to_smiles):

    for gene_id in bgc_to_substrates[bgc]:
        if gene_has_specificities(gene_id, bgc, bgc_to_substrates):
            gene = {'gene_id': gene_id, 'modules': []}
            for substrate in bgc_to_substrates[bgc][gene_id]:
                module = {'a_substr_spec': {}}
                mibig_substrate = module['a_substr_spec']
                spec_names, structures = get_names_and_structures(substrate, name_to_smiles)
                add_substrate_publications(substrate, mibig_substrate)
                add_substrate_evidence(substrate, mibig_substrate)
                mibig_substrate['substrates'] = []
                set_specificity(spec_names, structures, mibig_substrate)
                if substrate.module_specificities:
                    module['integrated_monomers'] = []
                    set_integrated_monomer(substrate, name_to_smiles, module)

                gene['modules'].append(module)
            bgc_data['cluster']['nrp']['nrps_genes'].append(gene)


def update_changelog(bgc_data, comment):
    entry_to_add_to = None
    version_value = 0
    for i, entry in enumerate(bgc_data['changelog']):
        version = entry['version']
        if float(version) > version_value:
            entry_to_add_to = i
            version_value = float(version)

    if float(version_value) > 2.1:
        if comment not in bgc_data['changelog'][entry_to_add_to]['comments']:
            bgc_data['changelog'][entry_to_add_to]['comments'].append(comment)
    else:
        bgc_data['changelog'].append({'comments': [comment],
                                      'contributors': ["AAAAAAAAAAAAAAAAAAAAAAAA"],
                                      'version': '2.3'})


def add_substrate_to_json(json_dir, substrate_sheet, smiles_sheet, json_out):
    bgc_to_substrates = get_substrates(substrate_sheet)
    name_to_smiles = get_structures(smiles_sheet)
    counter = 0
    counter_2 = 0
    if not os.path.exists(json_out):
        os.mkdir(json_out)
    for json_file in os.listdir(json_dir):

        if json_file.endswith('.json'):

            json_path = os.path.join(json_dir, json_file)
            out_path = os.path.join(json_out, json_file)
            if json_file == 'BGC0000357.json':
                print(f"Skipping {json_file}.")
                copy(json_path, out_path)
                continue
            bgc = json_file.split('.')[0]

            if bgc in bgc_to_substrates:
                with open(json_path, 'r') as json_data:
                    bgc_data = json.load(json_data)
                    if 'nrp' in bgc_data['cluster']:
                        if 'nrps_genes' in bgc_data['cluster']['nrp']:
                            changed_domain = False
                            for gene in bgc_data['cluster']['nrp']['nrps_genes']:
                                gene_id = gene['gene_id']
                                a_dom_modules = []
                                if 'modules' in gene:
                                    for module in gene['modules']:
                                        if 'a_substr_spec' in module:
                                            a_dom_modules.append(module)

                                module_nrs_recorded = True
                                for module in a_dom_modules:
                                    if "module_number" not in module:
                                        module_nrs_recorded = False

                                if module_nrs_recorded:
                                    try:
                                        a_dom_modules.sort(key=lambda x: int(x["module_number"]))
                                    except ValueError:
                                        a_dom_modules.sort(key=lambda x: x["module_number"])
                                    except KeyError:
                                        pass

                                if a_dom_modules and gene_id in bgc_to_substrates[bgc]:

                                    if len(a_dom_modules) == len(bgc_to_substrates[bgc][gene_id]):

                                        for i, module in enumerate(a_dom_modules):
                                            current_module = None
                                            for ori_module in gene['modules']:
                                                if module == ori_module:
                                                    current_module = ori_module
                                                    break

                                            substrate = bgc_to_substrates[bgc][gene_id][i]

                                            if substrate.specificities:

                                                spec_names, structures = get_names_and_structures(substrate, name_to_smiles)

                                                if spec_names:
                                                    mibig_substrate = current_module['a_substr_spec']

                                                    add_substrate_evidence(substrate, mibig_substrate)
                                                    add_substrate_publications(substrate, mibig_substrate)

                                                    mibig_substrate['substrates'] = []
                                                    set_specificity(spec_names, structures, mibig_substrate)
                                                    changed_domain = True

                                                    counter_2 += 1

                                elif not a_dom_modules and gene_id in bgc_to_substrates[bgc]:

                                    specificities_found = False
                                    for substrate in bgc_to_substrates[bgc][gene_id]:
                                        if substrate.specificities:
                                            specificities_found = True

                                    if 'modules' in gene and specificities_found:
                                        if len(bgc_to_substrates[bgc][gene_id]) == len(gene['modules']):
                                            try:
                                                sorted_modules = sorted(gene['modules'], key=lambda x: int(x['module_nr']))
                                            except ValueError:
                                                sorted_modules = sorted(gene['modules'], key=lambda x: x['module_nr'])
                                            except KeyError:
                                                sorted_modules = gene['modules'][:]

                                            for i, substrate in enumerate(bgc_to_substrates[bgc][gene_id]):
                                                mibig_module = sorted_modules[i]
                                                current_module = None
                                                for module in gene['modules']:
                                                    if module == mibig_module:
                                                        current_module = module
                                                        break

                                                assert current_module

                                                if substrate.specificities:

                                                    spec_names, structures = get_names_and_structures(substrate, name_to_smiles)

                                                    if spec_names:
                                                        current_module['a_substr_spec'] = {}
                                                        mibig_substrate = current_module['a_substr_spec']
                                                        add_substrate_evidence(substrate, mibig_substrate)
                                                        add_substrate_publications(substrate, mibig_substrate)

                                                        mibig_substrate['substrates'] = []
                                                        set_specificity(spec_names, structures, mibig_substrate)
                                                        changed_domain = True
                                                        counter_2 += 1

                                        else:
                                            print(bgc, gene_id)
                                            print(gene['modules'])

                            if changed_domain:
                                update_changelog(bgc_data, "Updated NRP substrate specificities.")
                        else:

                            if bgc_has_specificities(bgc, bgc_to_substrates):
                                counter += 1
                                bgc_data['cluster']['nrp']['nrps_genes'] = []
                                add_genes(bgc, bgc_to_substrates, bgc_data, name_to_smiles)
                                update_changelog(bgc_data, "Added NRP substrate specificities.")

                    else:

                        if bgc_has_specificities(bgc, bgc_to_substrates):
                            counter += 1
                            bgc_data['cluster']['nrp'] = {'nrps_genes': []}
                            add_genes(bgc, bgc_to_substrates, bgc_data, name_to_smiles)
                            update_changelog(bgc_data, "Added NRP substrate specificities.")

                    with open(out_path, 'w') as out:
                        json.dump(bgc_data, out, indent=4, ensure_ascii=False)
            else:
                copy(json_path, out_path)

    print(f"Substrate data for {counter} clusters added.")
    print(f"{counter_2} existing substrates updated.")


class Substrate:
    def __init__(self, line):
        line_info = line.split('\t')
        self.annotator = line_info[0]
        self.bgc_id = line_info[1]
        self.protein_id = line_info[2]
        self.uniprot_id = line_info[3]
        self.domain_nr = int(line_info[4])
        specificities = line_info[5].split('|')
        self.specificities = []
        for specificity in specificities:
            if specificity.strip():
                self.specificities.append(specificity.strip())

        module_specificities = line_info[6].split('|')
        self.module_specificities = []
        for specificity in module_specificities:
            if specificity.strip():
                self.module_specificities.append(specificity.strip())
        minor_specificities = line_info[7].split('|')
        self.minor_specificities = []
        for specificity in minor_specificities:
            if specificity.strip():
                self.minor_specificities.append(specificity.strip())
        self.evidence_1 = line_info[8]
        self.evidence_2 = line_info[9]
        self.evidence_3 = line_info[10]
        pmids = line_info[13].split('|')
        dois = line_info[14].split('|')
        self.pmid = []
        self.doi = []
        for pmid in pmids:
            if pmid:
                if pmid.isdigit():
                    self.pmid.append(pmid)
                else:
                    print(f"Invalid PMID: {pmid}")
        for doi in dois:
            if doi:
                if doi.startswith('http'):
                    self.doi.append(doi.split('doi.org/')[1])
                elif doi.startswith('10'):
                    self.doi.append(doi)
                else:
                    print(f"Invalid DOI: {doi}")
        self.has_annotation = line_info[17]
        self.sequence = line_info[18]


if __name__ == "__main__":
    add_substrate_to_json(argv[1], argv[2], argv[3], argv[4])

