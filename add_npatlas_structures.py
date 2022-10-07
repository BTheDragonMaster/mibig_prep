import os
import json
from argparse import ArgumentParser
from copy import deepcopy


def main():
    parser = ArgumentParser()
    parser.add_argument("-s", "--structure_data", help="Directory containing structure update information.", type=str, required=True)
    # parser.add_argument("-s", "--smiles", help="Tab-separated file containing NPAtlas SMILES", type=str, required=True)
    # parser.add_argument("-l", "--linking", help="Tab-separated file containing links between NPAtlas and MIBiG", type=str, required=True)
    parser.add_argument("-j", "--json", help="Directory containing json files", type=str, required=True)
    parser.add_argument("-o", "--out", help="Output directory for updated files", type=str, required=True)
    args = parser.parse_args()

    add_structures_to_json(args.json, args.out, args.structure_data)


def make_lowercase(npatlas_name):
    lower = True
    if npatlas_name[0].isalpha() and npatlas_name[0].isupper():
        for char in npatlas_name[1:3]:
            if not char.isalpha() or not char.islower():
                lower = False

    if lower:
        new_npatlas_name = npatlas_name[0].lower() + npatlas_name[1:]
        return new_npatlas_name
    else:
        return npatlas_name


def add_changelog(bgc_data, comment, version="next", contributor="AAAAAAAAAAAAAAAAAAAAAAAA"):
    version_found = False
    changelog_idx = None
    for i, changelog in enumerate(bgc_data["changelog"]):
        if changelog["version"] == version:
            version_found = True
            changelog_idx = i
            break

    if version_found:
        if comment not in bgc_data["changelog"][changelog_idx]["comments"]:
            bgc_data["changelog"][changelog_idx]["comments"].append(comment)
        if contributor not in bgc_data["changelog"][changelog_idx]["contributors"]:
            bgc_data["changelog"][changelog_idx]["contributors"].append(contributor)
    else:
        bgc_data["changelog"].append({"comments": [comment],
                                      "contributors": [contributor],
                                      "version": version})


def add_structures_to_json(json_dir, json_out, structure_data):
    links = Tabular(os.path.join(structure_data, "npatlas_links.txt"), [0, 1])
    npsmiles = Tabular(os.path.join(structure_data, "npatlas_smiles.txt"), [0], lowercase_ids=False)
    compounds_to_add = Tabular(os.path.join(structure_data, "compounds_to_add.txt"), [0, 1])
    name_mapping = Tabular(os.path.join(structure_data, "name_mapping.txt"), [4, 0])
    npaids_to_add = Tabular(os.path.join(structure_data, "npaids_to_add.txt"), [0, 1])
    npentries_to_ignore = Tabular(os.path.join(structure_data, "npatlas_entries_to_ignore.txt"), [0, 1])
    mibig_entries_to_overwrite = Tabular(os.path.join(structure_data, "overwrite_mibig.txt"), [0], lowercase_ids=False)

    bgcs_to_overwrite = mibig_entries_to_overwrite.get_column("bgc_id")
    bgcs_lacking_compounds = set(compounds_to_add.get_column("bgc_id"))
    bgcs_lacking_npaids = set(npaids_to_add.get_column("bgc_id"))

    bgcs = set()

    for compound_id in links.data:
        bgc = links.get_value(compound_id, "bgc_id")
        bgcs.add(bgc)

    nr_compounds_added = 0
    nr_compounds_removed = 0
    nr_compounds_linked = 0
    nr_compounds_renamed = 0
    nr_synonyms_added = 0
    nr_npaids_changed = 0
    nr_structures_changed = 0
    nr_structures_added = 0
    nr_lost_bioactivities = 0

    if not os.path.exists(json_out):
        os.mkdir(json_out)

    for json_file in os.listdir(json_dir):

        if json_file.endswith('.json'):

            json_path = os.path.join(json_dir, json_file)
            out_path = os.path.join(json_out, json_file)

            with open(json_path, 'r') as json_data:
                compounds_added = []
                compounds_removed = []
                compounds_linked = []
                compounds_renamed = []
                synonyms_added = []
                npaids_changed = []
                lost_bioactivities = []
                structures_changed = []
                soft_renamed = False
                structures_added = []

                bgc_data = json.load(json_data)
                original_bgc_data = deepcopy(bgc_data)
                bgc_id = bgc_data["cluster"]["mibig_accession"]
                if bgc_id in bgcs:
                    compounds = bgc_data["cluster"]["compounds"]

                    if bgc_id in bgcs_to_overwrite:
                        for compound in bgc_data["cluster"]["compounds"]:
                            compounds_removed.append(compound["compound"])
                            if "chem_acts" in compound or "chem_targets" in compound:
                                lost_bioactivities.append(compound["compound"])

                        bgc_data["cluster"]["compounds"] = []

                    np_compound_names = set()

                    for compound_id in links.data:
                        if links.get_value(compound_id, 'bgc_id') == bgc_id:
                            np_compound_name = links.get_value(compound_id, 'compound_name')
                            temp_id = npentries_to_ignore.id_separator.join([bgc_id.lower(), np_compound_name.lower()])
                            if temp_id not in npentries_to_ignore.data:
                                np_compound_names.add(np_compound_name)

                    compounds_not_in_mibig = []

                    for np_compound_name in np_compound_names:
                        temp_id = npentries_to_ignore.id_separator.join([bgc_id.lower(), np_compound_name.lower()])
                        npaid = links.get_value(temp_id, "npaid")
                        smiles = npsmiles.get_value(npaid, "compound_smiles")
                        official_np_compound_name = npsmiles.get_value(npaid, "compound_names")
                        use_mibig = True
                        mibig_name = None
                        synonym = ''

                        if name_mapping.id_separator.join([bgc_id.lower(), np_compound_name.lower()]) in name_mapping.data:

                            for name in name_mapping.data:
                                if name_mapping.get_value(name, 'npatlas_name').lower() == np_compound_name.lower():
                                    mibig_name = name_mapping.get_value(name, 'mibig_name')
                                    synonym = name_mapping.get_value(name, 'synonym')
                                    if name_mapping.get_value(name, "use").lower() == "npatlas":
                                        use_mibig = False
                                    break

                            assert mibig_name

                        elif np_compound_name.lower() in [c["compound"].lower() for c in compounds]:
                            mibig_name = np_compound_name

                        else:
                            if bgc_id not in bgcs_to_overwrite:
                                compounds_not_in_mibig.append(np_compound_name)

                        if mibig_name is not None:
                            compound_found = False

                            for compound in compounds:

                                if compound["compound"].lower() == mibig_name.lower():


                                    compound_found = True

                                    if not use_mibig:

                                        if compound["compound"].lower() != official_np_compound_name.lower() and compound["compound"].lower() != np_compound_name:
                                            new_compound_name = make_lowercase(np_compound_name)
                                            compound["compound"] = new_compound_name
                                            compounds_renamed.append((mibig_name, new_compound_name))

                                    if "database_id" in compound:

                                        has_npaid = False
                                        for database_id in compound["database_id"]:
                                            if database_id.startswith('npatlas'):

                                                if npaid != database_id.split(':')[1]:
                                                    print(f"Warning! Overwriting NPAID {database_id.split(':')[1]} with {npaid}")
                                                    npaids_changed.append(compound["compound"])
                                                has_npaid = True
                                        if not has_npaid:
                                            compound["database_id"].append(f"npatlas:{npaid}")
                                            compounds_linked.append(compound["compound"])

                                    else:
                                        compound["database_id"] = [f"npatlas:{npaid}"]
                                        compounds_linked.append(compound["compound"])

                                    if synonym:

                                        if "chem_synonyms" in compound:
                                            if synonym.lower() not in [c.lower() for c in compound["chem_synonyms"]]:
                                                compound["chem_synonyms"].append(synonym)
                                                synonyms_added.append((compound["compound"], synonym))
                                        else:
                                            compound["chem_synonyms"] = [synonym]
                                            synonyms_added.append((compound["compound"], synonym))

                                    if "chem_struct" in compound:
                                        if compound["chem_struct"] != smiles:
                                            compound["chem_struct"] = smiles
                                            structures_changed.append(compound["compound"])
                                    else:
                                        compound["chem_struct"] = smiles
                                        structures_added.append(compound["compound"])

                                    break

                            assert compound_found
                        else:
                            compound = {"compound": make_lowercase(official_np_compound_name),
                                        "chem_struct": smiles,
                                        "database_id": [f"npatlas:{npaid}"]}

                            bgc_data["cluster"]["compounds"].append(compound)
                            compounds_added.append(compound["compound"])
                            compounds_linked.append(compound["compound"])

                    if bgc_id in bgcs_lacking_compounds:
                        for compound_id in compounds_to_add.data:
                            if compounds_to_add.get_value(compound_id, "bgc_id") == bgc_id:
                                compound_name = compounds_to_add.get_value(compound_id, "compound_name")
                                smiles = compounds_to_add.get_value(compound_id, "smiles")
                                npaid = compounds_to_add.get_value(compound_id, "npaid")
                                compound = {"compound": make_lowercase(compound_name)}

                                if smiles:
                                    compound["chem_struct"] = smiles
                                if npaid:
                                    compound["database_id"] = [f"npatlas:{npaid}"]
                                    compounds_linked.append(compound["compound"])

                                compounds_added.append(compound["compound"])
                                bgc_data["cluster"]["compounds"].append(compound)

                    if bgc_id in bgcs_lacking_npaids:
                        for compound_id in npaids_to_add.data:
                            if bgc_id == npaids_to_add.get_value(compound_id, "bgc_id"):
                                compound_name = npaids_to_add.get_value(compound_id, "mibig_compound_name")
                                smiles = npaids_to_add.get_value(compound_id, "smiles")
                                npaid = npaids_to_add.get_value(compound_id, "npaid")
                                npaid_json = f"npatlas:{npaid}"
                                found_compound = False
                                for compound in compounds:
                                    if compound["compound"].lower() == compound_name.lower():
                                        found_compound = True
                                        compound["chem_struct"] = smiles
                                        if "database_id" in compound:
                                            if npaid_json not in compound["database_id"]:
                                                compound["database_id"].append(npaid_json)
                                                compounds_linked.append(compound["compound"])
                                        else:
                                            compound["database_id"] = [npaid_json]
                                            compounds_linked.append(compound["compound"])
                                assert found_compound

                    assert bgc_data["cluster"]["compounds"]
                    nr_compounds_added += len(compounds_added)
                    nr_compounds_linked += len(compounds_linked)
                    nr_compounds_removed += len(compounds_removed)
                    nr_compounds_renamed += len(compounds_renamed)
                    nr_npaids_changed += len(npaids_changed)
                    nr_synonyms_added += len(synonyms_added)
                    nr_structures_changed += len(structures_changed)
                    nr_structures_added += len(structures_added)
                    nr_lost_bioactivities += len(lost_bioactivities)

                    if lost_bioactivities:
                        print(bgc_id)
                        print(lost_bioactivities)

                    if compounds_removed:
                        compounds_removed_string = ', '.join(compounds_removed)
                        add_changelog(bgc_data, f"Removed compounds: {compounds_removed_string}.")

                    if structures_changed:
                        structures_changed_string = ', '.join(structures_changed)
                        add_changelog(bgc_data, f"Changed chemical structures for compounds: {structures_changed_string}.")

                    if structures_added:
                        structures_added_string = ', '.join(structures_added)
                        add_changelog(bgc_data, f"Added chemical structures for compounds: {structures_added_string}.")

                    if compounds_added:
                        compounds_added_string = ', '.join(compounds_added)
                        add_changelog(bgc_data, f"Added compounds: {compounds_added_string}.")

                    if compounds_linked:
                        compounds_linked_string = ', '.join(compounds_linked)
                        add_changelog(bgc_data, f"Added NPAtlas links to compounds: {compounds_linked_string}.")

                    if compounds_renamed:
                        compounds_renamed_string = ', '.join([f"{c[0]} to {c[1]}" for c in compounds_renamed])
                        add_changelog(bgc_data, f"Renamed compounds: {compounds_renamed_string}.")

                    if npaids_changed:
                        npaids_changed_string = ', '.join(npaids_changed)
                        add_changelog(bgc_data, f"Updated NPAtlas ID for compounds: {npaids_changed_string}.")

                    if synonyms_added:
                        synonyms_added_string = ', '.join([f"{c[1]} for {c[0]}" for c in synonyms_added])
                        add_changelog(bgc_data, f"Synonyms added: {synonyms_added_string}.")

                    if bgc_data != original_bgc_data:
                        assert compounds_added or compounds_linked or compounds_removed or npaids_changed or synonyms_added or soft_renamed or compounds_renamed or structures_changed or structures_added
                        with open(out_path, 'w') as out:
                            json.dump(bgc_data, out, indent=4, separators=(',', ': '), sort_keys=True, ensure_ascii=False)

    print(f"Compounds removed: {nr_compounds_removed}.")
    print(f"Compounds added: {nr_compounds_added}.")
    print(f"Compounds linked: {nr_compounds_linked}.")
    print(f"Compounds renamed: {nr_compounds_renamed}.")
    print(f"NPAIDS changed: {nr_npaids_changed}.")
    print(f"Synonyms added: {nr_synonyms_added}.")
    print(f"Structures changed: {nr_structures_changed}.")
    print(f"Structures added: {nr_structures_added}.")
    print(f"Lost bioactivities: {nr_lost_bioactivities}.")


class Tabular:
    def __init__(self, tabular_path, id_columns, separator='\t', id_separator='_', lowercase_ids=True):
        self.data = {}
        self.id_separator = id_separator
        with open(tabular_path, 'r') as tabular_file:
            self.header = tabular_file.readline()
            header = self.header.split('\t')
            self.categories = []
            for i, category in enumerate(header):
                self.categories.append(category.strip())
            for line in tabular_file:
                values = line.split(separator)
                row_id = []
                for id_column in id_columns:
                    if lowercase_ids:
                        row_id.append(values[id_column].strip(" \t\n").lower())
                    else:
                        row_id.append(values[id_column].strip(" \t\n"))

                row_id = self.id_separator.join(row_id)
                if row_id in self.data:
                    print(f"WARNING! Duplicate row ID found in file {tabular_path}: {row_id}")
                self.data[row_id] = {}

                for j, value in enumerate(values):
                    category = self.categories[j]
                    self.data[row_id][category] = value.strip(" '\"\t\n")

    def get_value(self, data_id, category):
        try:
            return self.data[data_id][category]
        except KeyError:
            if data_id in self.data:
                print(f"Cannot find category {category} in data.")
            else:
                print(f"Cannot find id {data_id} in data.")
            raise KeyError

    def get_column(self, category):
        column = []
        for data_id in self.data:
            column.append(self.get_value(data_id, category))

        return column

    def write_table(self, out_file):
        with open(out_file, 'w') as out:
            out.write(self.header)
            for seq_id in self.data:
                for i, category in enumerate(self.categories):
                    if i == len(self.categories) - 1:
                        out.write(f"{self.data[seq_id][category]}\n")
                    else:
                        out.write(f"{self.data[seq_id][category]}\t")


if __name__ == "__main__":
    main()
