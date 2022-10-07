import os
import json
from sys import argv
from shutil import copy


def format_publications(bioactivity, mibig_molecular_target):
    if bioactivity.pmid or bioactivity.doi:
        mibig_molecular_target["publications"] = []
        if bioactivity.pmid:
            for pmid in bioactivity.pmid:
                mibig_molecular_target["publications"].append(f"pubmed:{pmid}")
        elif bioactivity.doi:
            for doi in bioactivity.doi:
                mibig_molecular_target["publications"].append(f"doi:{doi}")


def parse_bioactivities(bioactivity_sheet):
    bgc_to_bioactivities = {}
    with open(bioactivity_sheet, 'r') as bioactivities:
        bioactivities.readline()
        for line in bioactivities:
            bioactivity = Bioactivity(line)
            if bioactivity.bgc_id not in bgc_to_bioactivities:
                bgc_to_bioactivities[bioactivity.bgc_id] = []
            bgc_to_bioactivities[bioactivity.bgc_id].append(bioactivity)

    return bgc_to_bioactivities


def add_changelog(bgc_data, comment, version=3.0, contributor="AAAAAAAAAAAAAAAAAAAAAAAA"):
    version_found = False
    changelog_idx = None
    for i, changelog in enumerate(bgc_data["changelog"]):
        if float(changelog["version"]) == version:
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

def add_bioactivities(json_dir, bioactivity_sheet, json_out):
    bgc_to_bioactivities = parse_bioactivities(bioactivity_sheet)
    counter = 0
    if not os.path.exists(json_out):
        os.mkdir(json_out)

    for json_file in os.listdir(json_dir):

        if json_file.endswith('.json'):

            json_path = os.path.join(json_dir, json_file)
            out_path = os.path.join(json_out, json_file)
            bgc = json_file.split('.')[0]

            if bgc in bgc_to_bioactivities:
                with open(json_path, 'r') as json_data:
                    bgc_data = json.load(json_data)
                    if 'compounds' not in bgc_data['cluster']:
                        print(f"No compounds found for {bgc}.")
                    else:
                        update = False
                        compound_names = []
                        for i, compound in enumerate(bgc_data['cluster']['compounds'][:]):
                            compound_name = compound['compound']
                            for bioactivity in bgc_to_bioactivities[bgc]:
                                annotated_names = []
                                for name in bioactivity.compound_names:
                                    annotated_names.append(name.lower())
                                check_1 = False
                                if compound_name.lower() in annotated_names and bioactivity.activities:
                                    update = True
                                    if 'chem_acts' not in bgc_data['cluster']['compounds'][i]:
                                        bgc_data['cluster']['compounds'][i]['chem_acts'] = []
                                    else:
                                        check_1 = True
                                    for activity in bioactivity.activities:
                                        if activity.lower() not in [x.lower() for x in bgc_data['cluster']['compounds'][i]['chem_acts']]:
                                            bgc_data['cluster']['compounds'][i]['chem_acts'].append(activity)

                                if check_1:
                                    print(bgc_data['cluster']['compounds'][i]['chem_acts'])

                                check = False

                                if compound_name.lower() in annotated_names and bioactivity.molecular_targets:
                                    update = True

                                    if 'chem_targets' not in bgc_data['cluster']['compounds'][i]:
                                        bgc_data['cluster']['compounds'][i]['chem_targets'] = []
                                    else:
                                        check = True
                                    for molecular_target in bioactivity.molecular_targets:
                                        targets = [x['target'] for x in bgc_data['cluster']['compounds'][i]['chem_targets']]
                                        if molecular_target.lower() not in [x.lower() for x in targets]:
                                            target = {"target": molecular_target}
                                            format_publications(bioactivity, target)
                                            bgc_data['cluster']['compounds'][i]['chem_targets'].append(target)
                                if check:
                                    print(bgc_data['cluster']['compounds'][i]['chem_targets'])
                        if update:
                            add_changelog(bgc_data, "Updated bioactivity data")
                            counter += 1

                    with open(out_path, 'w') as out:
                        json.dump(bgc_data, out, indent=4, ensure_ascii=False)

            else:
                copy(json_path, out_path)
    print(f"Activity data added for {counter} BGCs.")


class Bioactivity:
    def __init__(self, line):
        line_info = line.split('\t')
        self.annotator = line_info[0]
        self.bgc_id = line_info[1]
        compound_names = line_info[2].split('|')
        self.compound_names = []
        for compound_name in compound_names:
            if compound_name.strip():
                self.compound_names.append(compound_name.strip())
        self.activities = []
        activities = line_info[3:8]
        for activity in activities:
            if activity.strip():
                self.activities.append(activity)

        self.antibacterial = line_info[8]
        self.molecular_targets = []
        molecular_targets = line_info[9:13]
        for target in molecular_targets:
            if target.strip():
                self.molecular_targets.append(target)

        pmids = line_info[13].split('|')
        dois = line_info[14].split('|')

        self.pmid = []
        self.doi = []
        for pmid in pmids:
            pmid = pmid.strip()
            if pmid:
                if pmid.isdigit():
                    self.pmid.append(pmid)
                else:
                    print(f"Invalid PMID: {pmid}")
        for doi in dois:
            doi = doi.strip()
            if doi:
                if doi.startswith('http'):
                    try:
                        self.doi.append(doi.split('doi.org/')[1])
                    except IndexError:
                        print(f"Invalid DOI: {doi}")
                elif doi.startswith('10'):
                    self.doi.append(doi)
                else:
                    print(f"Invalid DOI: {doi}")


if __name__ == "__main__":
    add_bioactivities(argv[1], argv[2], argv[3])
