import json
import os
from sys import argv
from shutil import copy


def get_changelog(bgc_data, version=3.0):
    for changelog in bgc_data['changelog']:
        if float(changelog["version"]) == version:
            return changelog


def compare_bioactivity_data(bgc_data_brt, bgc_data_mibig):
    brt_changelog = get_changelog(bgc_data_brt)

    mistake = False

    if brt_changelog:
        if "Updated bioactivity data" in brt_changelog["comments"]:
            changed = False
            if bgc_data_brt['cluster']['compounds'] != bgc_data_mibig['cluster']['compounds']:
                changed = True
            if not changed:
                brt_changelog["comments"].remove("Updated bioactivity data")
                mistake = True

    if brt_changelog and not brt_changelog["comments"]:
        bgc_data_brt['changelog'].remove(brt_changelog)
        print(f"Removed changelog for {bgc_data_brt['cluster']['mibig_accession']}, bioactivity.")
        print('\n')

    return mistake


def compare_substrate_data(bgc_data_brt, bgc_data_mibig):
    brt_changelog = get_changelog(bgc_data_brt)
    mistake = False

    if brt_changelog:
        if "Updated NRP substrate specificities." in brt_changelog["comments"]:
            brt_changelog["comments"].remove("Updated NRP substrate specificities.")
            changed = False
            if bgc_data_brt['cluster']['nrp']['nrps_genes'] != bgc_data_mibig['cluster']['nrp']['nrps_genes']:
                changed = True
                brt_changelog["comments"].append("Updated NRP substrate specificities")
            if not changed:
                mistake = True

    if brt_changelog and not brt_changelog["comments"]:
        bgc_data_brt['changelog'].remove(brt_changelog)
        print(f"Removed changelog for {bgc_data_brt['cluster']['mibig_accession']}, updated substrates.")

    return mistake


def compare_substrate_addition(bgc_data_brt, bgc_data_mibig):
    brt_changelog = get_changelog(bgc_data_brt)
    mistake = False

    if brt_changelog:
        if "Added NRP substrate specificities." in brt_changelog["comments"]:
            brt_changelog["comments"].remove("Added NRP substrate specificities.")
            changed = False
            if 'nrp' in bgc_data_mibig['cluster'] and \
                    'nrps_genes' in bgc_data_mibig['cluster']['nrp'] and \
                    bgc_data_brt['cluster']['nrp']['nrps_genes'] != bgc_data_mibig['cluster']['nrp']['nrps_genes']:
                changed = True
            elif 'nrp' not in bgc_data_mibig['cluster']:
                changed = True
            elif 'nrps_genes' not in bgc_data_mibig['cluster']['nrp']:
                changed = True

            if not changed:

                mistake = True
            else:
                brt_changelog["comments"].append("Added NRP substrate specificities")

    if brt_changelog and not brt_changelog["comments"]:
        bgc_data_brt['changelog'].remove(brt_changelog)
        print(f"Removed changelog for {bgc_data_brt['cluster']['mibig_accession']}, added substrates.")

    return mistake


def correct_changelogs(json_dir_brt, json_dir_mibig, json_out):
    counter = 0
    counter_2 = 0
    if not os.path.exists(json_out):
        os.mkdir(json_out)
    for json_file in os.listdir(json_dir_brt):

        if json_file.endswith('.json'):
            bgc = json_file.split('.')[0]
            brt_path = os.path.join(json_dir_brt, json_file)
            mibig_path = os.path.join(json_dir_mibig, json_file)
            out_path = os.path.join(json_out, json_file)
            if os.path.exists(mibig_path):
                with open(brt_path, 'r') as brt_json:
                    with open(mibig_path, 'r') as mibig_json:
                        bgc_data_brt = json.load(brt_json)
                        bgc_data_mibig = json.load(mibig_json)
                        changelog_removed = compare_bioactivity_data(bgc_data_brt, bgc_data_mibig)
                        if changelog_removed:
                            counter += 1
                        substrate_changelog_removed = compare_substrate_data(bgc_data_brt, bgc_data_mibig)
                        if substrate_changelog_removed:
                            counter_2 += 1
                        added_changelog_removed = compare_substrate_addition(bgc_data_brt, bgc_data_mibig)
                        if added_changelog_removed:
                            counter_2 += 1

                        with open(out_path, 'w') as out:
                            json.dump(bgc_data_brt, out, indent=4, ensure_ascii=False)

            else:
                copy(brt_path, mibig_path)

    print(counter)
    print(counter_2)


if __name__ == "__main__":
    correct_changelogs(argv[1], argv[2], argv[3])