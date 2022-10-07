import os
import json
from sys import argv


def set_evidence(bgc_data):
    if 'evidence' not in bgc_data['cluster']['loci']:
        bgc_data['cluster']['loci']['evidence'] = []

    if not bgc_data['cluster']['loci']['evidence']:
        bgc_data['cluster']['loci']['evidence'].append('Sequence-based prediction')


def set_minimal(json_dir, json_out):
    counter = 0
    if not os.path.exists(json_out):
        os.mkdir(json_out)
    for json_file in os.listdir(json_dir):

        if json_file.endswith('.json'):

            json_path = os.path.join(json_dir, json_file)
            out_path = os.path.join(json_out, json_file)

            with open(json_path, 'r') as json_data:
                bgc_data = json.load(json_data)
                changelogs = bgc_data['changelog']
                for changelog in changelogs:
                    if float(changelog["version"]) > 2.1 and "Added NRP substrate specificities." in changelog["comments"]:
                        if bgc_data['cluster']['minimal']:
                            counter += 1
                        bgc_data['cluster']['minimal'] = False
                        set_evidence(bgc_data)

                    elif float(changelog["version"]) > 2.1 and "Updated NRP substrate specificities." in changelog["comments"]:
                        if bgc_data['cluster']['minimal']:
                            counter += 1
                        bgc_data['cluster']['minimal'] = False
                        counter += 1
                        set_evidence(bgc_data)

            with open(out_path, 'w') as out:
                json.dump(bgc_data, out, indent=4, ensure_ascii=False)
    print(f"Set {counter} clusters to non-minimal.")


if __name__ == "__main__":
    set_minimal(argv[1], argv[2])