import os
import json
from sys import argv


def merge_changelogs(json_dir, json_out):

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
                highest_version = 0.0
                highest_changelog_idx = None
                for i, changelog in enumerate(changelogs):
                    if float(changelog["version"]) > highest_version:
                        highest_changelog_idx = i
                        highest_version = float(changelog["version"])

                if highest_version == 3.0:
                    changelog_to_update = changelogs[highest_changelog_idx]
                else:
                    changelogs.append({"comments": [],
                                       "contributors": [],
                                       "version": "3.0"})
                    changelog_to_update = changelogs[-1]

                    updated = False

                changelogs_to_remove = []

                for i, changelog in enumerate(changelogs[:]):
                    if 3.0 > float(changelog["version"]) > 2.0:
                        updated = True
                        changelogs_to_remove.append(i)
                        for comment in changelog['comments']:
                            if comment not in changelog_to_update["comments"]:
                                changelog_to_update["comments"].append(comment)
                        for contributor in changelog['contributors']:
                            if contributor not in changelog_to_update["contributors"]:
                                changelog_to_update["contributors"].append(contributor)

                changelogs_to_remove.sort(reverse=True)
                for idx in changelogs_to_remove:
                    bgc_data['changelog'].pop(idx)

                if updated:
                    print(bgc_data['changelog'])

                with open(out_path, 'w') as out:
                    json.dump(bgc_data, out, indent=4, ensure_ascii=False)


if __name__ == "__main__":
    merge_changelogs(argv[1], argv[2])

