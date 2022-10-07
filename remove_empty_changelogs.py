import os
import json
from sys import argv


def remove_empty_changelogs(json_dir, json_out):

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
                changelogs_to_remove = []
                for i, changelog in enumerate(changelogs[:]):
                    if not changelog['comments']:
                        changelogs_to_remove.append(i)

                changelogs_to_remove.sort(reverse=True)
                for idx in changelogs_to_remove:
                    bgc_data['changelog'].pop(idx)

                with open(out_path, 'w') as out:
                    json.dump(bgc_data, out, indent=4, ensure_ascii=False)


if __name__ == "__main__":
    remove_empty_changelogs(argv[1], argv[2])

