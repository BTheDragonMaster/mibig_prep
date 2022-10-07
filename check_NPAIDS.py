from sys import argv


def get_npaids(npadb):
    npaids = set()
    with open(npadb, 'r') as npa_database:
        npa_database.readline()
        for line in npa_database:
            if line.strip():
                npaid = line.split('\t')[0]
                npaids.add(npaid)

    return npaids


def check_npaids(npadb, structure_sheet):
    npaids = get_npaids(npadb)
    has_structure = set()
    missing_npaids = set()
    with open(structure_sheet, 'r') as structures:
        structures.readline()
        for line in structures:
            if line.strip():
                npaid = line.split('\t')[2]
                structure = line.split('\t')[3]
                if structure.strip():
                    has_structure.add(npaid)

    with open(structure_sheet, 'r') as structures:
        structures.readline()
        for line in structures:
            if line.strip():
                npaid = line.split('\t')[2]
                if npaid not in npaids and npaid not in has_structure:
                    missing_npaids.add(npaid)

    for npaid in missing_npaids:
        print(npaid)


if __name__ == "__main__":
    check_npaids(argv[1], argv[2])
