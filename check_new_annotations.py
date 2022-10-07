

def add_substrate_to_json(json_dir, bgc_to_substrates):
    for json_file in os.listdir(json_dir):
        if json_file.endswith('.json'):
            json_path = os.path.join(json_dir, json_file)
            bgc = json_file.split('.')[0]
            if bgc in bgc_to_substrates:
                with open(json_path, 'r') as json_data:
                    bgc_data = json.load(json_data)
                    if 'nrp' in bgc_data['cluster']:
                        if 'nrps_genes' in bgc_data['cluster']['nrp']:
                            for gene in bgc_data['cluster']['nrp']['nrps_genes']:
                                if gene not in bgc_to_substrates[bgc][gene]:
                                    print(f"Gene {gene} missing from BGC.")
                    else:
                        pass
                    pprint(bgc_data)
                    break