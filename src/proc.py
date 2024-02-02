from src.utils import read_json, write_json
import pandas as pd
import numpy as np

def down_sample_regions(sample_size):
    counts_df = pd.read_csv("count_stats.csv")

    counts_dict = {}
    for record in counts_df.to_dict('records'):
        region = record['region']
        organism = record['organism']
        file = record['file']
        count = record['count']

        if organism not in counts_dict:
            counts_dict[organism] = {}
        
        if region not in counts_dict[organism]:
            counts_dict[organism][region] = []
        
        counts_dict[organism][region].append({
            'file': file,
            'size': count,
        })

    for organism in counts_dict:
        for region in counts_dict[organism]:
            datasets = counts_dict[organism][region]
            total_elements = 0
            for dataset in datasets:
                total_elements += dataset['size']

            if total_elements <= sample_size:
                total_samples = []
                for dataset in datasets:
                    total_samples += read_json(dataset['file'])

                write_json(f"./work/samples_{organism}_{region}.json", {
                    'organism': organism,
                    'region': region,
                    'samples': total_samples
                })

            else:
                samples = []
                for dataset in datasets:
                    proportion = dataset['size'] / total_elements
                    target_elements = round(sample_size * proportion)
                    print('loading file ' + dataset['file'])
                    data = read_json(dataset['file'])
                    print('sampling')
                    sampled_elements = np.random.choice(data, size=target_elements, replace=False)
                    samples += list(sampled_elements)

                print(f'Finished sampling for {organism} {region}')
                write_json(f"./work/samples_{organism}_{region}.json", {
                    'organism': organism,
                    'region': region,
                    'samples': samples[:sample_size]
                })
