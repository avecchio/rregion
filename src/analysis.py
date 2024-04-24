from src.utils import read_json
import glob
import pandas as pd

def organize_ranks():
    datapath = "./work/results/nmf_lstd_counts/"
    files = glob.glob(f"{datapath}nmf_lstd_results_*.json")
    print('Processing: ', len(files))
    all_ranked_features = []

    for file in files:
        print(f'Processing: {file}')
        data = read_json(file)
        region, feature_types = file.replace(datapath, "").replace(".json", "").split("_")[-2:]
        ranked_features = []
        for kmer, rank in data['feature_ranks'].items():
            ranked_features.append({
                'region': region,
                'feature_types': feature_types,
                'kmer': kmer,
                'rank': rank
            })

        sorted_ranked_features = sorted(ranked_features, key=lambda x: x["rank"])
        all_ranked_features += sorted_ranked_features[0:10]

    df = pd.DataFrame(all_ranked_features)
    df.to_csv("nmf_lstd_counts.csv")
    #print(df.head())

