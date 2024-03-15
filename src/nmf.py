from sklearn.decomposition import NMF
from itertools import product
from sklearn.decomposition import PCA
from src.utils import read_json, write_json
from src.stats import generate_combinations
from src.plot import lineplot, heatmap, clustermap
from multiprocessing import Pool
from scipy.spatial.distance import cosine
import pandas as pd
import numpy as np
import json
import matplotlib.pyplot as plt
import glob
import time
import random
import os

def rsme(matrix1, matrix2):
    squared_diff = (matrix1-matrix2) ** 2
    mean_squared_diff = np.mean(squared_diff)
    rsme_val = np.sqrt(mean_squared_diff)
    return rsme_val

def count_kmers(sequence, k):
    kmers_count = {}
    for i in range(len(sequence) - k + 1):
        kmer = sequence[i:i+k]
        kmers_count[kmer] = kmers_count.get(kmer, 0) + 1
    return kmers_count

def nmf_regions():
    work_dir = ".\\work\\stats\\log\\"
    # work_dir = "./work/stats/"
    stats_files = glob.glob(f"{work_dir}*.json")
    file_org_dict = {}

    kmer_lengths = [1, 2, 3, 4]
    kmer_combos = generate_combinations(kmer_lengths, 1, 4)
    
    for stats_file in stats_files:
        organism, region, file_idx, file_type, file_ext = stats_file.replace(work_dir, "").split(".")

        for kmer_combo in kmer_combos:
            data = []
            combo_name = '.'.join([str(k) for k in list(kmer_combo)])

            elements = read_json(stats_file)
            elements = random.sample(elements, 1000)
            print(organism, region, combo_name)
            element_counter = 0
            total_elements = len(elements)
            for element in elements:
                element_counter += 1
                print(f'Processing {element_counter}/{total_elements}')
                kmer_dict = {}
                for kmer_length in list(kmer_combo):
                    element_kmers = element['kmers'][str(kmer_length)]
                    kmer_dict.update(element_kmers)
                data.append(kmer_dict)

            print('creating dataframe')
            df = pd.DataFrame(data)
            df_filled = df.fillna(0)
            data_matrix = df_filled

            parameters = {
                #'num_components': range(5, 150),
                'num_components': list(range(20, 1500)),
                'init': ['random', 'nndsvd', 'nndsvda'],
                'solver': ['mu'],
                'beta_loss': ['kullback-leibler', 'itakura-saito'],
            }
            print('constructing parameters')
            values = parameters.values()
            parameter_keys = list(parameters.keys())
            combinations = list(product(*(values)))

            all_parameter_combinations = []
            total = len(combinations)
            counter = 0
            successes = 0
            complete_scores = []
            start = time.time()
            print('now looping through parameter combinations')
            for combination in combinations:
                #print(organism, region, len(data), combination, data[0])
                counter += 1
                parameter_dict = {}
                for index, item in enumerate(list(combination)):
                    parameter_key = parameter_keys[index]
                    parameter_dict[parameter_key] = item

                num_components = parameter_dict['num_components']
                init = parameter_dict['init']
                solver = parameter_dict['solver']
                beta_loss = parameter_dict['beta_loss']
                calc_start = time.time()
                print(f'Starting', combination)
                try:
                    nmf_model = None
                    if solver == 'cd':
                        nmf_model = NMF(
                            n_components = num_components,
                            init=init,
                            solver=solver,
                            #beta_loss=beta_loss,
                            random_state=42,
                            max_iter=200
                        )
                    else:
                        nmf_model = NMF(
                            n_components = num_components,
                            init=init,
                            solver=solver,
                            beta_loss=beta_loss,
                            random_state=42,
                            max_iter=200
                        )
                    nmf_model.fit(data_matrix)

                    W = nmf_model.transform(data_matrix)
                    H = nmf_model.components_

                    basis_vectors = H
                    node_embeddings = W

                    rsme_score = rsme(data_matrix, np.dot(node_embeddings, basis_vectors))
                    parameter_dict['rsme_score'] = rsme_score
                    complete_scores.append(parameter_dict)
                    successes += 1
                    print(organism, region, combo_name, f'{counter}/{total}', rsme_score, combination)
                except Exception as e:
                    #pass
                    print(e)
                with open(f'sampled_saved_params_{organism}_{region}_{combo_name}.json', 'w') as jf:
                    jf.write(json.dumps(complete_scores))
                calc_end = time.time()
                print(organism, region, combo_name, f'{counter}/{total}', (calc_end - calc_start))
            #print('final')
            end = time.time()
            sorted_results = sorted(complete_scores, key=lambda x:x['rsme_score'])
            #print(sorted_results[0], sorted_results[-1])
            print(successes, total)

            with open(f'sampled_best_params_{organism}_{region}_{combo_name}.json', 'w') as jf:
                jf.write(json.dumps({
                    'best_results': sorted_results[0],
                    'worst_results': sorted_results[-1],
                    'time': end - start
                }))

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def conduct_nmf_on_features(data):
    stats_file, kmer_combo, organism, region = data
    data = []
    combo_name = '.'.join([str(k) for k in list(kmer_combo)])
    fname  = f'nmf_matrix_{organism}_{region}_{combo_name}.json'
    if not (os.path.isfile(fname)):
        elements = read_json(stats_file)
        #elements = random.sample(elements, 1000)
        element_counter = 0
        total_elements = len(elements)
        for element in elements:
            element_counter += 1
            print(f'Processing {organism} {region} {element_counter}/{total_elements}')
            kmer_dict = {}
            for kmer_length in list(kmer_combo):
                element_kmers = element['kmers'][str(kmer_length)]
                kmer_dict.update(element_kmers)
            data.append(kmer_dict)
        print(f'Creating dataframe for {organism} {region} {combo_name}')
        calc_start = time.time()

        df = pd.DataFrame(data)
        df_filled = df.fillna(0)
        data_matrix = df_filled

        print(f'Starting NMF for {organism} {region} {combo_name}')
        parameters = {
            'num_components': 1500,
            'init': 'random',
            'solver': 'mu',
            'beta_loss': 'kullback-leibler',
            'random_state': 42,
            'max_iter': 200
        }
        try:
            nmf_model = nmf_model = NMF(
                n_components = parameters['num_components'],
                init = parameters['init'],
                solver = parameters['solver'],
                beta_loss = parameters['beta_loss'],
                random_state = parameters['random_state'],
                max_iter = parameters['max_iter']
            )

            nmf_model.fit(data_matrix)

            W = nmf_model.transform(data_matrix)
            H = nmf_model.components_
            print(f'Finished NMF for {organism} {region} {combo_name}')

            basis_vectors = H
            node_embeddings = W

            rsme_score = rsme(data_matrix, np.dot(node_embeddings, basis_vectors))
            calc_end = time.time()

            with open(fname, 'w') as jf:
                jf.write(json.dumps({
                    'columns': list(df.columns),
                    'rsme_score': rsme_score,
                    'parameters': parameters,
                    'H_basis_vectors': basis_vectors,
                    'W_node_embeddings': node_embeddings,
                    'time': (calc_end - calc_start)
                }, cls=NumpyEncoder))
            print(f'Saved NMF for {organism} {region} {combo_name}')
        except Exception as e:
            print(e)


def nmf_generate():
    work_dir = ".\\work\\stats\\log\\"
    # work_dir = "./work/stats/"
    stats_files = glob.glob(f"{work_dir}*.json")
    file_org_dict = {}

    kmer_lengths = [1, 2, 3, 4, 5] #, 6, 7]
    kmer_combos = generate_combinations(kmer_lengths, 1, 3)
    
    packets = []
    for stats_file in stats_files:
        organism, region, file_idx, file_type, file_ext = stats_file.replace(work_dir, "").split(".")
        for kmer_combo in kmer_combos:
            packets.append((stats_file, kmer_combo, organism, region))

    with Pool(2) as p:
        p.map(conduct_nmf_on_features, packets)

def get_metrics():
    parameter_files = glob.glob("sampled_saved_params*.json")
    all_parameters = []
    init_counts = {}
    beta_loss_counts = {}
    for parameter_file in parameter_files:
        combo = parameter_file.replace(".json", "").split("_")[-1]
        parameters = read_json(parameter_file)
        for parameter in parameters:
            parameter['combo'] = combo
            init = parameter['init']
            if init not in init_counts:
                init_counts[init] = 0
            init_counts[init] += 1

            beta_loss = parameter['beta_loss']
            if beta_loss not in beta_loss_counts:
                beta_loss_counts[beta_loss] = 0
            beta_loss_counts[beta_loss] += 1

            all_parameters.append(parameter)

    print(init_counts)
    print(beta_loss_counts)
    df = pd.DataFrame(all_parameters)

    plot_name = 'test_rsme.lineplot.png'
    lineplot(df, plot_name, {
        'xdata': 'num_components',
        'ydata': 'rsme_score',
        'xlabel': '# Components',
        'ylabel': 'RSME',
        'hue': 'combo',
        'title': '# Components vs RSME'
    })


def analyze_single_nmf(matrix_file):
    analysis_fname = matrix_file.replace("matrix", "results")
    #print(matrix_file.replace(".json", "").split("_"))
    region, kmer_combo = (matrix_file.replace(".json", "").split("_"))[-2:]
    print(region, kmer_combo)

    matrix_data = read_json(matrix_file)
    basis_vectors = matrix_data['H_basis_vectors']
    node_embeddings = matrix_data['W_node_embeddings']

    features = matrix_data['columns']
    df = pd.DataFrame(basis_vectors)
    df.columns = features

    columns = df.columns

    num_columns = len(columns)

    column_magnitudes = np.sum(df, axis=0)

    feature_ranks = np.argsort(-column_magnitudes)

    cosine_distances = np.zeros((num_columns, num_columns))

    for i in range(num_columns):
        for j in range(num_columns):
            cosine_distances[i, j] = cosine(df[columns[i]], df[columns[j]])

    for i in range(num_columns):
        for j in range(num_columns):
            print(f'Distance between {columns[i]} and {columns[j]}:', cosine_distances[i, j])

    distance_df = pd.DataFrame(cosine_distances)
    distance_df.columns = columns
    distance_df['feature'] = columns

    def series_to_json(series):
        return json.loads(pd.Series(series).to_json())

    with open(analysis_fname, 'w') as mf:
        mf.write(json.dumps({
            'column_magnitudes': series_to_json(column_magnitudes),
            'feature_ranks': series_to_json(feature_ranks),
            'cosine_distances': distance_df.to_dict('records')
        }))

    nedf = pd.DataFrame(node_embeddings)

    kmer_title = kmer_combo.replace(".", ",")
    heatmap(nedf, f'./work/analysis/node_embeddings_heatmap/{region}_{kmer_combo}_heatmap.png', {
        'xlabel': 'Topics',
        'ylabel': 'Samples',
        'title': f'NMF Node Embeddings for {region} (k=[{kmer_title}])',
    })

    clustermap(nedf, f'./work/analysis/node_embeddings_heatmap/{region}_{kmer_combo}_clustermap.png', {
        'xlabel': 'Topics',
        'ylabel': 'Samples',
        'title': f'NMF Node Embeddings for {region} (k=[{kmer_title}])',
    })

    distance_df = distance_df.drop('feature', axis=1)
    print(distance_df.columns)
    distance_df.index = distance_df.columns
    heatmap(distance_df, f'./work/analysis/cosine_distances/{region}_{kmer_combo}_heatmap.png', {
        'xlabel': 'Kmers',
        'ylabel': 'Kmers',
        'title': f'Cosine distance of Kmers for {region} (k=[{kmer_title}])',
    })

    clustermap(distance_df, f'./work/analysis/cosine_distances/{region}_{kmer_combo}_clustermap.png', {
        'xlabel': 'Kmers',
        'ylabel': 'Kmers',
        'title': f'Cosine distance of Kmers for {region} (k=[{kmer_title}])',
    })


def analyze_nmf():
    matrix_files = glob.glob("nmf_matrix_human_*")
    with Pool(4) as p:
        p.map(analyze_single_nmf, matrix_files)
