from sklearn.decomposition import NMF
import pandas as pd
import numpy as np
import json
from itertools import product
from sklearn.decomposition import PCA
from src.utils import read_json, write_json
import matplotlib.pyplot as plt
import glob
import time
import random
from src.stats import generate_combinations
from src.plot import lineplot

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

import json
import numpy as np
class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)

def nmf_generate():
    work_dir = ".\\work\\stats\\log\\"
    # work_dir = "./work/stats/"
    stats_files = glob.glob(f"{work_dir}*.json")
    file_org_dict = {}

    kmer_lengths = [1, 2, 3, 4, 5, 6, 7]
    kmer_combos = generate_combinations(kmer_lengths, 1, 3)
    
    for stats_file in stats_files:
        organism, region, file_idx, file_type, file_ext = stats_file.replace(work_dir, "").split(".")

        for kmer_combo in kmer_combos:
            data = []
            combo_name = '.'.join([str(k) for k in list(kmer_combo)])

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

                with open(f'nmf_matrix_{organism}_{region}_{combo_name}.json', 'w') as jf:
                    jf.write(json.dumps({
                        'rsme_score': rsme_score,
                        'parameters': parameters,
                        'H_basis_vectors': basis_vectors,
                        'W_node_embeddings': node_embeddings,
                        'time': (calc_end - calc_start)
                    }, cls=NumpyEncoder))
                print(f'Saved NMF for {organism} {region} {combo_name}')
            except Exception as e:
                print(e)

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
    
def analyze_nmf():

    data_org_dict = {}
    stats_dir = "./work/stats/log/"
    #stats_dir = ".\\work\\stats\\log\\"
    # work_dir = "./work/stats/"
    stats_files = glob.glob(f"{stats_dir}*.json")
    print(stats_files)
    for stats_file in stats_files:
        organism, region, file_idx, file_type, file_ext = stats_file.replace(stats_dir, "").split(".")
        if organism not in data_org_dict:
            data_org_dict[organism] = {}
        if region not in data_org_dict[organism]:
            data_org_dict[organism][region] = {
                'parameter_files': []
            }
        data_org_dict[organism][region]['stats_file'] = stats_file
    print(data_org_dict)
    parameter_files = glob.glob("sampled_best_params*.json")
    for parameter_file in parameter_files:
        organism, region, combo = parameter_file.replace(".json", "").split("_")[3:]
        print(organism)
        data_org_dict[organism][region]['parameter_files'].append(parameter_file)

    for organism in data_org_dict:
        for region in data_org_dict[organism]:
            parameter_files = data_org_dict[organism][region]['parameter_files']
            stats_file = data_org_dict[organism][region]['stats_file']
            #stats_data = read_json(stats_file)
            for parameter_file in parameter_files:
                combo = parameter_file.replace(".json", "").split("_")[-1]
                parameters = read_json(parameter_file)['best_results']
                print()

                for element in elements:
                    element_counter += 1
                    kmer_dict = {}
                    for kmer_length in list(combo.split(".")):
                        element_kmers = element['kmers'][str(kmer_length)]
                        kmer_dict.update(element_kmers)
                    data.append(kmer_dict)
                    df = pd.DataFrame(kmer_data)
                    #'num_components'
                    #'init'
                    #'solver'
                    #'beta_loss'
                    #'kullback-leibler'
                    #print(combo, parameters)

    '''
    parameter_file = parameter_files[0]
    
    '''
    '''
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
                print(f'Processing {e
                lement_counter}/{total_elements}')
                kmer_dict = {}
                for kmer_length in list(kmer_combo):
                    element_kmers = element['kmers'][str(kmer_length)]
                    kmer_dict.update(element_kmers)
                data.append(kmer_dict)
            df = pd.DataFrame(kmer_data)

nmf_model = NMF(
                            n_components = num_components,
                            init=init,
                            solver=solver,
                            beta_loss=beta_loss,
                            random_state=42,
                            max_iter=200
                        )

            W = nmf_model.transform(data_matrix)
            H = nmf_model.components_

            basis_vectors = H
            node_embeddings = W

            df = pd.DataFrame(data)
            column_distances = []

            columns = df.columns

            for i in range(len(columns)):
                for j in range(len(columns)):
                    distance = cosine(df[columns[i]], df[columns[j]])
                    column_distances.append((columns[i], columns[j], distance))

    '''
'''
# Perform grid search with cross-validation
grid_search = GridSearchCV(estimator=nmf_model, param_grid=param_grid, scoring='neg_mean_squared_error', cv=5)
grid_search.fit(data_matrix)  # Assuming 'data_matrix' is your input data matrix

# Get the best model
best_nmf_model = grid_search.best_estimator_

# Fit the best model to your data
best_nmf_model.fit(data_matrix)

# Obtain the optimized factor matrices
basis_vectors = best_nmf_model.components_  # Basis vectors (W)
node_embeddings = best_nmf_model.transform(data_matrix)  # Coefficients (H)

# Compute RMSE for the best model


print("Optimized NMF Model RMSE:", rmse_score)
print("Best Parameters:", grid_search.best_params_)



# Create an array to store the labels of the embedding vectors
embedding_vector_labels = []

# Iterate over the embedding vectors and assign labels based on the associated features
for i, embedding_vector in enumerate(node_embeddings):
    # Get the indices of the features associated with the embedding vector
    feature_indices = feature_embedding_mapping[i]
    # Get the labels of the associated features
    associated_labels = [feature_labels[idx] for idx in feature_indices]
    # Assign the labels to the embedding vector
    embedding_vector_labels.append(associated_labels)

# Print the labels of the embedding vectors
for i, labels in enumerate(embedding_vector_labels):
    print(f"Embedding Vector {i+1} Labels:", labels)




import numpy as np

# Assuming 'basis_vectors' contains the basis vectors (W matrix) obtained from NMF
# 'feature_names' contains the names or labels of the features

# Compute the average magnitude of each basis vector
basis_vector_magnitudes = np.linalg.norm(basis_vectors, axis=1)
# Sort the basis vectors based on their magnitudes (descending order)
sorted_indices = np.argsort(basis_vector_magnitudes)[::-1]

# Get the top ten features from the sorted basis vectors
top_features_indices = sorted_indices[:10]
top_features = [feature_names[i] for i in top_features_indices]

# Print the top ten features
print("Top Ten Features:")
for i, feature in enumerate(top_features, start=1):
    print(f"{i}. {feature}")

import numpy as np
from scipy.spatial.distance import cosine

# Assuming 'basis_vectors' contains the basis vectors (W matrix) obtained from NMF

# Compute the cosine distance between each pair of basis vectors
num_basis_vectors = basis_vectors.shape[0]
cosine_distances = np.zeros((num_basis_vectors, num_basis_vectors))

for i in range(num_basis_vectors):
    for j in range(num_basis_vectors):
        cosine_distances[i, j] = cosine(basis_vectors[i], basis_vectors[j])

# Print the cosine distances
print("Cosine Distances between Basis Vectors:")
for i in range(num_basis_vectors):
    for j in range(num_basis_vectors):
        print(f"({i}, {j}): {cosine_distances[i, j]}")






import numpy as np
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt

# Assuming 'node_embeddings' contains the embedding vectors learned from NMF

# Reduce the dimensionality of the embedding vectors using PCA
pca = PCA(n_components=2)  # Reduce to 2 dimensions for visualization
embedding_2d = pca.fit_transform(node_embeddings)

# Plot the reduced-dimensional embeddings
plt.figure(figsize=(8, 6))
plt.scatter(embedding_2d[:, 0], embedding_2d[:, 1], alpha=0.5)
plt.title('2D Visualization of NMF Embeddings')
plt.xlabel('Principal Component 1')
plt.ylabel('Principal Component 2')
plt.grid(True)
plt.show()
'''