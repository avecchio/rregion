from sklearn.decomposition import NMF
import pandas as pd
import numpy as np
import json
from itertools import combinations
from itertools import product
from sklearn.decomposition import PCA
from src.utils import read_json, write_json
import matplotlib.pyplot as plt
import glob
import time
import random

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
    work_dir = ".\\work\\stats\\"
    # work_dir = "./work/stats/"
    stats_files = glob.glob(f"{work_dir}*.json")
    file_org_dict = {}

    kmer_lengths = [1, 2, 3, 4, 5, 6, 7]
    kmer_combos = generate_combinations(kmer_lengths)

    for stats_file in stats_files:
        organism, region, file_idx, file_type, file_ext = stats_file.replace(work_dir, "").split(".")
        if organism not in file_org_dict:
            file_org_dict[organism] = {}
        if region not in file_org_dict[organism]:
            file_org_dict[organism][region] = []
        file_org_dict[organism][region].append(stats_file)
        
    for organism in file_org_dict:
        for region in file_org_dict[organism]:
            if organism == 'human' and region in ['promoters']:
                stats_files = file_org_dict[organism][region]
                #for stats_file in stats_files:
                for kmer_combo in kmer_combos:
                    data = []
                    combo_name = '.'.join([str(k) for k in list(kmer_combo)])
                    for stat_file in stats_files:
                        elements = read_json(stat_file)
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
                        #print(organism, region, combo_name)

                    original_length = len(data)
                    if region in ['enhancers', 'silencers', 'promoters']:
                        data = random.choice(data, k=50000)
                    sampled_elements = len(data)
                    write_json(f'{}_{}')
                    # , original_length, sampled_elements
                    print('creating dataframe', original_length, sampled_elements)
                    df = pd.DataFrame(data)
                    df_filled = df.fillna(0)
                    data_matrix = df_filled

                    #pca = PCA(n_components=2)
                    #pca_result = pca.fit_transform(data_matrix)
                    '''
                    # Visualize PCA plot
                    plt.figure(figsize=(8, 6))
                    plt.scatter(pca_result[:, 0], pca_result[:, 1])
                    plt.title('PCA Plot of k-mer Counts')
                    plt.xlabel('Principal Component 1')
                    plt.ylabel('Principal Component 2')
                    for i, seq in enumerate(sequences):
                        plt.annotate(seq, (pca_result[i, 0], pca_result[i, 1]))
                    plt.grid(True)
                    plt.savefig(f'nmf_pca_plot_{organism}_{region}_{combo_name}.png')
                    '''

                    parameters = {
                        #'num_components': range(5, 150),
                        'num_components': list(range(2, 750)),
                        'init': ['random', 'nndsvd', 'nndsvda'],
                        'solver': ['cu','mu'],
                        'beta_loss': ['frobenius', 'kullback-leibler', 'itakura-saito'],
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
                        counter += 1
                        parameter_dict = {}
                        for index, item in enumerate(list(combination)):
                            parameter_key = parameter_keys[index]
                            parameter_dict[parameter_key] = item

                        num_components = parameter_dict['num_components']
                        init = parameter_dict['init']
                        solver = parameter_dict['solver']
                        beta_loss = parameter_dict['beta_loss']
                        try:
                            nmf_model = NMF(
                                n_components = num_components,
                                init=init,
                                solver=solver,
                                beta_loss=beta_loss,
                                random_state=42
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
                            print(organism, region, combo_name, f'{counter}/{total}', rsme_score)
                        except Exception as e:
                            pass
                            #print(e)
                    #print('final')
                    end = time.time()
                    sorted_results = sorted(complete_scores, key=lambda x:x['rsme_score'])
                    #print(sorted_results[0], sorted_results[-1])
                    print(successes, total)

                    with open(f'best_params_{organism}_{region}_{combo_name}.json', 'w') as jf:
                        jf.write(json.dumps({
                            'best_results': sorted_results[0],
                            'worst_results': sorted_results[-1],
                            'time': end - start
                        }))

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