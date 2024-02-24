from sklearn.decomposition import NMF

# Create an instance of the NMF model
num_components = 10  # Number of components (embedding dimensions)
nmf_model = NMF(n_components=num_components, init='random', random_state=42)

# Fit the model to your data
# Assuming 'data_matrix' is your input data matrix (non-negative)
nmf_model.fit(data_matrix)

# Obtain the factor matrices (basis vectors and coefficients)
basis_vectors = nmf_model.components_  # Basis vectors (W)
node_embeddings = nmf_model.transform(data_matrix)  # Coefficients (H)









# Define hyperparameters to tune
param_grid = {
    'n_components': [5, 10, 15],  # Number of components (embedding dimensions)
    'alpha': [0.01, 0.1, 1.0],     # Regularization parameter
    'l1_ratio': [0.0, 0.5, 1.0]     # Sparsity control parameter
}

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
rmse_score = rmse(data_matrix, np.dot(node_embeddings, basis_vectors))

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