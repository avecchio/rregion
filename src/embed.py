from sklearn.decomposition import NMF
import numpy as np

# Generate sample data
# Replace this with your actual data
#data = np.random.rand(1000000, 2000)  # 100 samples, 20 features
data = np.random.rand(1300000, 2**7)  # 100 samples, 20 features

# Specify the number of components (i.e., the reduced dimensionality)
for n in range(5, 30):
    n_components = n

    # Initialize the NMF model
    nmf = NMF(n_components=n_components, init='random', random_state=0)

    # Fit the model to the data
    W = nmf.fit_transform(data)  # W contains the transformed data (lower-dimensional representation)
    H = nmf.components_  # H contains the components (basis vectors)

    # Reconstruction of the original data
    data_reconstructed = np.dot(W, H)
    print(data_reconstructed)

    ww = json.dumps(W, cls=NumpyEncoder)
    hh = json.dumps(H, cls=NumpyEncoder)

    with open(f'nnmf.{n}.json', 'w') as jf:
        jf.write(json.dumps({'w': ww, 'h': hh}))