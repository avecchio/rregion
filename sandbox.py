'''
import pandas as pd
data = {
    'Seq1': 'ATCGATCGATCG',
    'Seq2': 'CGATCGATCGAT',
    'Seq3': 'GATCGATCGATC',
    # ... add more sequences as needed
}

df = pd.DataFrame(data)

from collections import Counter

def generate_kmer_counts(sequence, k):
    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    return dict(Counter(kmers))

k = 3  # Specify the desired k-mer length
df_kmer_counts = df.applymap(lambda seq: generate_kmer_counts(seq, k)).fillna(0)

print(df_kmer_counts)
'''
'''
import numpy as np
import nimfa

X = np.random.rand(100, 50)

rank = 5

nmtf = nimfa.Nmtf(X, seed="random_vcol", rank=rank, max_iter=100)

nmtf_fit = nmtf()
basis_matrix = nmtf_fit.basis()
coeff_matrix1 = nmtf_fit.coef(idx=0)
coeff_matrix2 = nmtf_fit.coef(idx=1)
'''
'''
import numpy as np
from nimfa import TriFactorization

# Generate sample data
# Replace this with your actual data
data = np.random.rand(100, 20)  # 100 samples, 20 features

# Specify the number of components
n_components = 5

# Initialize the NTF model
model = TriFactorization(data, rank=n_components)

# Fit the model to the data
fit = model()

# Get the factor matrices
W = fit.basis()  # Basis matrix (size: n_samples x n_components)
H = fit.coef()   # Coefficient matrix (size: n_features x n_components)
G = fit.gain()   # Gain matrix (size: n_components x n_components)

# Reconstruction of the original data
data_reconstructed = np.dot(np.dot(W, G), H.T)
'''
import pandas as pd
import numpy as np
from scipy.spatial.distance import cosine
import json

class NumpyEncoder(json.JSONEncoder):
    def default(self, obj):
        if isinstance(obj, np.ndarray):
            return obj.tolist()
        return json.JSONEncoder.default(self, obj)


data = {
    'col1': [1, 2, 3, 4, 5],
    'col2': [6, 7, 8, 9, 10],
    'col3': [11, 12, 13, 14, 15],
    'col4': [16, 17, 18, 19, 20],
    'col5': [21, 22, 23, 24, 25],
}

df = pd.DataFrame(data)

columns = df.columns

num_columns = len(columns)

column_magnitudes = np.sum(df, axis=0)

#print(column_magnitudes)
feature_ranks = np.argsort(-column_magnitudes)
#print(feature_ranks)

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
#distance_df.to_csv(f"sandbox.csv", index=False)

#print(cosine_distances)

def series_to_json(series):
    return json.loads(pd.Series(series).to_json())

with open('nmf_metrics.json', 'w') as mf:
    mf.write(json.dumps({
        'column_magnitudes': series_to_json(column_magnitudes),
        'feature_ranks': series_to_json(feature_ranks),
        'cosine_distances': distance_df.to_dict('records')
    }))
