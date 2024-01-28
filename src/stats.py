import numpy as np
from sklearn.manifold import TSNE
import random

def seq_len(seq):
    return len(seq)

def gc_content(seq):
    if len(seq) > 0:
        counter = 0
        for nuc in seq:
            if nuc == 'g' or nuc == 'c':
                counter += 1
        return counter / len(seq)
    return 0

def kmer_counter(kmers):
    kmer_dict = {}

    for kmer in kmers:
        if kmer not in kmer_dict:
            kmer_dict[kmer] = 0
        kmer_dict[kmer] += 1

    return kmer_dict

def dict_to_list(dict, keyname, valuename):
    items = []
    for key in dict:
        value = dict[key]
        items.append({
            keyname: key,
            valuename: value
        })
    return items

def sample(l, sample_size):
    if sample_size <= len(l):
        return random.sample(l, sample_size)
    return l

def proportionally_sample(l, sample_size):
    total_elements = sum(len(arr) for arr in l)
    proportions = [len(arr) / total_elements for arr in l]
    target_elements = [round(sample_size * prop) for prop in proportions]
    sampled_elements = [np.random.choice(arr, size=target, replace=False) for arr, target in zip(l, target_elements)]
    result = np.concatenate(sampled_elements)
    return result[:sample_size]

def eda(data, labels):
    tsne = TSNE(n_components=2)
    data_2d = tsne.fit_transform(data)
    return (data_2d[:,0], data_2d[:,1])

def kmerize(sequence, size):
    seq = sequence.lower()
    return [seq[x:x+size] for x in range(len(seq) - size + 1)]

def jaccard_similarity(set1, set2):
    try:
        intersection = len(set1.intersection(set2))
        union = len(set1.union(set2))
        return intersection / union
    except:
        return 0

def jaccard_distance(sequence1, sequence2, k):
    kmers1 = kmerize(sequence1, k)
    kmers2 = kmerize(sequence2, k)
    similarity = jaccard_similarity(set(kmers1), set(kmers2))
    distance = 1 - similarity
    #print(len(kmers1), len(kmers2), distance)
    return distance

def conduct_pca():
    import pandas as pd    
    # Example dataset with DNA sequences
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

    from sklearn.preprocessing import StandardScaler

    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(df_kmer_counts)

    from sklearn.decomposition import PCA

    n_components = 2
    pca = PCA(n_components=n_components)

    principal_components = pca.fit_transform(standardized_data)

    pc_df = pd.DataFrame(data=principal_components, columns=[f'PC{i}' for i in range(1, n_components + 1)])

    import seaborn as sns
    import matplotlib.pyplot as plt

    # Assuming pc_df is already created
    # pc_df = pd.DataFrame(data=principal_components, columns=[f'PC{i}' for i in range(1, n_components + 1)])

    # Scatter plot
    sns.scatterplot(x='PC1', y='PC2', data=pc_df)
    plt.xlabel('Principal Component 1')
    plt.ylabel('Principal Component 2')
    plt.title('PCA of k-mers')
    plt.show()