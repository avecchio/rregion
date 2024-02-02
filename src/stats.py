import numpy as np
from sklearn.manifold import TSNE
import random

def seq_len(seq):
    return len(seq)

def gc_content(seq):
    counter = 0
    if len(seq) > 0:
        for nuc in seq.lower():
            if nuc == 'g' or nuc == 'c':
                counter += 1
        return counter / len(seq)
    return counter

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
