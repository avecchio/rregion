import numpy as np
from sklearn.manifold import TSNE
import random
import glob
from itertools import combinations
from src.utils import read_json, write_json
from multiprocessing import Pool
import numpy as np

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

def normalize_kmer_counts(kmer_counts, sequence_length):
    normalized_counts = {}
    mode = 'log10'
    for kmer, count in kmer_counts.items():
        if mode == 'len_std':
            normalized_counts[kmer] = (count / sequence_length)
        elif mode == 'log':
            if count > 0:
                normalized_counts[kmer] = np.log(count)
            else:
                normalized_counts[kmer] = count
        elif mode == 'log10':
            if count > 0:
                normalized_counts[kmer] = np.log10(count)
            else:
                normalized_counts[kmer] = count
        elif mode == 'log_len_std':
            if count > 0:
                normalized_counts[kmer] = np.log(count / sequence_length)
            else:
                normalized_counts[kmer] = count
    return normalized_counts

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

def generate_combinations(lst, min, max):
    all_combinations = []
    for r in range(min, max+1):  # Generating combinations of lengths 1, 2, and 3
        all_combinations.extend(combinations(lst, r))
    return all_combinations

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
random.seed(42)

def sample_elements(data):
    store_path, organism, region, max_sample_size, stats_files = data
    total_elements = 0
    counter = 0
    total_files = len(stats_files)
    for stat_file in stats_files:
        counter += 1
        print(f'Counting {region} {counter}/{total_files}')
        elements = read_json(stat_file)
        total_elements += len(elements)

    if total_elements < max_sample_size:
        counter = 0
        total_files = len(stats_files)
        data = []
        for stat_file in stats_files:
            counter += 1
            print(f'Reading {region} {counter}/{total_files}')
            elements = read_json(stat_file)
            data += elements
        write_json(f'{store_path}sampled_{organism}_{region}_{max_sample_size}.json', data)

    else:
        data = []
        counter = 0
        sample_frac = max_sample_size/total_elements
        for stat_file in stats_files:
            counter += 1
            print(f'Sampling {region} {counter}/{total_files} {sample_frac}')
            elements = read_json(stat_file)
            print(len(elements))
            sampled_elements = random.sample(elements, int(len(elements) * sample_frac))
            data += sampled_elements
        write_json(f'{store_path}sampled_{organism}_{region}_{max_sample_size}.json', data)

def sample_regions(max_sample_size):
    sequences_dir = ".\\work\\sequences\\"
    samples_dir = ".\\work\\samples\\"
    # work_dir = "./work/stats/"
    stats_files = glob.glob(f"{sequences_dir}*.json")
    print(len(stats_files))
    file_org_dict = {}

    for stats_file in stats_files:
        if 'human' in stats_file:
            print(stats_file.replace(sequences_dir, "").split("."))
            organism, region, file_idx, file_type, file_ext = stats_file.replace(sequences_dir, "").split(".")
            if organism not in file_org_dict:
                file_org_dict[organism] = {}
            if region not in file_org_dict[organism]:
                file_org_dict[organism][region] = []
            file_org_dict[organism][region].append(stats_file)
        
    packets = []
    for organism in file_org_dict:
        for region in file_org_dict[organism]:
            if organism == 'human' and region != 'nrcna':
                stats_files = file_org_dict[organism][region]
                packets.append((samples_dir, organism, region, max_sample_size, stats_files))
                print(samples_dir, organism, region)

    with Pool(4) as p:
        p.map(sample_elements, packets)
    #for key, value in file_org_dict['human'].items():
    #    print(key, len(value))
