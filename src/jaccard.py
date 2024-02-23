import random
from src.extract import read_fasta_file

from src.stats import jaccard_distance
from src.utils import write_json

def sample_items(l, sample_size):
    if len(l) <= sample_size:
        return l
    return random.sample()

def sparse_clustering(sequence_references, sample_size, k):
    group_counter = 1
    clusters = {
        'group1': [sequence_references[0]]
    }
    for sequence_reference in sequence_references[1:]:
        sequence = read_fasta_file()
        for group in clusters:
            been_added = False
            references = clusters[group]
            samples = sample_items(references, sample_size)
            distances = []
            for s in samples:
                sample_sequence = read_fasta_file
                dist = jaccard_distance(sequence, sample_sequence, k)
                distances.append(dist)
            dist_avg = sum(distances) / len(distances)
            if dist_avg >= 0.95:
                clusters[group].append(sequence_reference)
                been_added = True
            if (been_added == False):
                group_counter += 1
                clusters[str(f'group{group_counter}')] = [sequence_reference]
    write_json(f'', clusters)
