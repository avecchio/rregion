from os import write
from src.stats import gc_content, seq_len, kmerize, kmer_counter, jaccard_distance
from src.utils import read_json, write_json, write_csv, mkdir, write_fasta
import pandas as pd
import glob
from multiprocessing import Pool
from itertools import combinations

def f(x):
    return x*x

def nuc_enc(nuc, binarize):
    bin_transform = {'a': [0,0,0,1], 't': [0,0,1,0], 'g': [0,1,0,0], 'c': [1,0,0,0], 'n': [1, 1, 1, 1]}
    dec_transform = {'a': 1, 't': 3, 'g': 7, 'c': 15, 'n': 0}
    if binarize:
        return bin_transform[nuc]
    else:
        return dec_transform[nuc]

def one_hot_enc(seq, binarize=True):
    return [nuc_enc(nuc, binarize) for nuc in seq]

def iupac_clean(seq):
    for char in 'ryswkmbdhv':
        seq = seq.replace(char, 'n')
    return seq

def get_sequence_stats(region_name, organism, sequence, kmer_lengths):
    clean_sequence = iupac_clean(sequence)

    stats = {
        'region_name': region_name,
        'organism': organism,
        'gc_content': gc_content(clean_sequence),
        'len': seq_len(clean_sequence),
        'kmers': {}
    }

    for k in kmer_lengths:
        stats['kmers'][k] = {}
        kmers = kmerize(clean_sequence, k)
        kmer_counts = kmer_counter(kmers)
        for kmer in kmer_counts:
            counts = kmer_counts[kmer]
            stats['kmers'][k][kmer] = counts
    return stats

def exec_jaccard_dist(data_package):
    #print(data_package)
    ref_entry, query_entry, k = data_package
    distance = jaccard_distance(ref_entry['sequence'], query_entry['sequence'], k)
    return {
        'reference_id': ref_entry['id'],
        'query_id': query_entry['id'],
        'distance': distance
    }

def get_region_counts():
    counts_dict = {}
    organism_sequence_files = glob.glob(f"./work/sequences/*.json")
    total = len(organism_sequence_files)
    counter = 0
    for file in organism_sequence_files:
        counts_dict[file] = {}
        counter += 1
        print(counter, total)
        elements = read_json(file)
        for element in elements:
            region_type = element['region']
            organism = element['organism']
            if organism not in counts_dict[file]:
                counts_dict[file][organism] = {}
            if region_type not in counts_dict[file][organism]:
                counts_dict[file][organism][region_type] = 0
            counts_dict[file][organism][region_type] += 1

    array_counts = []
    for file in counts_dict:
        for organism in counts_dict[file]:
            for region_type in counts_dict[file][organism]:
                count = counts_dict[file][organism][region_type]
                array_counts.append({
                    'organism': organism,
                    'file': file,
                    'region': region_type,
                    'count': count
                })

    write_csv(f'count_stats.csv', array_counts)


def calc_all_jaccard_distances():
    sampled_sequence_files = glob.glob(f"./work/samples/*.json")
    region_file_dict = {}
    for sampled_sequence_file in sampled_sequence_files:
        organism, region = sampled_sequence_file.replace("./work/samples/", "").replace(".json", "").split("_")[1:3]
        if region not in region_file_dict:
            region_file_dict[region] = []
        region_file_dict[region].append(sampled_sequence_file)

    print(region_file_dict)
    #for region in region_file_dict:
    #if region_type not in region_file_dict:
    #    print(f'Region {region_type} does not exist')
    #    print(f'Available regions are: ' + '|'.join(list(region_file_dict.keys())))
    for region in region_file_dict:
        sample_paths = []
        for file in region_file_dict[region]:
            file_data = read_json(file)
            organism = file_data['organism']
            region = file_data['region']
            sample_dir = f'./work/jaccard/data/{region}/{organism}/'
            mkdir(sample_dir)
            for sample in file_data['samples']:
                #print(sample)
                sample_id = sample['id']
                sample_fasta_file = f'{sample_dir}/{sample_id}.fasta'
                sample_paths.append(sample_fasta_file.replace("work/jaccard/", "").replace("//", "/"))

                write_fasta(f'{sample_dir}/{sample_id}.fasta', [sample])

        with open(f'./work/jaccard/{region}.filelist.txt', 'w') as f:
            # Convert each element to a string and write to the file
            for sample_path in sample_paths:
                f.write(sample_path + '\n')
    '''
    all_samples = file_data['samples']
    jaccard_distances = []
    counter = 0
    total = len(all_samples) ** 2
    for  in all_samples:
        ref['sequence']
        #mkdir


        for query in all_samples:
            counter += 1
            ref_id = ref['id']
            query_id = query['id']
            print(f'Processing {region} {ref_id}:{query_id} {counter}/{total}')
            distance = jaccard_distance(ref['sequence'], query['sequence'], k)
            jaccard_distances.append({
                'reference_id': ref_id,
                'query_id': query_id,
                'distance': distance
            })

            if counter % 10000000 == 0:

                fname = f'./work/jaccard/jaccard.{region}.{k}.{counter}.json'
                write_json(fname, {
                    'region': region,
                    'distances': jaccard_distances
                })

                jaccard_distances = []

    fname = f'./work/jaccard/jaccard.{region}.{k}.{counter}.json'
    write_json(fname, {
        'region': region,
        'distances': jaccard_distances
    })
    '''

def get_stats(organism):
    organism_sequence_files = glob.glob(f"./work/sequences/{organism}*.json")
    for index, organism_sequence_file in enumerate(organism_sequence_files):
        print(organism_sequence_file.replace("./work/sequences/", "").split("."))
        organism_name, type_iter = organism_sequence_file.replace("./work/sequences/", "").split(".")[0:2]
        sequence_type = type_iter.split("_")[0]
        print('Loading sequences for ' + organism_name)
        sequence_entries = read_json(organism_sequence_file)
        total_entries = len(sequence_entries)
        print(f'Processing {total_entries} entries')
        stats = []
        counter = 0

        for sequence_entry in sequence_entries:
            counter += 1
            #print(f'Processing {counter}/{total_entries} for {organism}')
            region = sequence_entry['region']
            sequence = sequence_entry['sequence']
            sequence_stats = get_sequence_stats(region, organism, sequence, [1, 2, 3, 4])
            stats.append(sequence_stats)
        print('Finished processing... saving')
        #df = pd.DataFrame(stats)
        fname = f'./work/stats/{organism_name}.{sequence_type}.{index}.stats.json'
        write_json(fname, stats)
        #print(fname)
        #df.to_csv(fname, index=False)
        