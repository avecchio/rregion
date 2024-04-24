from src.stats import gc_content, seq_len, kmerize, kmer_counter, jaccard_distance, normalize_kmer_counts
from src.utils import chunk_array, read_json, write_json, write_csv, mkdir, write_fasta
import pandas as pd
import glob
import os
from multiprocessing import Pool
from itertools import combinations

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
        kmers = [kmer for kmer in kmers if 'n' not in kmer]
        kmer_counts = normalize_kmer_counts(kmer_counter(kmers), len(sequence))        
        for kmer in kmer_counts:
            counts = kmer_counts[kmer]
            stats['kmers'][k][kmer] = counts
    return stats

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


def get_stats(organism):
    work_path = ".\\work\\samples\\"
    #organism_sequence_files = glob.glob(f"./work/sequences/{organism}*.json")
    organism_sequence_files = glob.glob(f"{work_path}sampled_{organism}*.json")
    for index, organism_sequence_file in enumerate(organism_sequence_files):
        print(organism_sequence_file.replace(work_path, "").split("_"))
        organism_name, sequence_type = organism_sequence_file.replace(work_path, "").split("_")[1:3]
        fname = f'.\\work\\stats\\lstd\\{organism_name}.{sequence_type}.{index}.stats.json'
        if not (os.path.isfile(fname)):
            print('Loading sequences for ' + organism_name)
            sequence_entries = read_json(organism_sequence_file)
            total_entries = len(sequence_entries)
            print(f'Processing {total_entries} entries')
            stats = []
            counter = 0

            for sequence_entry in sequence_entries:
                counter += 1
                print(f'Processing {counter}/{total_entries} for {sequence_type}')
                region = sequence_entry['region']
                sequence = sequence_entry['sequence']
                sequence_stats = get_sequence_stats(region, organism, sequence, [1, 2, 3, 4, 5, 6, 7])
                stats.append(sequence_stats)
            print('Finished processing... saving')
            write_json(fname, stats)
        else:
            print('has processed ' + fname)

def split_data():
    sequence_files = glob.glob("./work/sequences/*.json")
    metadata_reference = []
    main_path = "./work/split_sequences"
    file_index = 0
    for sequence_file in sequence_files:
        file_index += 1
        sequence_entries = read_json(sequence_file)
        chunk_index = 0
        for chunk in chunk_array(sequence_entries, 10000):
            chunk_index += 1
            entry_index = 0
            for entry in chunk:
                entry_index += 1
                seq_id = entry['id']

                folder_path = '/'.join([
                    main_path,
                    str(file_index),
                    str(chunk_index),
                ])

                mkdir(folder_path)

                fasta_path = folder_path + '/' + f'{str(entry_index)}.{seq_id}.fasta'
                #print(folder_path)
                #print(main_path, str(file_index), str(chunk_index))
                #print(fasta_path)
                print(f'Creating ' + fasta_path)
                write_fasta(fasta_path, [{
                    'id': seq_id,
                    'sequence': entry['sequence']
                }])
                del entry['sequence']
                entry['path'] = fasta_path
                metadata_reference.append(entry)

    write_json(f'split_sequence_metadata.json', metadata_reference)
