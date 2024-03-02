from src.utils import mkdir, read_json, write_fasta, write_json
from src.plot import boxen
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pandas as pd
import glob
import random
#import seaborn as sns
#import matplotlib.pyplot as plt

def name_filter(name, items):
    for item in items:
        if item in name:
            return True
    return False

def load_stats():
    stats_files = glob.glob("./work/stats/*.json")
    total = len(stats_files)
    simple_stats_data = []
    counter = 0
    sampled_stats_data = []
    for stats_file in stats_files:
        if name_filter(stats_file, ['enhancers', 'insulators', 'silencers', 'promoters']):
            counter += 1
            print(f'Loading {counter} of {total}')
            stats_data = read_json(stats_file)
            extracted_stats_data = []
            for entry in stats_data:
                del entry['kmers']
                extracted_stats_data.append(entry)
            simple_stats_data += extracted_stats_data
            sampled_stats_data += random.sample(simple_stats_data, round(len(extracted_stats_data) / 1000))
    write_json(f'stats_comb.samp.json', sampled_stats_data)
    write_json(f'stats_comb.json', simple_stats_data)
    return pd.DataFrame(read_json(f'stats_comb.samp.json'))

def plot_gc_content(stats_data):

    print(f'Plotting: ' + 'GC Content Per Region')
    boxen(stats_data, f'./work/plots/gc_content_region.png', {
        'xdata': 'gc_content',
        'hdata': 'region_name',
        'xlabel': 'GC Content (%)',
        'xlogscale': False,
        'ylogscale': False,
        'ylabel': 'Frequency',
        'title': 'GC Content Per Region',
    })

    print(f'Plotting: ' + 'GC Content Per Organism')
    boxen(stats_data, f'./work/plots/gc_content_organism.png', {
        'xdata': 'gc_content',
        'hdata': 'organism',
        'xlabel': 'GC Content (%)',
        'xlogscale': False,
        'ylogscale': False,
        'ylabel': 'Frequency',
        'title': 'GC Content Per Organism',
    })

def plot_seq_len(stats_data):

    print(f'Plotting: ' + 'Sequence Length Per Region')
    boxen(stats_data, f'./work/plots/seq_len_region.png', {
        'xdata': 'len',
        'hdata': 'region_name',
        'xlabel': 'Sequence Len (bp)',
        'xlogscale': True,
        'ylogscale': False,
        'ylabel': 'Frequency',
        'title': 'Sequence Length Per Region',
    })

    print(f'Plotting: ' + 'Sequence Length Per Organism')
    boxen(stats_data, f'./work/plots/seq_len_organism.png', {
        'xdata': 'len',
        'hdata': 'organism',
        'xlabel': 'Sequence Len (bp)',
        'xlogscale': True,
        'ylogscale': False,
        'ylabel': 'Frequency',
        'title': 'Sequence Length Per Organism',
    })

def bin_numbers(data_file, save_file):
    bin_counts = {}
    with open(data_file) as dashf:
        counter = 0
        for line in dashf:
            counter += 1
            print(data_file, counter)
            if counter > 1:
                for num_str in line.strip().split("\t"):
                    try:
                        num = str(round(float(num_str), 2))
                        if num not in bin_counts:
                            bin_counts[num] = 0
                        bin_counts[num] += 1
                    except:
                        pass
    write_json(save_file, bin_counts)

    #bin_num = str(round(num, 2))
    
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

def conduct_pca():
    stats_files = glob.glob("sampled_human_*.json")

    total = len(stats_files)
    for k in [1, 2, 3, 4]:
        counter = 0
        kstats = []
        for stats_file in stats_files:
            counter += 1
            print(f'Loading {counter} of {total}')
            stats_data = read_json(stats_file)
            entry_kstats = []
            for entry in stats_data:
                kdata = entry['kmers'][str(k)]
                #print(list(entry.keys()))
                kdata['region'] = entry['region_name']
                kdata['organism'] = entry['organism']
                entry_kstats.append(kdata)
            kstats += random.sample(entry_kstats, round(len(entry_kstats) / 1000))
        write_json(f'k{k}_stats.samp.json', kstats)
        
    # Example dataset with DNA sequences
    #data = {
    #    'Seq1': 'ATCGATCGATCG',
    #    'Seq2': 'CGATCGATCGAT',
    #    'Seq3': 'GATCGATCGATC',
    #    # ... add more sequences as needed
    #}

    #df = pd.DataFrame(data)

    #from collections import Counter

    #def generate_kmer_counts(sequence, k):
    #    kmers = [sequence[i:i+k] for i in range(len(sequence) - k + 1)]
    #    return dict(Counter(kmers))

    #k = 3  # Specify the desired k-mer length
    #df_kmer_counts = df.applymap(lambda seq: generate_kmer_counts(seq, k)).fillna(0)
    '''
    scaler = StandardScaler()
    standardized_data = scaler.fit_transform(df_kmer_counts)

    n_components = 2
    pca = PCA(n_components=n_components)

    principal_components = pca.fit_transform(standardized_data)

    pc_df = pd.DataFrame(data=principal_components, columns=[f'PC{i}' for i in range(1, n_components + 1)])
    '''
    # Scatter plot
    #sns.scatterplot(x='PC1', y='PC2', data=pc_df)
    #plt.xlabel('Principal Component 1')
    #plt.ylabel('Principal Component 2')
    #plt.title('PCA of k-mers')
    #plt.show()