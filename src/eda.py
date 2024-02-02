from src.utils import read_json, write_json
from src.plot import boxen
from sklearn.preprocessing import StandardScaler
from sklearn.decomposition import PCA
import pandas as pd
import glob
import random
#import seaborn as sns
#import matplotlib.pyplot as plt

def load_stats():
    '''
    stats_files = glob.glob("./work/stats/*.json")
    total = len(stats_files)
    simple_stats_data = []
    counter = 0
    sampled_stats_data = []
    for stats_file in stats_files:
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
    '''
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


def conduct_pca():

    stats_files = glob.glob("./work/stats/*.json")
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