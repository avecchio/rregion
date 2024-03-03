from src.utils import mkdir, read_json, write_fasta, write_json
from src.plot import boxen, scatterplot
from src.stats import generate_combinations
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

def conduct_pca(n_components, features, labels, plot_name, title, scales = []):
    pca = PCA(n_components=n_components)
    principal_components = pca.fit_transform(features)

    pc_df = pd.DataFrame(data=principal_components, columns=[f'PC{i}' for i in range(1, n_components + 1)])
    pc_df['Region'] = labels
    scatterplot(pc_df, plot_name, {
        'xdata': 'PC1',
        'ydata': 'PC2',
        'xlabel': 'PC1',
        'ylabel': 'PC2',
        'title': title,
        'xscalelog': ('xscalelog' in scales)
        'yscalelog': ('yscalelog' in scales)
    }, datalabels='Region')


def distribution_analysis():
    #stats_files = glob.glob("sampled_human_*.json")
    kmer_lengths = [1, 2, 3, 4, 5, 6, 7]
    kmer_combos = generate_combinations(kmer_lengths, 1, 3)

    '''
    data = []
    for stat_file in stats_files:
        elements = read_json(stat_file)
        subset = random.sample(elements, 200)
        data += subset
    write_json(f'test_subsample_200.json', data)
    '''
    data = read_json(f'test_subsample_200.json')
    for kmer_combo in kmer_combos:
        combo_name = '.'.join([str(k) for k in list(kmer_combo)])
        print(combo_name)
        parsed_data = []
        for element in data:
            restructured_dict = {
                'organism': element['organism'],
                'region': element['region_name']
            }
            for kmer_length in kmer_combo:
                ks = element['kmers'][str(kmer_length)]
                restructured_dict.update(ks)
            parsed_data.append(restructured_dict)
        df = pd.DataFrame(parsed_data)
        df = df.fillna(0)
        # remove unused column
        df = df.drop('organism', axis=1)
        # return all features but the labels
        features = df.drop('region', axis=1)

def pca_analysis():
    stats_files = glob.glob("sampled_human_*.json")
    kmer_lengths = [1, 2, 3, 4, 5, 6, 7]
    kmer_combos = generate_combinations(kmer_lengths, 1, 3)

    data = []
    for stat_file in stats_files:
        elements = read_json(stat_file)
        #subset = random.sample(elements, 200)
        data += elements
        #data += subset
    #write_json(f'test_subsample_200.json', data)
    
    #data = read_json(f'test_subsample_200.json')
    for kmer_combo in kmer_combos:
        combo_name = '.'.join([str(k) for k in list(kmer_combo)])
        print(combo_name)
        parsed_data = []
        for element in data:
            region_name = element['region_name']
            if region_name not in ['UTR5', 'UTR3']:
                restructured_dict = {
                    'organism': element['organism'],
                    'region': region_name
                }
                for kmer_length in kmer_combo:
                    ks = element['kmers'][str(kmer_length)]
                    restructured_dict.update(ks)
                parsed_data.append(restructured_dict)
        df = pd.DataFrame(parsed_data)
        df = df.fillna(0)        # remove unused column
        df = df.drop('organism', axis=1)
        # return all features but the labels
        features = df.drop('region', axis=1)
        labels = df['region']
        for num_components in [2, 4, 6, 8, 10]:
            try:
                title_combo_name = ','.join([str(k) for k in list(kmer_combo)])
                plot_title = f'PCA of Kmer Distributions (k={title_combo_name},ncomp={num_components})'

                plot_name = f'./work/plots/pca_n{num_components}_{combo_name}_scatter.png'
                conduct_pca(num_components, features, labels, plot_name, plot_title)

                scaler = StandardScaler()
                standardized_features = scaler.fit_transform(features)

                standardized_plot_title = f'PCA of Standardized Kmer Distributions (k={title_combo_name},ncomp={num_components})'
                standardized_plot_name = f'./work/plots/pca_standardized_n{num_components}_{combo_name}_scatter.png'
                conduct_pca(num_components, standardized_features, labels, standardized_plot_name, standardized_plot_title)

                log_plot_title = f'Log10 PCA of Kmer Distributions (k={title_combo_name},ncomp={num_components})'
                log_plot_name = f'./work/plots/pca_xscalelog_n{num_components}_{combo_name}_scatter.png'
                conduct_pca(num_components, standardized_features, labels, log_plot_name, log_plot_title, ['xscalelog'])

                log_plot_title = f'Log10 PCA of Kmer Distributions (k={title_combo_name},ncomp={num_components})'
                log_plot_name = f'./work/plots/pca_logscale_n{num_components}_{combo_name}_scatter.png'
                conduct_pca(num_components, standardized_features, labels, log_plot_name, log_plot_title, ['xscalelog', 'yscalelog'])
                return
            except:
                pass