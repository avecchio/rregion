from src.generate import gen_options
from src.extract import extract_sequences
from src.preproc import get_stats, get_region_counts, split_data
from src.proc import down_sample_regions
from src.eda import load_stats, plot_gc_content, plot_seq_len, conduct_pca, bin_numbers, calc_all_jaccard_distances
from src.nmf import nmf_regions
from src.stats import sample_regions
import sys
import yaml
from multiprocessing import Pool
from src.utils import mkdir

def f(x):
    return x*x

cmd = sys.argv[1]

if cmd == 'gen':
    gen_options("./configs", [1e-5, 1e-4, 1e-3])
elif cmd == 'plot':
    pass
elif cmd == 'extract':

    with open("datamap.yml", "r") as stream:
        data_dict = yaml.safe_load(stream)
    mkdir('./work/sequences')
    extract_sequences(data_dict, 'mouse')
    extract_sequences(data_dict, 'human')

elif cmd == 'getstats':
    mkdir('./work/stats')
    get_stats('mouse')
    get_stats('human')
elif cmd == 'sample':
    sample_regions(50000)
elif cmd == 'counts':
    get_region_counts()
elif cmd == 'nmf':
    nmf_regions()

elif cmd == 'sample':
    down_sample_regions(20000)

elif cmd == 'jaccard':
    #for k in [4]:#, 7, 8, 9]:
    calc_all_jaccard_distances()

elif cmd == 'binnum':
    def bf(config):
        bin_numbers(config[0], config[1])

    paths = [
        #("./work/jaccard/enhancers.k5.dist.txt", 'enhancers_dashing_vals.json'),
        #("./work/jaccard/UTR3.dist.k4.txt", 'utr3_k4_dashing_vals.json'),
        #("./work/jaccard/UTR5.dist.k4.txt", 'utr5_k4_dashing_vals.json'),
        #("./work/jaccard/insulators.k4.dist.txt", 'insulators_k4_dashing_vals.json'),
        #("./work/jaccard/ncrna.k4.dist.txt", 'ncrna_k4_dashing_vals.json'),
        ("./work/jaccard/promoters.k3.dist.txt", 'promoters_k3_dashing_vals.json'),
        #("./work/jaccard/silencers.k4.dist.txt", 'silencers_k4_dashing_vals.json'),
        #("./work/jaccard/silencers.k5.dist.txt", 'silencers_dashing_vals.json'),
        #("./work/jaccard/silencers.k3.dist.txt", 'silencers_k3_dashing_vals.json'),
        #("./work/jaccard/UTR3.k5.dist.txt", 'UTR3_dashing_vals.json'),
        #("./work/jaccard/UTR5.k5.dist.txt", 'UTR5_dashing_vals.json'),
        #("./work/jaccard/ncrna.k5.dist.txt", 'ncrna_dashing_vals.json')
    ]

    with Pool(6) as p:
        print(p.map(bf, paths))


    #bin_numbers("./work/jaccard/silencers.k5.dist.txt", 'silencers_dashing_vals.json')
    
elif cmd == 'eda':
    stats_data = load_stats()
    #plot_gc_content(stats_data)
    #plot_seq_len(stats_data)

elif cmd == 'pca':
    conduct_pca()

elif cmd == 'split_sequences':
    split_data()