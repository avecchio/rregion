from src.generate import gen_options
from src.extract import extract_sequences
from src.preproc import get_stats, calc_all_jaccard_distances, get_region_counts
from src.proc import down_sample_regions
from src.eda import load_stats, plot_gc_content, plot_seq_len, conduct_pca
import sys
import yaml

cmd = sys.argv[1]

if cmd == 'gen':
    gen_options("./configs", [1e-5, 1e-4, 1e-3])
elif cmd == 'plot':
    pass
elif cmd == 'extract':

    with open("datamap.yml", "r") as stream:
        data_dict = yaml.safe_load(stream)
    extract_sequences(data_dict, 'mouse')
    extract_sequences(data_dict, 'human')

elif cmd == 'getstats':
    get_stats('mouse')
    get_stats('human')

elif cmd == 'counts':
    get_region_counts()

elif cmd == 'sample':
    down_sample_regions(20000)

elif cmd == 'jaccard':
    #for k in [4]:#, 7, 8, 9]:
    calc_all_jaccard_distances()

elif cmd == 'eda':
    stats_data = load_stats()
    plot_gc_content(stats_data)
    plot_seq_len(stats_data)

elif cmd == 'pca':
    conduct_pca()