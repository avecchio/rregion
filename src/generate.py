from itertools import product
from src.utils import write_json

def gen_options(path, options_range):
    opt_options = ['Adamax', 'Adagrad', 'Adam', 'SGD']
    l2_options = options_range
    lr_options = options_range
    eps_options = options_range
    all_options = [opt_options, l2_options, lr_options, eps_options]
    combinations = list(product(*all_options))

    counter = 0
    for combination in combinations:
        counter += 1
        configuration = {
            'opt': combination[0],
            'l2': combination[1],
            'lr': combination[2],
            'eps': combination[3]
        }

        write_json(f'{path}/configuration_{counter}.json', configuration)
