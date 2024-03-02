import pandas as pd
from Bio import SeqIO
from src.utils import write_json, read_json, write_file, write_fasta
import glob

def read_fasta_file(fasta_file):
    fasta_data = {}
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            fasta_data[record.id] = str(record.seq)
    return fasta_data

#def jsontofasta(workdir):
#    json_files = glob.glob(f'./work/sequences/*.json')
#    for json_file in json_files:
#        sequence_entries = read_json(json_file)
#        for entry in sequence_entries:
#            region = entry["region"]
#        write_file
#        write_fasta


def load_genome(chromosome_map):
    genome_dict = {}
    for chr in chromosome_map:
        fasta_file = chromosome_map[chr]
        try:
            fasta_dict = read_fasta_file(fasta_file)
            keyid = list(fasta_dict.keys())[0]
            genome_dict[chr] = fasta_dict[keyid]
        except:
            print('could not load ' + chr)
    return genome_dict

def extract_coords(coord_str):
    region_name, start_end = coord_str.split(":")
    start, end = start_end.split("-")
    return region_name, int(start), int(end)

def extract_untranslated_regions(fasta_file):
    print(fasta_file)
    fasta_dict = read_fasta_file(fasta_file)
    utr_sequences = []
    for entry_header in fasta_dict:
        sequence = fasta_dict[entry_header]
        utr3_coord_str = None
        utr5_coord_str = None
        for component in entry_header.split("|"):
            if 'UTR3' in component:
                utr3_coord_str = component
            elif 'UTR5' in component:
                utr5_coord_str = component
        
        if utr3_coord_str is not None:
            utr3_name, utr3_start, utr3_end = extract_coords(utr3_coord_str)
            utr3_seq = sequence[utr3_start - 1: utr3_end + 1]
            utr_sequences.append({'region': utr3_name, 'sequence': utr3_seq})
        
        if utr5_coord_str is not None:
            utr5_name, utr5_start, utr5_end = extract_coords(utr5_coord_str)
            utr5_seq = sequence[utr5_start - 1: utr5_end + 1]
            utr_sequences.append({'region': utr5_name, 'sequence': utr5_seq})

    return utr_sequences

def read_bed_file(bed_file, sep='\t'):
    with open(bed_file, 'r') as bfp:
        lines = bfp.readlines()

    coords = []
    for line in lines:
        line_entry = line.split(sep)
        coords.append({
            'chr': line_entry[0],
            'start': int(line_entry[1]),
            'end': int(line_entry[2])
        })
    return coords

def read_insulator_file(insulator_file, loc_col, organism_filter):
    with open(insulator_file, 'r') as bfp:
        lines = bfp.readlines()

    coords = []
    for line in lines:
        line_entry = line.split("\t")
        organism = line_entry[1]
        if organism == organism_filter:
            loc = line_entry[loc_col]
            #print(line_entry[0], loc)
            if ':' in loc and loc != ":0-0":
                chr, coord = loc.split(":")
                start, end = coord.split("-")

                coords.append({
                    'chr': chr,
                    'start': int(start),
                    'end': int(end)
                })
    return coords

def divide_chunks(l, n):     
    # looping till length l 
    for i in range(0, len(l), n):  
        yield l[i:i + n] 
  
def extract_sequence(genome_dict, region, coords):
    sequences = []
    print('Extracting sequences for region:' + region)
    total = 0
    for coord in coords:
        chr = coord['chr']
        start = coord['start']
        end = coord['end']
        if chr in genome_dict:
            total += (end - start + 1)
            print(chr, start, end, total)
            chr_seq = genome_dict[chr]
            region_seq = chr_seq[start:end + 1]
            sequences.append({'region': region, 'sequence': region_seq})
        else:
            pass
            #print(f'{chr} not in genome dictionary')
    return sequences    

def assign_uuids(coords, reference, counter = 0):
    updated_coords = []
    for coord in coords:
        counter += 1
        region = coord['region']
        coord['id'] = f'{reference}_{region}_{str(counter)}'
        coord['organism'] = reference
        updated_coords.append(coord)
    return counter, updated_coords

def extract_sequences(datamap, ref_name):
    work_dir = "./work/sequences"
    print('Extracting for ' + ref_name)
    data_ref = datamap[ref_name]

    genome_dict = load_genome(data_ref['genome'])

    counter = 0
    utr_sequences = extract_untranslated_regions(data_ref['transcripts'])
    c, utr_seqs = assign_uuids(utr_sequences, ref_name, counter)
    write_json(f'{work_dir}/{ref_name}.utrs.1.sequences.json', utr_seqs)
    counter = c
    
    #print({
    #    'promoters': len(promoter_coords),
    #    'silencers': len(silencer_coords),
    #    'ncrna': len(ncrna_coords),
    #    'insulators': len(insulator_comp_coords) + len(insulator_exp_coords)
    #})

    promoter_counter = 0
    promoter_coords = read_bed_file(data_ref['promoters'], sep=' ')
    print(f'Processing {len(promoter_coords)} promoters')
    for i, coords_set in enumerate(divide_chunks(promoter_coords, 10000)):
        promoter_sequences = extract_sequence(genome_dict, 'promoters', coords_set)
        c, promotoer_seqs = assign_uuids(promoter_sequences, ref_name, promoter_counter)
        promoter_counter = c
        write_json(f'{work_dir}/{ref_name}.promoter.{i}.sequences.json', promotoer_seqs)
    promoter_coords = []

    silencer_counter = 0
    silencer_coords = read_bed_file(data_ref['silencers'])
    print(f'Processing {len(silencer_coords)} silencers')
    for i, coords_set in enumerate(divide_chunks(silencer_coords, 10000)):
        silencer_sequences = extract_sequence(genome_dict, 'silencers', coords_set)
        c, promotoer_seqs = assign_uuids(silencer_sequences, ref_name, silencer_counter)
        silencer_counter = c
        write_json(f'{work_dir}/{ref_name}.silencer.{i}.sequences.json', promotoer_seqs)
    silencer_coords = []

    ncrna_coords = read_bed_file(data_ref['ncrna'])
    print(f'Processing {len(ncrna_coords)} ncrna')
    ncrna_counter = 0
    for i, coords_set in enumerate(divide_chunks(ncrna_coords, 10000)):
        ncrna_sequences = extract_sequence(genome_dict, 'ncrna', coords_set)
        c, ncrna_seqs = assign_uuids(ncrna_sequences, ref_name, ncrna_counter)
        ncrna_counter = c
        write_json(f'{work_dir}/{ref_name}.nrcna.{i}.sequences.json', ncrna_seqs)
    ncrna_coords = []

    insulator_exp_coords = read_insulator_file("./data/insulators/allexp.txt", 3, ref_name.capitalize())
    insulator_comp_coords = read_insulator_file("./data/insulators/allcomp.txt", 2, ref_name.capitalize())
    insulator_sequences = []
    insulator_sequences += extract_sequence(genome_dict, 'insulators', insulator_exp_coords)
    insulator_sequences += extract_sequence(genome_dict, 'insulators', insulator_comp_coords)
    c, insulator_seqs = assign_uuids(insulator_sequences, ref_name)
    write_json(f'{work_dir}/{ref_name}.insulator.sequences.json', insulator_seqs)
    insulator_comp_coords = []
    insulator_exp_coords = []
    insulator_sequences = []

    counter = 0
    enhancer_counter = 0
    for enhancer_file in data_ref['enhancers']:
        counter += 1
        enhancer_coords = read_bed_file(enhancer_file)
        print(f'Processing {len(enhancer_coords)} enhancers')
        enhancer_sequences = extract_sequence(genome_dict, 'enhancers', enhancer_coords)
        c, enhancer_seqs = assign_uuids(enhancer_sequences, ref_name, enhancer_counter)
        enhancer_counter = c
        enhancer_coords = []
        write_json(f'{work_dir}/{ref_name}.enhancers.{counter}.sequences.json', enhancer_seqs)
    