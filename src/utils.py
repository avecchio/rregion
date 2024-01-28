from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import json
import os
import pandas as pd

def mkdir(path):
    os.mkdir(path)

def write_json(fname, data):
    with open(fname, 'w') as jf:
        jf.write(json.dumps(data))

def read_json(fname):
    with open(fname, 'r') as jf:
        data = json.loads(jf.read())
    return data

def write_fasta(filename, sequence_entries):
    records = []
    for sequence_entry in sequence_entries:
        entry_id = sequence_entry['id']
        sequence = Seq(sequence_entry['sequence'])
        records.append(sequence, id=entry_id, description="")
    SeqIO.write(records, filename, "fasta")

def write_file(filename, lines):
    with open(filename, 'w') as fp:
        for line in lines:
            fp.write(f'{line}\n')

def write_csv(filename, data):
    df = pd.DataFrame(data)
    df.to_csv(filename, index=False)