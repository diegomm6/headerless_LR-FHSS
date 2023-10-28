from hdlib_plot import *
from hdlib_generate import *
from hdlib_solve import *
from hdlib_load import *
import random
import time

################
# Main
################

params = {
        'num_chn': 35,      # number of channels [0,1,2,3...,num_chn-1]
        't_slots': 1000,    # number of time slots
        'num_seq': 512,     # number of sequences
        'num_frg': 10,      # number of fragments / sequence lenght
        'num_trx': 500,     # number of transmissions
        'rd_seed': 0,       # random seed
        }

# Generate sequences
seqs = generate_Seq(params)
file_name = 'headerless-06-Seqs.csv'
print_seq(seqs, file_name)

# Generate transmissions
Tt = generate_T(params)
file_name = 'headerless-06-Ttrue.csv'
print_T(Tt, file_name)

# Generate M matrix
m = generate_m(seqs, Tt, params)

