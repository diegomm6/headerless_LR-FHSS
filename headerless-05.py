from hdlib_plot import *
from hdlib_generate import *
from hdlib_solve import *
from hdlib_load import *
import random
import time

################
# Main
################

num_runs = 10
print('num_runs: {}'.format(num_runs))
num_trxs = range(500, 3400, 100)
print('num_trxs: {}'.format(num_trxs))
# num_frgs = [10, 20, 40, 80]
num_frgs = range(10, 100, 20)
print('num_frgs: {}'.format(num_frgs))

with open('metrics-05-matr-occ.csv', 'a+') as file:
    header = 'chn,slots,seqs,frgs,trxs,seed,m_occ,m_col,tx_dec_rate_1_3,tx_dec_rate_2_3'
    file.write(header + '\n')
    print(header)
    for num_trx in num_trxs:
        for num_frg in num_frgs:
            for run in range(num_runs):
                params = {
                        'num_chn': 35,      # number of channels [0,1,2,3...,num_chn-1]
                        't_slots': 1000,    # number of time slots
                        'num_seq': 512,     # number of sequences
                        'num_frg': num_frg, # number of fragments / sequence lenght
                        'num_trx': num_trx, # number of transmissions
                        'rd_seed': run*2,   # random seed
                        }
                
                random.seed(params['rd_seed'])

                # Generate sequences
                seqs = generate_Seq(params)

                # Generate transmissions
                Tt = generate_T(params)

                # Generate Mc matrix (0=emtpy, 1=rx, 2=coll)
                # m = generate_m(seqs, Tt, params)
                mc = generate_mc(seqs, Tt, params)

                m_occ = get_matrix_occupation(mc, params)
                m_col = get_matrix_collision(mc, params)
                tx_dec_rate_1_3 = get_decode_rate(mc, seqs, Tt, 0.333, params)
                tx_dec_rate_2_3 = get_decode_rate(mc, seqs, Tt, 0.666, params)

                string = "{},{},{},{},{},{},".format(params['num_chn'],params['t_slots'],params['num_seq'],params['num_frg'],params['num_trx'],params['rd_seed'],)
                string += "{},{},{},{}".format(m_occ, m_col, tx_dec_rate_1_3, tx_dec_rate_2_3)
                file.write(string + '\n')

                # assumption:
                # Tp is always equal to Tt







