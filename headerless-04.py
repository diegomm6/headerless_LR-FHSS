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
num_trxs = range(500, 10000, 100)
print('num_trxs: {}'.format(num_trxs))
# num_frgs = [10, 20, 40, 80]
num_frgs = range(10, 100, 20)
print('num_frgs: {}'.format(num_frgs))


with open('metrics-04-frgs.csv', 'a+') as file:
    header = 'chn,slots,seqs,frgs,trxs,seed,TP,FP,FN,len(T),len(T\'),dup(T),dup(T\'),time[s]'
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
                # Print sequences
                # file_name = 'trx{:04d}-run{:04d}-Seqs.csv'.format(num_trx, run)
                # print_seq(seqs, file_name)

                # Plot all sequences
                # file_name = 'run{:04d}-Seq'.format(run)
                # for s, seq in enumerate(seqs):
                #     plot_seq(seq, s, params, file_name + '{:02d}'.format(s))

                # Generate transmissions
                Tt = generate_T(params)
                # Print Tp
                # file_name = 'trx{:04d}-run{:04d}-Ttrue.csv'.format(num_trx, run)
                # print_T(Tt, file_name)
                # Generate M matrix
                m = generate_m(seqs, Tt, params)

                # Plot transmissions
                # file_name = 'run{:04d}-Ttrue'.format(run)
                # plot_T(seqs, Tt, params, file_name)
                # file_name = 'run{:04d}-M'.format(run)

                # # Solve by milp
                # start_time = time.process_time()
                # Tp = solve_by_milp1(seqs, m, params)
                # solve_time = time.process_time() - start_time
                # # Print Tpred
                # # file_name = 'trx{:04d}-run{:04d}-Tpred.csv'.format(num_trx, run)
                # # print_T(Tp, file_name)

                # # Plot metrics
                # # file_name = 'trx{:04d}-run{:04d}-Metrics.csv'.format(num_trx, run)
                # string = print_metrics(Tt, Tp, solve_time, params)
                # file.write(string + ',milp\n')
                # file.flush()

                # Solve by oana's
                start_time = time.process_time()
                Tp = solve_by_oana1(seqs, m, params)
                solve_time = time.process_time() - start_time

                # Plot metrics
                string = print_metrics(Tt, Tp, solve_time, params)
                file.write(string + ',oana\n')
                file.flush()

