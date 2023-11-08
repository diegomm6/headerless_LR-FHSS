from hdlib_plot import *
#from hdlib_generate import *
#from hdlib_solve import *
#from hdlib_load import *
import random
import time

import seaborn as sns
import matplotlib.pyplot as plt

from generate import *
from solve import *
from solve2 import *

################
# Main
################

num_runs = 1
print('num_runs: {}'.format(num_runs))

num_trxs = [500] #range(100, 3100, 300)
print('num_trxs: {}'.format(num_trxs))

# num_frgs = [10, 20, 40, 80]
num_frgs = [30]
print('num_frgs: {}'.format(num_frgs))


with open('metrics-03-frgs.csv', 'a+') as file:
    header = 'chn,slots,seqs,frgs,trxs,seed,TP,FP,FN,len(T),len(T\'),dup(T),dup(T\'),time[s]'
    file.write(header + '\n')
    print(header)
    for num_trx in num_trxs:
        for num_frg in num_frgs:
            for run in range(num_runs):
                params = {
                        'num_chn': 280,      # number of channels [0,1,2,3...,num_chn-1]
                        't_slots': 7000,    # number of time slots
                        'num_seq': 384,     # number of sequences
                        'num_frg': num_frg, # number of fragments / sequence lenght
                        'num_trx': num_trx, # number of transmissions
                        'rd_seed': run*2,   # random seed
                        'num_hdr': 2,
                        'granularity': 6
                        }

                random.seed(params['rd_seed'])

                # Generate sequences
                #seqs = generate_FHSfamily(params)

                # Generate transmissions
                #Tt = generate_traffic(params)

                # Generate M matrix
                #m = generate_m(seqs, Tt, params)

                seqs, Tt, m = generate_extended_m(params)
                params['num_seq'] = len(seqs)

                
                def find_repeated_elements_in_2d_array(X):
                    k=[]
                    repeated_elements = set()
                    X_tuples = [tuple(subarray) for subarray in X]

                    for i in range(len(X_tuples)):
                        for j in range(i + 1, len(X_tuples)):
                            if X_tuples[i][0] == X_tuples[j][0]:
                                repeated_elements.add(X_tuples[i])
                                repeated_elements.add(X_tuples[j])

                    for element in repeated_elements:
                        k.append(list(element))
                        #print(list(element))  # Convert back to list for printing

                    return k

                k = find_repeated_elements_in_2d_array(Tt)

                t = [tx[0] for tx in k]
                Z = [x for _,x in sorted(zip(t,k))]

                #print(len(Z))
                t = [tx[0] for tx in Tt]
                Z = [x for _,x in sorted(zip(t,Tt))]
                #print("Tt")
                #for i, tx in enumerate(Z):
                #    print(f"{i} t= {tx[0]}, s= {tx[1]//23}, l= {tx[2]}") 

                """
                print(m.shape)
                plt.figure(figsize=(18,12))
                sns.heatmap(m.T)
                plt.title('matrix of observed transmissions, equal to 1 if tx>=1 else -1')
                plt.xlabel('time slots')
                plt.ylabel('channels')
                plt.show()
                #plt.savefig('transmissions-driver.png')
                plt.close('all')
                
                """

                #l=[1]
                #l[1]=1

                # Solve by milp
                start_time = time.process_time()
                Tp = solve_by_milp2(seqs, m, params)
                solve_time = time.process_time() - start_time

                print(f"Tp = {len(Tp)} tx")
                Mp = generate_m(seqs, Tp, params)

                idMp = np.argwhere(Mp > 0)
                idm = np.argwhere(m > 0)

                #print(f"Mp = {len(idMp)}, {m.shape}")
                #print(f"m = {len(idm)}, {Mp.shape}")
                #print(np.sum(m == Mp))
                #print((Mp - m).any())

                #print("\n\nTp")
                #for i, tx in enumerate(Tp):
                #    print(f"{i} t= {tx[0]}, s= {tx[1]//23}, l= {tx[2]}") 


                # Plot metrics
                # file_name = 'trx{:04d}-run{:04d}-Metrics.csv'.format(num_trx, run)
                string = print_metrics(Tt, Tp, solve_time, params)
                file.write(string + ',milp\n')
                file.flush()

