import numpy as np
import random


#######################
# generate sequence file
#######################
def generate_Seq(params):

    num_seq = params['num_seq']
    num_frg = params['num_frg']
    num_chn = params['num_chn']

    seqs = []
    for _sn in range(num_seq):
        seqs.append([])
        for _sl in range(num_frg):
            seqs[-1].append(random.randrange(0, num_chn, 1))

    return seqs

            
#######################
# Generate traffic file
#######################
def generate_T(params):

    num_frg = params['num_frg']
    num_trx = params['num_trx']
    t_slots = params['t_slots']
    num_seq = params['num_seq']

    T = [] 
    # number of transmissions
    for tx in range(num_trx):
        t = random.randrange(0, t_slots - num_frg, 1)
        s = random.randrange(0, num_seq, 1) 
        T.append((t, s, num_frg))

    return T


#######################
# Generate m
#######################
def generate_m(seqs, Tt, params):

    num_chn = params['num_chn']
    t_slots = params['t_slots']
    num_frg = params['num_frg']

    # create matrix m (all -1)
    m = []  
    for t in range(t_slots):
        m.append([])
        for c in range(num_chn):
            m[-1].append(-1)
    m = np.matrix(m)
    
    # populate m (0=emtpy, 1=rx/coll)
    for tx in Tt:
        t = tx[0]
        s = tx[1]
        for ts, p in enumerate(range(num_frg)): 
            if m[t + ts, seqs[s][p]] > -1:
                # colission
                m[t + ts, seqs[s][p]] = 1
            else:
                m[t + ts, seqs[s][p]] = 1

    return m

#######################
# Generate mc (collissions marked)
#######################
def generate_mc(seqs, Tt, params):

    num_chn = params['num_chn']
    t_slots = params['t_slots']
    num_frg = params['num_frg']

    # create matrix m (all -1)
    m = []  
    for t in range(t_slots):
        m.append([])
        for c in range(num_chn):
            m[-1].append(-1)
    m = np.matrix(m)
    
    # populate m (0=emtpy, 1=rx, 2=coll)
    for tx in Tt:
        t = tx[0]
        s = tx[1]
        for ts, p in enumerate(range(num_frg)): 
            if m[t + ts, seqs[s][p]] > -1:
                # colission
                m[t + ts, seqs[s][p]] = 2
            else:
                m[t + ts, seqs[s][p]] = 1

    return m