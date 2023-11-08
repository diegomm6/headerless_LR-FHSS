import numpy as np
import random
from LR_FHSS_DriverMethod import LR_FHSS_DriverFamily
from LiFanMethod import LiFanFamily


#####################
# generate FHS family
#####################
def generate_FHSfamily(params):

    num_frg = params['num_frg']
    num_hdr = params['num_hdr']

    driverFHSfam = LR_FHSS_DriverFamily(q=num_hdr+num_frg, regionDR="EU137")

    return driverFHSfam.FHSfam


#######################
# Generate traffic file
#######################
def generate_traffic(params):

    num_frg = params['num_frg']
    num_trx = params['num_trx']
    num_seq = params['num_seq']
    num_hdr = params['num_hdr']
    t_slots = params['t_slots']
    gran = params['granularity']

    seq_len = num_hdr * round(gran * 233 / 102.4) + gran * num_frg

    T = [] 
    # number of transmissions
    for tx in range(num_trx):
        t = random.randrange(0, t_slots - seq_len, 1)
        s = random.randrange(0, num_seq, 1) 
        T.append((t, s, num_frg))

    return T


###########################
# Generate collision matrix
###########################
def generate_m(seqs, Tt, params):

    num_chn = params['num_chn']
    num_frg = params['num_frg']
    num_hdr = params['num_hdr']
    t_slots = params['t_slots']
    gran = params['granularity']

    collision_matrix = np.ones((t_slots, num_chn)) * -1
    
    # populate m (0=emtpy, 1=rx/coll)
    for tx in Tt:

        time = tx[0]
        sequence = seqs[tx[1]]

        for fh, chn in enumerate(sequence):

            used_slots = gran # write fragment
            if fh < num_hdr:  # write header
                used_slots = round(gran * 233 / 102.4)
                
            for g in range(used_slots):
                collision_matrix[time + g][chn] = 1

            time += used_slots

    return collision_matrix


##########################
# Generate extended matrix
##########################
def generate_extended_m(params):

    num_frg = params['num_frg']
    num_trx = params['num_trx']
    num_hdr = params['num_hdr']
    num_chn = params['num_chn']
    t_slots = params['t_slots']
    gran = params['granularity']

    length = num_hdr + num_frg
    driverFHSfam = LR_FHSS_DriverFamily(q=length, regionDR="EU137")

    extendedFamily = []
    for fhs in driverFHSfam.FHSfam:
        for i in range(10, length+1, 1):
            extendedFamily.append(fhs[:i])

    num_seq = len(extendedFamily)
    max_seq_len = num_hdr * round(gran * 233 / 102.4) + gran * num_frg

    T = [] 
    # number of transmissions
    for tx in range(num_trx):
        t = random.randrange(0, t_slots - max_seq_len, 1)
        s = random.randrange(0, num_seq, 1) 
        T.append((t, s, len(extendedFamily[s])))


    collision_matrix = np.ones((t_slots, num_chn)) * -1
    
    # populate m (0=emtpy, 1=rx/coll)
    for tx in T:

        time = tx[0]
        sequence = extendedFamily[tx[1]]

        for fh, chn in enumerate(sequence):

            used_slots = gran # write fragment
            if fh < num_hdr:  # write header
                used_slots = round(gran * 233 / 102.4)
                
            for g in range(used_slots):
                collision_matrix[time + g][chn] = 1

            time += used_slots

    return extendedFamily, T, collision_matrix


##########################
# Generate extended matrix
##########################
def generate_extended_lifan_m(params):

    num_frg = params['num_frg']
    num_trx = params['num_trx']
    num_hdr = params['num_hdr']
    num_chn = params['num_chn']
    t_slots = params['t_slots']
    gran = params['granularity']

    length = num_hdr + num_frg
    lifanFHSfam = LiFanFamily(q=length, maxfreq=280, mingap=8)

    extendedFamily = []

    max_seq_len = num_hdr * round(gran * 233 / 102.4) + gran * num_frg

    T = []
    # number of transmissions
    for tx in range(num_trx):
        t = random.randrange(0, t_slots - max_seq_len, 1)
        s = random.randrange(0, 540, 1)
        l = random.randint(10, length)
        extendedFamily.append(lifanFHSfam.get_subseq(280, s, l))
        T.append((t, s, l))

    collision_matrix = np.ones((t_slots, num_chn)) * -1
    
    # populate m (0=emtpy, 1=rx/coll)
    for tx in T:

        time = tx[0]
        sequence = extendedFamily[tx[1]]

        for fh, chn in enumerate(sequence):

            used_slots = gran # write fragment
            if fh < num_hdr:  # write header
                used_slots = round(gran * 233 / 102.4)
                
            for g in range(used_slots):
                collision_matrix[time + g][chn] = 1

            time += used_slots

    return extendedFamily, T, collision_matrix
