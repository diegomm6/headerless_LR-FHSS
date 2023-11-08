import numpy as np
from matplotlib import colors
from matplotlib.patches import Patch
from matplotlib.patches import Circle
import matplotlib.pyplot as plt


#######################
# Plot single sequence
#######################
def plot_seq(seq, s, params, file_name):

    num_chn = params['num_chn']
    num_frg = params['num_frg'] 

    seqs_cmap = plt.cm.Set2(np.linspace(0,1,params['num_seq']))
    color = seqs_cmap[s]
    
    # prepare sequence in m 
    m = []  
    for t in range(num_frg):
        m.append([])
        for c in range(num_chn):
            if seq[t] == c:
                m[-1].append(0)
            else:
                m[-1].append(-1)
    m = np.matrix(m)
    # print(m)

    # create font and plot 
    font = {'family': 'serif',
            'weight': 'normal',
            'size': 8}
    plt.rc('font', **font)
    _, axs = plt.subplots(figsize=(6, 6))  

    # plot
    axs.set_title("Sequence: {} {}".format(file_name, seq))
    cmap = colors.ListedColormap(['#00000000', color])

    # transponse to have time in x axis
    axs.imshow(m.transpose(), cmap=cmap, zorder=1.0)

    # limis
    axs.set_xlim(-0.5, num_frg - 0.5)
    axs.set_ylim(-0.5, num_chn - 0.5)

    # set sticks
    time_labels = ['t' + str(x) for x in range(num_frg)]
    axs.set_xticks(np.arange(num_frg), labels=time_labels)
    channel_labels = ['c' + str(x) for x in range(num_chn)]
    axs.set_yticks(np.arange(num_chn), labels=channel_labels)

    # text annotations
    for t in range(num_frg):
        for c in range(num_chn):
            if m[t, c] > -1:
                axs.text(t, c, 'c' + str(c), ha="center", va="center", color="black")

    # highlight starting slot
    circle = Circle((0, seq[0]), radius=0.5, ec='white', fill=False, ls=':', lw=3)
    axs.add_patch(circle)

    # grids
    for c in range(num_chn + 1):
        axs.hlines(c - 0.5, -0.5, num_frg, color='lightgray', ls='-', lw=1, zorder=-100.0)
    for t in range(num_frg + 1):
        axs.vlines(t - 0.5, -0.5, num_chn, color='lightgray', ls='-', lw=1, zorder=-100.0)

    # export
    plt.tight_layout()
    plt.savefig("{}.png".format(file_name), format='png', dpi=300)
    plt.savefig("{}.pdf".format(file_name), format='pdf')
    # plt.savefig("{}.svg".format(file_name), format='svg')

#######################
# Plot T 
#######################
def plot_T(seqs, T, params, file_name):

    num_chn = params['num_chn']
    t_slots = params['t_slots']
    seqs_cmap = plt.cm.Set2(np.linspace(0,1,params['num_seq']))

    # prepare sequence in m 
    m = []  
    for t in range(t_slots):
        m.append([])
        for c in range(num_chn):
            m[-1].append(-1)
    m = np.matrix(m)
    
    for tx in T:
        t = tx[0]
        s = tx[1]
        for ts, p in enumerate(range(len(seqs[s]))): 
            if m[t + ts, seqs[s][p]] > -1:
                # colission
                m[t + ts, seqs[s][p]] = len(seqs_cmap) + 1
            else:
                m[t + ts, seqs[s][p]] = s

    # create font and plot 
    font = {'family': 'serif',
            'weight': 'normal',
            'size': 8}
    plt.rc('font', **font)
    _, axs = plt.subplots(figsize=(6 * t_slots/10, 6))  

    # plot
    axs.set_title("T")

    cmap = colors.ListedColormap(['#00000000'] + [s for s in seqs_cmap] + ['red'])

    # transponse to have time in x axis
    axs.imshow(m.transpose(), cmap=cmap, zorder=1.0)

    # limis
    axs.set_xlim(-0.5, t_slots - 0.5)
    axs.set_ylim(-0.5, num_chn - 0.5)

    # set sticks
    time_labels = ['t' + str(x) for x in range(t_slots)]
    axs.set_xticks(np.arange(t_slots), labels=time_labels)
    channel_labels = ['c' + str(x) for x in range(num_chn)]
    axs.set_yticks(np.arange(num_chn), labels=channel_labels)

    # text annotations
    for t in range(t_slots):
        for c in range(num_chn):
            if m[t, c] == len(seqs_cmap) + 1:
                # collision
                axs.text(t, c, 'collision\nc{}'.format(str(c)), ha="center", va="center", color="black")
            elif m[t, c] > -1:
                # single tx
                axs.text(t, c, 's{:02d}\nc{}'.format(m[t, c],str(c)), ha="center", va="center", color="black")

    # highlight starting slot
    for tx in T:
        t = tx[0]; s = tx[1]
        circle = Circle((t, seqs[s][0]), radius=0.5, ec='white', fill=False, ls=':', lw=3)
        axs.add_patch(circle)

    # grids
    for c in range(num_chn + 1):
        axs.hlines(c - 0.5, -0.5, t_slots, color='lightgray', ls='-', lw=1, zorder=-100.0)
    for t in range(t_slots + 1):
        axs.vlines(t - 0.5, -0.5, num_chn, color='lightgray', ls='-', lw=1, zorder=-100.0)

    # legend
    legends = []
    for s, seq_color in enumerate(seqs_cmap):
        legends.append(Patch(facecolor=seq_color, label='Seq-{:02d}'.format(s)))
    legends.append(Patch(facecolor='red', label='Collision'))
    axs.legend(handles=legends)

    # export
    plt.tight_layout()
    plt.savefig("{}.png".format(file_name), format='png', dpi=300)
    #plt.savefig("{}.svg".format(file_name), format='svg')
    plt.savefig("{}.pdf".format(file_name), format='pdf')


#######################
# Plot M 
#######################
def plot_M(m, params, file_name):

    num_chn = params['num_chn']
    t_slots = params['t_slots']

    # create font and plot 
    font = {'family': 'serif',
            'weight': 'normal',
            'size': 8}
    plt.rc('font', **font)
    _, axs = plt.subplots(figsize=(6 * t_slots/10, 6))  

    # plot
    axs.set_title("M")

    cmap = colors.ListedColormap(['#00000000'] + ['orange'])

    # transponse to have time in x axis
    axs.imshow(m.transpose(), cmap=cmap, zorder=1.0)

    # limis
    axs.set_xlim(-0.5, t_slots - 0.5)
    axs.set_ylim(-0.5, num_chn - 0.5)

    # set sticks
    time_labels = ['t' + str(x) for x in range(t_slots)]
    axs.set_xticks(np.arange(t_slots), labels=time_labels)
    channel_labels = ['c' + str(x) for x in range(num_chn)]
    axs.set_yticks(np.arange(num_chn), labels=channel_labels)

    # text annotations
    for t in range(t_slots):
        for c in range(num_chn):
            if m[t, c] > -1:
                axs.text(t, c, 'c' + str(c), ha="center", va="center", color="black")
    
    # grids
    for c in range(num_chn + 1):
        axs.hlines(c - 0.5, -0.5, t_slots, color='lightgray', ls='-', lw=1, zorder=-100.0)
    for t in range(t_slots + 1):
        axs.vlines(t - 0.5, -0.5, num_chn, color='lightgray', ls='-', lw=1, zorder=-100.0)

    # legend
    legends = []
    legends.append(Patch(facecolor='orange', label='Transmission'))
    legends.append(Patch(facecolor='red', label='Collision'))
    axs.legend(handles=legends)

    # export
    plt.tight_layout()
    plt.savefig("{}.png".format(file_name), format='png', dpi=300)
    #plt.savefig("{}.svg".format(file_name), format='svg')
    plt.savefig("{}.pdf".format(file_name), format='pdf')


#######################
# Print sequences
#######################   
def print_seq(seqs, file_name):

    with open(file_name, 'w') as file:
        for seq in seqs:
            for s in seq:
                file.write('{},'.format(s))
            file.write('\n')

            
#######################
# Print T
#######################   
def print_T(T, file_name):

    with open(file_name, 'w') as file:
        for t in T:
            file.write('{},{},{}\n'.format(t[0], t[1], t[2]))


#######################
# Print metrics
#######################
def print_metrics(Tt, Tp, solve_time, params, file_name = ""):

    tp = 0 # (s, t, p) in T     and  in T'
    fp = 0 # (s, t, p) not in T but  in T'
    fn = 0 # (s, t, p) in T     but  not in T'
    
    for t in Tt:
        if t in Tp:
            tp += 1
            # print('TP:', t)
        else:
            fn += 1
            #print('FN:', t[0], t[1]//23, t[2])
    for t in Tp:
        if t not in Tt:
            fp += 1
            #print('FP:', t[0], t[1]//23, t[2])

    Tt_set = set(Tt)
    Tp_set = set(Tp)

    header = 'chn,slots,seqs,frgs,trxs,seed,'
    string = '{},{},{},{},{},{},'.format(params['num_chn'], params['t_slots'], params['num_seq'], params['num_frg'], params['num_trx'], params['rd_seed'])
    header += 'TP,FP,FN,len(T),len(T\'),dup(T),dup(T\'),time[s]\n'
    string += '{},{},{},{},{},{},{},{:.2f}'.format(tp, fp, fn, len(Tt), len(Tp), len(Tt) - len(Tt_set), len(Tp) - len(Tp_set), solve_time)
    
    if file_name != "":
        with open(file_name, 'w') as file:
            file.write(header)
            file.write(string)
 
    print(string)
    return string

#######################
# Get matrix occupation rate
#######################
def get_matrix_occupation(m, params):

    num_chn = params['num_chn']
    t_slots = params['t_slots']

    slots_occp = 0
    slots_free = 0

    # (0=emtpy, 1=rx, 2=coll)
    for t in range(t_slots):
        for c in range(num_chn):
            if m[t, c] > -1:
                slots_occp += 1
            else:
                slots_free += 1

    return slots_occp/(slots_occp + slots_free)

#######################
# Get matrix collicion rate
#######################
def get_matrix_collision(mc, params):

    num_chn = params['num_chn']
    t_slots = params['t_slots']

    slots_coll = 0
    slots_free_or_occ = 0

    # (0=emtpy, 1=rx, 2=coll)
    for t in range(t_slots):
        for c in range(num_chn):
            if mc[t, c] == 2: 
                slots_coll += 1
            else:
                slots_free_or_occ += 1

    return slots_coll/(slots_coll + slots_free_or_occ)

#######################
# Compute frame decode rate
#######################
def get_decode_rate(mc, seqs, Tt, code_rate, params):

    num_chn = params['num_chn']
    t_slots = params['t_slots']
    num_frg = params['num_frg']

    tx_decoded = 0
    tx_collided = 0
    
    # explore mc (0=emtpy, 1=rx, 2=coll)
    for tx in Tt:
        collisions = 0
        t = tx[0]
        s = tx[1]
        for ts, p in enumerate(range(num_frg)): 
            if mc[t + ts, seqs[s][p]] == 0:
                print("error! matrix cannot be 0 here!")
                quit()
            if mc[t + ts, seqs[s][p]] == 2:
                collisions += 1
        
        if collisions <= num_frg * (1 - code_rate):
            tx_decoded += 1
        else:
            tx_collided += 1

    return tx_decoded/(tx_decoded + tx_collided)