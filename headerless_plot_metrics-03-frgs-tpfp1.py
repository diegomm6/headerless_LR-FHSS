import numpy as np
from matplotlib import colors
from matplotlib.patches import Patch
from matplotlib.patches import Circle
import matplotlib.pyplot as plt
import csv
from statistics import mean

num_chn = []
t_slots = []
num_seq = []
num_frg = []
num_trx = []
rd_seed = []
tp = []
fp = []
fn = []
lenTt = []
lenTp = []
dupTt = []
dupTp = []
time = []

file_name = 'metrics'
with open(file_name + '.csv', 'r') as file:
    csv_reader = csv.reader(file, delimiter=',')
    for r, row in enumerate(csv_reader):
        if row[0] == 'chn':
            continue

        num_chn.append(int(row[0]))
        t_slots.append(int(row[1]))
        num_seq.append(int(row[2]))
        num_frg.append(int(row[3]))
        num_trx.append(int(row[4]))
        rd_seed.append(int(row[5]))

        tp.append(int(row[6]))
        fp.append(int(row[7]))
        fn.append(int(row[8]))
        lenTt.append(int(row[9]))
        lenTp.append(int(row[10]))
        dupTt.append(int(row[11]))
        dupTp.append(int(row[12]))
        time.append(float(row[13]))


            
# create font and plot 
font = {'family': 'serif',
        'weight': 'normal',
        'size': 8}
plt.rc('font', **font)
_, axs = plt.subplots(figsize=(6, 6)) 

#axs.set_title("Headerless LR-FHSS MILP model 1")
axs.set_title("Headerless LR-FHSS Detection Effectiveness (TP, FP and FN)") # (ILP model / Heuristic)
axs.set_xlabel('Frame Transmissions')
axs.set_ylabel('TP and FP (log)')
# axs.set_xscale('log')
axs.set_yscale('log')
# axs.grid(True)
# limis
# axs.set_xlim( , )
axs.set_ylim(bottom=1, top=10000000)

axs.grid(axis='x', which='both', color='0.95')
axs.grid(axis='y', which='major', color='0.93')
axs.grid(axis='y', which='minor', color='0.96')

num_frg_set = sorted(list(set(num_frg)))
cmap = plt.cm.Set1(np.linspace(0, 1, len(num_frg_set)))
for f, frgs in enumerate(num_frg_set):

        num_trx_set = sorted(list(set(num_trx)))
        tp_min = []
        tp_avg = []
        tp_max = []
        fp_min = []
        fp_avg = []
        fp_max = []
        fn_min = []
        fn_avg = []
        fn_max = []

        for trx_set in num_trx_set:
            tp_tmp = []
            fp_tmp = []
            fn_tmp = []
            for t, trx in enumerate(num_trx):
                if trx == trx_set:
                    if frgs == num_frg[t]:
                        tp_tmp.append(tp[t])
                        fp_tmp.append(fp[t])
                        fn_tmp.append(fn[t])

            # to be able to plot while still populating the csv
            if len(tp_tmp) == 0:
                tp_min.append(tp_min[-1])
                tp_avg.append(tp_avg[-1])
                tp_max.append(tp_max[-1])
                fp_min.append(fp_min[-1])
                fp_avg.append(fp_avg[-1])
                fp_max.append(fp_max[-1])
                fn_min.append(fn_min[-1])
                fn_avg.append(fn_avg[-1])
                fn_max.append(fn_max[-1])
                continue
            
            tp_min.append(min(tp_tmp))
            tp_avg.append(mean(tp_tmp))
            tp_max.append(max(tp_tmp))
            fp_min.append(min(fp_tmp))
            fp_avg.append(mean(fp_tmp))
            fp_max.append(max(fp_tmp))
            fn_min.append(min(fn_tmp))
            fn_avg.append(mean(fn_tmp))
            fn_max.append(max(fn_tmp))

        max_num_trx = num_chn[0] * t_slots[0] / frgs
        #axs.vlines(max_num_trx, ls='-.', lw=1, color=cmap[f], ymin=0, ymax=10000000)
        #if f == len(num_frg_set) - 1:
        #    axs.text(max_num_trx - 85, 1.2, 'slots = fragments', rotation='vertical', color='gray')

        axs.fill_between(num_trx_set, tp_min, tp_max, alpha=0.3, color=cmap[f], zorder=10)
        axs.fill_between(num_trx_set, fp_min, fp_max, alpha=0.3, color=cmap[f], zorder=10)
        axs.fill_between(num_trx_set, fn_min, fn_max, alpha=0.3, color=cmap[f], zorder=10)

        if f == 0: # only print once
            axs.plot(num_trx_set, num_trx_set, marker='', ls=':', lw=2, label='Tx ()', color='black', zorder=20)
        axs.plot(num_trx_set, tp_avg, marker='P', label='TP ({})'.format(''), color=cmap[f], zorder=15)
        axs.plot(num_trx_set, fp_avg, marker='v', label='FP ({})'.format(''), color=cmap[f], zorder=15)
        axs.plot(num_trx_set, fn_avg, marker='o', label='FN ({})'.format(''), color=cmap[f])

text = "TP for all curves\nare overlapped here\n(TP=?)"
#axs.annotate(text, xy=(num_trx_set[12], num_trx_set[12]), xycoords='data',
#            xytext=(2200, 10), textcoords='data', # textcoords='axes fraction',
#            arrowprops=dict(arrowstyle="->",
#                            connectionstyle="arc3, rad=0.2"),
#            horizontalalignment='center', verticalalignment='top',
#            )

text = ""
text += "Channels: {}\n".format(num_chn[0])
text += "Time slots: {}\n".format(t_slots[0])
text += "Sequences: {}\n".format(num_seq[0])
text += "Fragments: 8 to 30\n"
text += "Headers: 2\n"
text += "Runs: {}".format(int(max(rd_seed)/2 + 1))
# text += "transmissions: {}\n".format(num_trx[0])
axs.text(0.03, 0.98, text, ha='left', va='top', transform=axs.transAxes,
        bbox={'facecolor':'white', 'edgecolor':'lightgray', 'alpha': 0.8, 'pad': 3})

text = ""
text += "Shadowed areas indicates min/max spread\n"
text += "True Positives (TP) = ??\n"
text += "False Negatives (FN) = ??"
axs.text(0.98, 0.98, text, ha='right', va='top', transform=axs.transAxes,
        bbox={'facecolor':'white', 'edgecolor':'lightgray', 'alpha': 0.8, 'pad': 3})

# print(axs.get_xlim())

axs.legend(loc='lower right')
# plt.tight_layout()
plt.savefig("{}-tpfp1.png".format(file_name), format='png', dpi=300)
plt.savefig("{}-tpfp1.pdf".format(file_name), format='pdf')

print("xlim: {}".format(axs.get_xlim()))





# # boxplot data
# box_data = []
# positions = []
# for tx_set in sorted(list(set(trx))):
#     box_data.append([])
#     positions.append(tx_set)
#     for t, tr in enumerate(trx):
#         if tx_set == tr:
#             box_data[-1].append(fp[t])

# bplot = axs.boxplot(box_data, positions=positions, widths=10)