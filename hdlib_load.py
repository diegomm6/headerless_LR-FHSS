import csv

#######################
# Load sequence file
#######################
def load_seqs(file_name):

    # sequences (9 bits -> up to 512)
    seqs = []

    with open(file_name, 'r') as file:
        csv_reader = csv.reader(file, delimiter=',')
        for r, row in enumerate(csv_reader):
            seqs.append([])
            for c, col in enumerate(row):
                if col != '':
                    seqs[-1].append(int(col))

    return seqs


#######################
# Load traffic file
#######################
def load_T(file_name):

    T = [] 
    with open(file_name, 'r') as file:
        csv_reader = csv.reader(file, delimiter=',')
        for r, row in enumerate(csv_reader):
            T.append((int(row[0]), int(row[1]), int(row[2])))
            
    return T