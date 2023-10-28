import gurobipy as gp
from gurobipy import GRB
from gurobipy import quicksum

########################
# Solve by Oana's Algo 2
########################
def solve_by_oana1(seqs, m, params):

    t_slots = params['t_slots']
    num_frg = params['num_frg']

    # force complete matching
    # matching_threshold = num_frg

    # output
    T = []

    # iterate over time
    for t in range(t_slots - num_frg):

        # find sequences at t
        # found_seqs = []
        for s, seq in enumerate(seqs):
            # matching = 0
            
            # # explore this seq into the future
            # for ts, p in enumerate(range(num_frg)): 
            #     if m[t + ts, seq[p]] == 1:
            #         matching += 1

            # # store if matching higher than threshold
            # if matching >= matching_threshold:
            #     found_seqs.append([s, matching])

            # if seq can start at t
            seq_fits_m = True
            for ts, p in enumerate(range(num_frg)): 
                if m[t + ts, seq[p]] == -1:
                    seq_fits_m = False
                    break
            if seq_fits_m:
                T.append((t, s, num_frg))

        # # get best starting sequences
        # if found_seqs:
        #     best_seq = max(found_seqs, key=lambda x: x[1])
        #     print(t, best_seq)

        #     T.append((t, best_seq[0], num_frg))

    # minDur, maxDur, observedMatrix
    # observedTraffic = []
    # for t in m until t - minDur
    #     foundSequences = []
    #         for seq in seqs:
    #             matching = 0
    #             for j to j < maxDur
    #                 if observedMatrix[time+j][seq[time+j]] == 1
    #                     matching ++
    #             if matching/maxDuration >= threshold
    #                 foundSequences.append(seq) [(s, maxDuration, matching)]
    #     if foundSequences:
    #         bestSequence = [] <- single foundSequences with highest matching
    #         observedTraffic.append()
    # return observedTraffic

    return T    

#######################
# Solve by MILP1
#######################
def solve_by_milp1(seqs, m, params):

    num_chn = params['num_chn']
    t_slots = params['t_slots']
    num_frg = params['num_frg']

    # output
    T = []

    # prevent gurobi output
    env = gp.Env(empty=True)
    env.setParam('OutputFlag', 0)
    env.start()

    # solve milp model
    try:
        # create model
        model = gp.Model("milp1", env=env)

        # variable: create M_{t,c}
        Mvars = {}
        for t in range(t_slots):
            for c in range(num_chn):
                if m[t, c] > -1:
                    Mvars[t, c] = model.addVar(vtype=GRB.BINARY, name="M.{}.{}".format(t, c))

        # variable: create T_{t,s} 
        Tvars = {}
        for t in range(t_slots):
            for s, seq in enumerate(seqs):
                # if seq can start at t
                seq_fits_m = True
                for ts, p in enumerate(range(num_frg)): 
                    if m[t + ts, seq[p]] == -1:
                        seq_fits_m = False
                        break
                if seq_fits_m:
                    # add T_{t,s} variable
                    Tvars[t, s] = model.addVar(vtype=GRB.BINARY, name="T.{}.{}".format(t, s))

        # constraint: FRG * T_{t,s} = Sum M_{t,c} corresponding to t,s
        for key in Tvars:
            t = key[0]
            s = key[1]
            expr1 = num_frg * Tvars[t, s]
            expr2 = 0
            for ts, p in enumerate(range(num_frg)): 
                # use M variables
                # expr2 += Mvars[t + ts, seqs[s][p]]
                # use m variables (all 1)
                expr2 += m[t + ts, seqs[s][p]]
            model.addConstr(expr1 == expr2, "cT.{}.{}".format(t, s))
        
        # set objective
        model.setObjective(quicksum(list(Tvars.values())), GRB.MINIMIZE)

        # optimize model
        model.optimize()

        # print results
        for v in model.getVars():
            # print("{} {}".format(v.VarName, v.X))
            var_list = v.VarName.split('.')
            if var_list[0] == 'T':
                t = int(var_list[1])
                s = int(var_list[2])
                if v.X == 1:
                    T.append((t, s, num_frg))
                    # print('added {}'.format(T[-1]))

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    except AttributeError:
        print('Encountered an attribute error')

    return T
