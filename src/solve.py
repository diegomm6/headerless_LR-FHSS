import gurobipy as gp
from gurobipy import GRB
from gurobipy import quicksum

#######################
# Solve by MILP1
#######################
def solve_by_milp1(seqs, m, params):

    num_chn = params['num_chn']
    t_slots = params['t_slots']
    num_frg = params['num_frg']
    num_hdr = params['num_hdr']
    gran = params['granularity']

    max_seq_len = num_hdr * round(gran * 233 / 102.4) + gran * num_frg

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
                if m[t][c] > -1:
                    Mvars[t, c] = model.addVar(vtype=GRB.BINARY, name="M.{}.{}".format(t, c))


        Tvars = {}
        for t in range(t_slots - max_seq_len):
            for s, seq in enumerate(seqs):

                time = t
                seq_fits_m = True # if seq can start at t

                for fh, chn in enumerate(seq):

                    used_slots = gran # write fragment
                    if fh < num_hdr:  # write header
                        used_slots = round(gran * 233 / 102.4)
                    
                    for g in range(used_slots):
                        if m[time + g][chn] == -1:
                            seq_fits_m = False
                            break

                    if not seq_fits_m:
                        break

                    time += used_slots


                if seq_fits_m:
                    # add T_{t,s} variable
                    Tvars[t, s] = model.addVar(vtype=GRB.BINARY, name="T.{}.{}".format(t, s))

        # constraint: FRG * T_{t,s} = Sum M_{t,c} corresponding to t,s
        for key in Tvars:
            t = key[0]
            s = key[1]

            # expr1 = seq_len * Tvars[t, s]   # modify for extended model

            s_len = num_hdr * round(gran * 233 / 102.4) + gran * (len(seqs[s]) - num_hdr)
            expr1 = s_len * Tvars[t, s]

            expr2 = 0
            time = t
            for fh, chn in enumerate(seqs[s]): 

                used_slots = gran # fragment
                if fh < num_hdr:  # header
                    used_slots = round(gran * 233 / 102.4)

                for g in range(used_slots):
                    expr2 += m[time + g][chn]
                
                time += used_slots

            model.addConstr(expr1 == expr2, "cT.{}.{}".format(t, s))

        # set objective
        model.setObjective(quicksum(list(Tvars.values())), GRB.MINIMIZE)

        # optimize model
        model.optimize()

        # print results
        for v in model.getVars():
            #print("{} {}".format(v.VarName, v.X))
            var_list = v.VarName.split('.')
            if var_list[0] == 'T':
                t = int(var_list[1])
                s = int(var_list[2])
                if v.X == 1:
                    # T.append((t, s, num_frg))  # modify for extended model
                    T.append((t, s, len(seqs[s])))
                    # print('added {}'.format(T[-1]))

    except gp.GurobiError as e:
        print('Error code ' + str(e.errno) + ': ' + str(e))

    #except AttributeError:
    #    print('Encountered an attribute error')

    return T
