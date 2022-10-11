def func(market, ev, ess, backup_rate):
    import pandas as pd
    import numpy as np
    import gurobipy as gp
    from gurobipy import GRB

    model = gp.Model('SBM bid aggregator')
    # model.setParam('OutputFlag', False)
    ### Define variable:
    N = {'market':3, 'sc':1, 'pool':4}
    T = 24
    ev = ev.reset_index()
    ess = ess.reset_index()
    pools = ['FCR','aFRR_p','aFRR_n','sc']

    # ev Power {vdx,m,t}
    p_ev = [[model.addVars(range(ev.duration[vdx]), ub=ev.pcs[vdx], lb=0, vtype = GRB.INTEGER) for m in range(N['pool'])] for vdx in range(len(ev))]
    # ev UC {vdx,m,t}
    uc_ev = [[model.addVars(range(ev.duration[vdx]), ub=1, lb=0, vtype = GRB.INTEGER) for m in range(N['pool'])] for vdx in range(len(ev))]
    # ev goal soc
    goalsoc_ev = model.addVars(range(len(ev)), lb=0, vtype = GRB.INTEGER)

    # ess Power {vdx,m,t}
    p_ess = [[model.addVars(range(ess.duration[vdx]), ub=ess.pcs[vdx], lb=0, vtype = GRB.INTEGER) for m in range(N['pool'])] for vdx in range(len(ess))]
    # ess UC {vdx,m,t}
    uc_ess = [[model.addVars(range(ess.duration[vdx]), ub=1, lb=0, vtype = GRB.INTEGER) for m in range(N['pool'])] for vdx in range(len(ess))]
    # ess goal soc
    goalsoc_ess = model.addVars(range(len(ess)), lb=0, vtype = GRB.INTEGER)

    ### Constraints:
    # uc: ev
    for vdx in range(len(ev)):
        for idx in range(ev.duration[vdx]):
            uc = 0
            for m in range(N['pool']):
                uc += uc_ev[vdx][m][idx]
            model.addConstr(uc <= 1)
    # uc: ess
    for vdx in range(len(ess)):
        for idx in range(ess.duration[vdx]):
            uc = 0
            for m in range(N['pool']):
                uc += uc_ess[vdx][m][idx]
            model.addConstr(uc <= 1)

    # uc-p: ev
    for vdx in range(len(ev)):
        for idx in range(ev.duration[vdx]):
            for m in range(N['pool']):
                model.addConstr(p_ev[vdx][m][idx] <= uc_ev[vdx][m][idx]*ev.pcs[vdx])
    # uc-p: ess
    for vdx in range(len(ess)):
        for idx in range(ess.duration[vdx]):
            for m in range(N['pool']):
                model.addConstr(p_ess[vdx][m][idx] <= uc_ess[vdx][m][idx]*ess.pcs[vdx])

    # soc: ev
    for vdx in range(len(ev)):
        soc = ev.initialSOC[vdx]*ev.capacity[vdx]
        for idx in range(ev.duration[vdx]):
            t = ev.arr[vdx] + idx
            for m in range(N['market']):
                soc += p_ev[vdx][0][idx] * market['S_'+pools[m]][t]
            soc += p_ev[vdx][N['market']][idx]
            model.addConstr(soc <= ev.maxSOC[vdx]*ev.capacity[vdx])
            model.addConstr(soc >= ev.minSOC[vdx]*ev.capacity[vdx])
        # goal soc: ev
        model.addConstr(goalsoc_ev[vdx] >= ev.goalSOC[vdx]*ev.capacity[vdx] - soc)
        model.addConstr(goalsoc_ev[vdx] >= soc - ev.goalSOC[vdx]*ev.capacity[vdx])

    # soc: ess
    for vdx in range(len(ess)):
        soc = ess.initialSOC[vdx]*ess.capacity[vdx]
        for idx in range(ess.duration[vdx]):
            t = ess.arr[vdx] + idx
            for m in range(N['market']):
                soc += p_ess[vdx][0][idx] * market['S_'+pools[m]][t]
            soc += p_ess[vdx][N['market']][idx]
            model.addConstr(soc <= ess.maxSOC[vdx]*ess.capacity[vdx])
            model.addConstr(soc >= ess.minSOC[vdx]*ess.capacity[vdx])
        # goal soc: ess
        model.addConstr(goalsoc_ess[vdx] >= ess.goalSOC[vdx]*ess.capacity[vdx] - soc)
        model.addConstr(goalsoc_ess[vdx] >= soc - ess.goalSOC[vdx]*ess.capacity[vdx])
    
    # cbl: ev
    for vdx in range(len(ev)):
        for idx in range(ev.duration[vdx]-1):
            for m in range(N['pool']):                
                cbl_up = uc_ev[vdx][m][idx+1] - uc_ev[vdx][m][idx]
                m_list = list(range(N['pool']))
                m_list.remove(m)
                for mm in m_list:
                    cbl_dn = uc_ev[vdx][m][idx] - uc_ev[vdx][m][idx+1]
                    model.addConstr(cbl_up + cbl_dn <= 1)
    # cbl: ess
    for vdx in range(len(ess)):
        for idx in range(ess.duration[vdx]-1):
            for m in range(N['pool']):                
                cbl_up = uc_ess[vdx][m][idx+1] - uc_ess[vdx][m][idx]
                m_list = list(range(N['pool']))
                m_list.remove(m)
                for mm in m_list:
                    cbl_dn = uc_ess[vdx][m][idx] - uc_ess[vdx][m][idx+1]
                    model.addConstr(cbl_up + cbl_dn <= 1)

    # 4 hour bidding slot: ev
    for m in range(N['market']):
        for s in range(6):
            i = 0
            t = 4*s + i
            p_front = 0
            for vdx in range(len(ev)):
                if ev.arr[vdx] <= t <= ev.dep[vdx]:
                    idx = t - ev.arr[vdx]
                    p_front += p_ev[vdx][m][idx]
            for i in range(1,4):
                t = 4*s + i
                p_behind = 0
                for vdx in range(len(ev)):
                    if ev.arr[vdx] <= t <= ev.dep[vdx]:
                        idx = t - ev.arr[vdx]
                        p_behind += p_ev[vdx][m][idx]
                model.addConstr(p_front == p_behind)
                p_front = p_behind
    # 4 hour bidding slot: ess
    for m in range(N['market']):
        for s in range(6):
            i = 0
            t = 4*s + i
            p_front = 0
            for vdx in range(len(ess)):
                idx = t - ess.arr[vdx]
                p_front += p_ess[vdx][m][idx]
            for i in range(1,4):
                t = 4*s + i
                p_behind = 0
                for vdx in range(len(ess)):
                    idx = t - ess.arr[vdx]
                    p_behind += p_ess[vdx][m][idx]
                model.addConstr(p_front == p_behind)
                p_front = p_behind

    # back-up resource slot

    ### objective function
    # energy bid
    obj_E = 0
    m = 0
    for m in range(N['market']):
        for vdx in range(len(ev)):
            for idx in range(ev.duration[vdx]):
                t = ev.arr[vdx] + idx
                obj_E += p_ev[vdx][m][idx] * market['S_'+pools[m]][t] * market['E_'+pools[m]][t]
    # capacity bid
    obj_C = 0
    m = 0
    for m in range(N['market']):
        for vdx in range(len(ev)):
            for idx in range(ev.duration[vdx]):
                t = ev.arr[vdx] + idx
                obj_E += p_ev[vdx][m][idx] * market['C_'+pools[m]][t]

    # charging cost
    
    model.ModelSense = GRB.MINIMIZE
    model.setObjective(-obj_E)    
    model.optimize()
    ###############################
    # Power
    header_p = [] 
    header_uc = []
    for vdx in range(len(ev)):
        for m in range(N['pool']):
            header_p.append('p{}_'.format(vdx)+pools[m])
            header_uc.append('uc{}_'.format(vdx)+pools[m])
    # ev
    Power_ev = pd.DataFrame(np.zeros([T,len(header_p)]),columns=header_p)
    UC_ev = pd.DataFrame(np.zeros([T,len(header_uc)]),columns=header_uc)
    for vdx in range(len(ev)):
        for m in range(N['pool']):
            for idx in range(ev.duration[vdx]):
                t = ev.arr[vdx]+idx
                Power_ev['p{}_'.format(vdx)+pools[m]][t] = p_ev[vdx][m][idx].x
                UC_ev['uc{}_'.format(vdx)+pools[m]][t] = uc_ev[vdx][m][idx].x

    # ess
    Power_ess = pd.DataFrame(np.zeros([T,len(header_p)]),columns=header_p)
    UC_ess = pd.DataFrame(np.zeros([T,len(header_uc)]),columns=header_uc)
    for vdx in range(len(ess)):
        for m in range(N['pool']):
            for idx in range(ess.duration[vdx]):
                t = ess.arr[vdx]+idx
                Power_ess['p{}_'.format(vdx)+pools[m]][t] = p_ess[vdx][m][idx].x
                UC_ess['uc{}_'.format(vdx)+pools[m]][t] = uc_ess[vdx][m][idx].x
                

    return Power_ev, UC_ev, Power_ess, UC_ess