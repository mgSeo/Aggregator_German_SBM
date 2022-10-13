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
    goalsoc_ev = model.addVars(range(len(ev)), lb=0)

    # ess Power {vdx,m,t}
    p_ess = [[model.addVars(range(ess.duration[vdx]), ub=ess.pcs[vdx], lb=0, vtype = GRB.INTEGER) for m in range(N['pool'])] for vdx in range(len(ess))]
    # ess UC {vdx,m,t}
    uc_ess = [[model.addVars(range(ess.duration[vdx]), ub=1, lb=0, vtype = GRB.INTEGER) for m in range(N['pool'])] for vdx in range(len(ess))]
    # ess goal soc
    goalsoc_ess = model.addVars(range(len(ess)), lb=0)

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
    for m in range(N['market']):
        for t in range(T):
            backup_ess = 0
            for vdx in range(len(ess)):
                backup_ess += p_ess[vdx][m][t]
            backup_ev = 0
            for vdx in range(len(ev)):
                if ev.arr[vdx] <= t <= ev.dep[vdx]:
                    idx = t - ev.arr[vdx]
                    backup_ev += p_ev[vdx][m][idx]
            model.addConstr(backup_ess == backup_ev*backup_rate)

    ### objective function
    # energy bid
    obj_profit_E = 0    
    for m in range(N['market']):
        for vdx in range(len(ev)):
            for idx in range(ev.duration[vdx]):
                t = ev.arr[vdx] + idx
                obj_profit_E += p_ev[vdx][m][idx] * market['S_'+pools[m]][t] * market['E_'+pools[m]][t]    
        for vdx in range(len(ess)):
            for idx in range(ess.duration[vdx]):
                t = ess.arr[vdx] + idx
                obj_profit_E += p_ess[vdx][m][idx] * market['S_'+pools[m]][t] * market['E_'+pools[m]][t]

    # capacity bid
    obj_profit_C = 0
    for m in range(N['market']):
        for vdx in range(len(ev)):
            for idx in range(ev.duration[vdx]):
                t = ev.arr[vdx] + idx
                obj_profit_C += p_ev[vdx][m][idx] * market['C_'+pools[m]][t]
        for vdx in range(len(ess)):
            for idx in range(ess.duration[vdx]):
                t = ess.arr[vdx] + idx
                obj_profit_C += p_ess[vdx][m][idx] * market['C_'+pools[m]][t]

    # charging cost
    obj_cost_charge = 0
    m = N['market']
    for vdx in range(len(ev)):
        for idx in range(ev.duration[vdx]):
            t = ev.arr[vdx] + idx
            obj_cost_charge += p_ev[vdx][m][idx] * market['sc'][t]
    for vdx in range(len(ess)):
        for idx in range(ess.duration[vdx]):
            t = ess.arr[vdx] + idx
            obj_cost_charge += p_ess[vdx][m][idx] * market['sc'][t]

    # goal soc
    obj_slack_goalsoc = 0
    for vdx in range(len(ev)): obj_slack_goalsoc += goalsoc_ev[vdx]
    for vdx in range(len(ess)): obj_slack_goalsoc += goalsoc_ess[vdx]    

    model.ModelSense = GRB.MINIMIZE
    model.setObjective(-obj_profit_E - obj_profit_C + obj_cost_charge + obj_slack_goalsoc*9999)
    model.optimize()
    ###############################
    ### header: power, uc
    header_p = []
    header_uc = []
    header_soc = []
    for vdx in range(len(ev)):
        for m in range(N['pool']):
            header_p.append('evp{}_'.format(vdx)+pools[m])
            header_uc.append('evuc{}_'.format(vdx)+pools[m])
        header_soc.append('ev{}'.format(vdx))
    for vdx in range(len(ess)):
        for m in range(N['pool']):
            header_p.append('essp{}_'.format(vdx)+pools[m])
            header_uc.append('essuc{}_'.format(vdx)+pools[m])
        header_soc.append('ess{}'.format(vdx))
    # power & uc
    Power = pd.DataFrame(np.zeros([T,len(header_p)]),columns=header_p)
    UC = pd.DataFrame(np.zeros([T,len(header_uc)]),columns=header_uc)
    SoC = pd.DataFrame(np.zeros([T+1,len(header_soc)]),columns=header_soc)
    for vdx in range(len(ev)):
        for m in range(N['pool']):
            for idx in range(ev.duration[vdx]):
                t = ev.arr[vdx]+idx
                Power['evp{}_'.format(vdx)+pools[m]][t] = p_ev[vdx][m][idx].x
                UC['evuc{}_'.format(vdx)+pools[m]][t] = uc_ev[vdx][m][idx].x
    for vdx in range(len(ess)):
        for m in range(N['pool']):
            for idx in range(ess.duration[vdx]):
                t = ess.arr[vdx]+idx
                Power['essp{}_'.format(vdx)+pools[m]][t] = p_ess[vdx][m][idx].x
                UC['essuc{}_'.format(vdx)+pools[m]][t] = uc_ess[vdx][m][idx].x
    # SoC
    for vdx in range(len(ev)):
        soc = ev.initialSOC[vdx]*ev.capacity[vdx]
        SoC['ev{}'.format(vdx)][ev.arr[vdx]] = ev.initialSOC[vdx] * 100
        for idx in range(ev.duration[vdx]):
            t = ev.arr[vdx] + idx
            for m in range(N['market']):
                soc += p_ev[vdx][0][idx].x * market['S_'+pools[m]][t]
            soc += p_ev[vdx][N['market']][idx].x
            SoC['ev{}'.format(vdx)][t+1] = soc/ev.capacity[vdx]*100

    # soc: ess
    for vdx in range(len(ess)):
        soc = ess.initialSOC[vdx]*ess.capacity[vdx]
        SoC['ess{}'.format(vdx)][ess.arr[vdx]] = ess.initialSOC[vdx] * 100
        for idx in range(ess.duration[vdx]):
            t = ess.arr[vdx] + idx
            for m in range(N['market']):
                soc += p_ess[vdx][0][idx].x * market['S_'+pools[m]][t]
            soc += p_ess[vdx][N['market']][idx].x
            SoC['ess{}'.format(vdx)][t+1] = soc/ess.capacity[vdx]*100        

    ### header: bid
    header_bid = []    
    for m in range(N['pool']): header_bid.append('ev_'+pools[m])
    for m in range(N['pool']): header_bid.append('ess_'+pools[m])
    Bid = pd.DataFrame(np.zeros([T,len(header_bid)],))
    # bid-ev
    Bid = pd.DataFrame(np.zeros([T,len(header_bid)]),columns=header_bid)    
    for vdx in range(len(ev)):
        for m in range(N['pool']):
            Bid['ev_'+pools[m]] += Power['evp{}_'.format(vdx)+pools[m]]
    # bid-ess
    for vdx in range(len(ess)):
        for m in range(N['pool']):
            Bid['ess_'+pools[m]] += Power['essp{}_'.format(vdx)+pools[m]]

    return Power, UC, Bid, SoC