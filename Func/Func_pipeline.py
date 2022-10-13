def func(market, resource, signal):
    # [%] -> [p.u.]    
    resource['initialSOC'] /= 100
    resource['goalSOC'] /= 100
    resource['minSOC'] /= 100
    resource['maxSOC'] /= 100

    # def ev, ess
    ev  = resource[resource['type']==1]
    ess = resource[resource['type']==2]

    # plug-time, now = 1 not 0 
    idx = ev[ev.dep<ev.arr].index.to_list()
    ev['dep'][idx] += 24
    idx = ev[ev.dep>23].index.to_list()
    ev['dep'][idx] = 23
    ev['duration'] = ev.dep - ev.arr + 1        
    ev.drop(ev[ev['duration'] < 3].index, inplace=True)

    # EV goalSOC feasibility
    goalSOC_feasibility = ev.duration * ev.pcs - ev.capacity * (ev.goalSOC - ev.initialSOC)
    idx = goalSOC_feasibility[goalSOC_feasibility < 0].index.to_list()
    if idx:
        ev['goalSOC'][idx] = ev.goalSOC[idx] + goalSOC_feasibility[idx]/ev.capacity[idx]

    T = 24; # 스케줄링의 총 슬롯 수
    ess['duration'] = T

    # market unit, % -> p.u.
    market['S_FCR'] /= 100
    market['Sr_FCR'] /= 100
    market['S_aFRR_p'] /= 100
    market['S_aFRR_n'] /= 100

    signal /= 100


    

    return market, ev, ess