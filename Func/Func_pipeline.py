def func(market, resource):
    # [%] -> [p.u.]    
    resource['initialSOC'] /= 100
    resource['goalSOC'] /= 100
    resource['minSOC'] /= 100
    resource['maxSOC'] /= 100

    # def ev, ess
    ev  = resource[resource['type']==1]
    ess = resource[resource['type']==2]

    # plug-time, now = 1 not 0    
    ev['duration'] = ev.dep - ev.arr + 1

    # EV goalSOC feasibility
    goalSOC_feasibility = ev.duration * ev.pcs - ev.capacity * (ev.goalSOC - ev.initialSOC)
    idx = goalSOC_feasibility[goalSOC_feasibility < 0].index.to_list()
    if idx:
        ev['goalSOC'][idx] = ev.goalSOC[idx] + goalSOC_feasibility[idx]/ev.capacity[idx]

    T = 24; # 스케줄링의 총 슬롯 수
    ess['duration'] = T

    return market, ev, ess