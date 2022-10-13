import pandas as pd
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))

market = pd.read_csv('input/german_SBM_data_py.csv') # read_csv
resource = pd.read_csv('input/resources_py.csv')
signal = pd.read_csv('input/scenario.csv')

from Func import Func_pipeline # data processing
market, ev, ess = Func_pipeline.func(market, resource)

from Func import SBM_bid_model # SBM scheduling
backup_rate = 0.2
power, uc, bid, soc = SBM_bid_model.func(market, ev, ess, backup_rate)

from Func import SBM_robust_bid_model # SBM scheduling
backup_rate = 0.2
r_power, r_uc, r_bid, r_soc = SBM_robust_bid_model.func(market, ev, ess, backup_rate, signal)

power.to_csv('power.csv',index=False)
uc.to_csv('uc.csv',index=False)
bid.to_csv('bid.csv',index=False)
soc.to_csv('soc.csv',index=False)

r_power.to_csv('robust_power.csv',index=False)
r_uc.to_csv('robust_uc.csv',index=False)
r_bid.to_csv('robust_bid.csv',index=False)
r_soc.to_csv('robust_soc.csv',index=False)