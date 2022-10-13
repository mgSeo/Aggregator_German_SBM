import pandas as pd
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))

market = pd.read_csv('input/german_SBM_data_py.csv') # read_csv
resource = pd.read_csv('input/resources_py.csv')

from Func import Func_pipeline # data processing
market, ev, ess = Func_pipeline.func(market, resource)

from Func import SBM_bid_model # SBM scheduling
backup_rate = 0.2
power, uc, bid, soc = SBM_bid_model.func(market, ev, ess, backup_rate)

power.to_csv('power.csv',index=False)
uc.to_csv('uc.csv',index=False)
bid.to_csv('bid.csv',index=False)
soc.to_csv('soc.csv',index=False)