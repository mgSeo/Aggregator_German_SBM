from pydoc import resolve
import pandas as pd
import sys, os
sys.path.append(os.path.dirname(os.path.abspath(os.path.dirname(__file__))))

market = pd.read_csv('input/german_SBM_data_py.csv') # read_csv
resource = pd.read_csv('input/resources_py.csv')

from Func import Func_pipeline # data processing
market, ev, ess = Func_pipeline.func(market, resource)

from Func import SBM_bid_model # SBM scheduling
backup_rate = 0.2
Power_ev, UC_ev, Power_ess, UC_ess = SBM_bid_model.func(market, ev, ess, backup_rate)

Power_ev.to_csv('power_ev.csv',index=False)
UC_ev.to_csv('uc_ev.csv',index=False)
Power_ess.to_csv('power_ess.csv',index=False)
UC_ess.to_csv('uc_ess.csv',index=False)