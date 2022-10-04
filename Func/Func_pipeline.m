function [market, ev, ess] = Func_pipeline()
% read
market = readtable('german_SBM_data.csv'); % market data
resource = readtable('resources.csv'); % resource data

%% preprocess
% def ev and ess
ev = resource(resource.type==1,:);
ess = resource(resource.type==2,:);

% [%] -> [p.u.]    
ev.initialSOC = ev.initialSOC / 100;
ev.goalSOC = ev.goalSOC / 100;
ev.minSOC = ev.minSOC / 100;
ev.maxSOC = ev.maxSOC / 100;
ess.initialSOC = ess.initialSOC / 100;
ess.goalSOC = ess.goalSOC / 100;
ess.minSOC = ess.minSOC / 100;
ess.maxSOC = ess.maxSOC / 100;

% plug-time, now = 1 not 0
ev.in = ev.in + 1;
ev.out = ev.out + 1;
ev.duration = ev.out - ev.in + 1;

% EV goalSOC feasibility
goalSOC_feasibility = ev.duration.*ev.pcs - ev.capacity .* (ev.goalSOC - ev.initialSOC);
idx = goalSOC_feasibility < 0;
ev.goalSOC(idx) = ev.goalSOC(idx) + goalSOC_feasibility(idx)./ev.capacity(idx);

T = 24; % 스케줄링의 총 슬롯 수
ess.duration = ones(height(ess),1)*T;
end