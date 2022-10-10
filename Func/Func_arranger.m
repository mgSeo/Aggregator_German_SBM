function [SoC, R, Bid] = Func_arranger(x,market,ev,ess)
%% parameters:
T = 24;
N.market = 3;
N.sc = 1;
N.ev = height(ev);
N.ess = height(ess);
N.dur_ev = sum(ev.duration);
N.dur_ess = sum(ess.duration);
N.uc_ev = N.dur_ev * (N.market+N.sc);
N.uc_ess = N.dur_ess * (N.market+N.sc);
N.p_ev = N.dur_ev * (N.market+N.sc);
N.p_ess = N.dur_ess * (N.market+N.sc);
lenN = N.p_ev + N.p_ess + N.uc_ev + N.uc_ess + N.dur_ev + N.dur_ess;
%%
SoC = table();
% soc-ev
A_soc_ev = zeros(N.dur_ev,lenN);
B_soc_ev_ub = repelem((ev.maxSOC-ev.initialSOC).*ev.capacity,ev.duration,1);
B_soc_ev_lb = repelem((ev.initialSOC-ev.minSOC).*ev.capacity,ev.duration,1);
temp.d = 0; temp.v = 0;
for i = 1:N.ev
    dur = ev.in(i):ev.out(i);
    mat = tril(ones(ev.duration(i)));
    dev = [mat.*market.S_FCR(dur) mat.*market.S_aFRR_p(dur) mat.*market.S_aFRR_n(dur) mat];
    A_soc_ev(temp.d+1 : temp.d+ev.duration(i), temp.v+1 : temp.v+ev.duration(i)*(N.market+N.sc))  = dev;
    temp.d = temp.d + ev.duration(i);
    temp.v = temp.v + ev.duration(i)*(N.market+N.sc);
end
soe = A_soc_ev*x + repelem(ev.initialSOC.*ev.capacity,ev.duration,1);
for i = 1:N.ev
    mat = zeros(24,1);
    mat(ev.in(i):ev.out(i))= soe(1:ev.duration(i))/ev.capacity(i)*100;
    id = ['ev_id_',num2str(ev.id(i))];
    SoC(:,i) = array2table(mat);
    SoC.Properties.VariableNames{i} = id;
    soe(1:ev.duration(i)) = [];
end
% soc-ess
A_soc_ess = zeros(N.dur_ess,lenN);
temp.d = 0; temp.v = N.uc_ev;
for i = 1:N.ess
    dur = 1:T;
    mat = tril(ones(ess.duration(i)));
    dev = [mat.*market.S_FCR(dur) mat.*market.S_aFRR_p(dur) mat.*market.S_aFRR_n(dur) mat];
    A_soc_ess(temp.d+1 : temp.d+ess.duration(i), temp.v+1 : temp.v+ess.duration(i)*(N.market+N.sc))  = dev;
    temp.d = temp.d + ess.duration(i);
    temp.v = temp.v + ess.duration(i)*(N.market+N.sc);
end
soe = A_soc_ess*x + repelem(ess.initialSOC.*ess.capacity,ess.duration,1);
for i = 1:N.ess
    mat = soe(1:ess.duration(i))/ess.capacity(i)*100;
    id = ['ess_id_',num2str(ess.id(i))];
    SoC(:,N.ev+i) = array2table(mat);
    SoC.Properties.VariableNames{N.ev+i} = id;
    soe(1:ess.duration(i)) = [];
end

%% UC-P
R = table();
% power
Market = {'FCR_ev_id_', 'aFRRp_ev_id_', 'aFRRn_ev_id_', 'sc_ev_id_'};
temp.i = 0;

for i = 1:N.ev
    temp.d = 0;
    for r = 1:N.market + N.sc
        mat = zeros(24,1);
        mat(ev.in(i):ev.out(i)) = round(x(temp.d+1 : temp.d+ev.duration(i)));
        temp.d = temp.d + ev.duration(i);

        id = [Market{r},num2str(ev.id(i))];
        R(:,temp.i+r) = array2table(mat);
        R.Properties.VariableNames{temp.i+r} = id;
    end
    x(1:ev.duration(i)*4) = [];
    temp.i = temp.i + 4;
end

Market = {'FCR_ess_id_', 'aFRRp_ess_id_', 'aFRRn_ess_id_', 'sc_ess_id_'};
for i = 1:N.ess
    temp.d = 0;
    for r = 1:N.market + N.sc
        mat = round(x(temp.d+1 : temp.d+ess.duration(i)));
        temp.d = temp.d + ess.duration(i);

        id = [Market{r},num2str(ess.id(i))];
        R(:,temp.i+r) = array2table(mat);
        R.Properties.VariableNames{temp.i+r} = id;
    end
    x(1:ess.duration(i)*4) = [];
    temp.i = temp.i + 4;
end


% uc
Market = {'uc_FCR_ev_id_', 'uc_aFRRp_ev_id_', 'uc_aFRRn_ev_id_', 'uc_sc_ev_id_'};
for i = 1:N.ev
    temp.d = 0;
    for r = 1:N.market + N.sc
        mat = zeros(24,1);
        mat(ev.in(i):ev.out(i)) = round(x(temp.d+1 : temp.d+ev.duration(i)));
        temp.d = temp.d + ev.duration(i);

        id = [Market{r},num2str(ev.id(i))];
        R(:,temp.i+r) = array2table(mat);
        R.Properties.VariableNames{temp.i+r} = id;
    end
    x(1:ev.duration(i)*4) = [];
    temp.i = temp.i + 4;
end

Market = {'uc_FCR_ess_id_', 'uc_aFRRp_ess_id_', 'uc_aFRRn_ess_id_', 'uc_sc_ess_id_'};
for i = 1:N.ess
    temp.d = 0;
    for r = 1:N.market + N.sc
        mat = round(x(temp.d+1 : temp.d+ess.duration(i)));
        temp.d = temp.d + ess.duration(i);
        id = [Market{r},num2str(ess.id(i))];
        R(:,temp.i+r) = array2table(mat);
        R.Properties.VariableNames{temp.i+r} = id;
    end
    x(1:ess.duration(i)*4) = [];
    temp.i = temp.i + 4;
end

% reset sequence
Market = {'FCR_ev_id_', 'aFRRp_ev_id_', 'aFRRn_ev_id_', 'sc_ev_id_', 'uc_FCR_ev_id_', 'uc_aFRRp_ev_id_', 'uc_aFRRn_ev_id_', 'uc_sc_ev_id_'};
for i = 1:N.ev % table column sorting
    for n = 1:7
        id_left = [Market{n},num2str(ev.id(i))];
        id_right = [Market{n+1},num2str(ev.id(i))];
        R = movevars(R, id_right, 'after', id_left);
    end
end

%% Bid
%% Bid
% load('Bid.mat')
% result Bid
bids = zeros(T,4);
for m = 1:4
    for vdx = 1:N.ev
        bids(:,m) = bids(:,m) + R{:,(vdx-1)*8+m};
    end
end

FCR = bids(:,1);
aFRRp = bids(:,2);
aFRRn = bids(:,3);
sc = bids(:,4);
Bid = table(FCR,aFRRp,aFRRn,sc);

end