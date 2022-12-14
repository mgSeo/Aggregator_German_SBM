function [x, exitflag] = SBM_bid_model(market, ev, ess)
%% Parameters
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
%     cbl = sum(ESS.num*(T-1)*(M-SC)) + sum(EV.dur-1)*(M-SC); % CBL 제약에 필요한 변수 수
% N_cbl = sumESS + sumEV + EV.num + maxREid*T*2; % CBL 제약 적용을 용이하게 하기위한 간이 식
lenN = N.p_ev + N.p_ess + N.uc_ev + N.uc_ess; % 필요한 결정 변수의 수를 미리 산정 후 코딩 (아래는 결정 변수 목록)
% ev UC / ess UC / ev total duration / ess total duration / CBL 추가필요(start-up, shut-down)

lb = zeros(lenN,1);
ub = ones(lenN,1);
ub(1 : N.p_ev) = repelem(ev.pcs,ev.duration*(N.market+N.sc),1); % capacity bid
ub(N.p_ev+1 : N.p_ev+N.p_ess) = repelem(ess.pcs,ess.duration*(N.market+N.sc),1); % capacity bid
intcon = 1:lenN;
%% Constraints for market participation
%%% UC
% UC-ev
A_uc_ev = zeros(N.dur_ev,lenN);
temp.d = 0; temp.v = N.p_ev + N.p_ess;
for vdx = 1:N.ev
    A_uc_ev(temp.d+1 : temp.d+ev.duration(vdx), temp.v+1 : temp.v+ev.duration(vdx)*(N.market+N.sc)) = repmat(eye(ev.duration(vdx)),1,(N.market+N.sc));
    temp.d = temp.d + ev.duration(vdx);
    temp.v = temp.v + ev.duration(vdx)*(N.market+N.sc);
end
B_uc_ev = ones(N.dur_ev,1);

% UC-ess
A_uc_ess = zeros(N.dur_ess,lenN);
temp.d = 0; temp.v = N.p_ev + N.p_ess + N.uc_ev;
for vdx = 1:N.ess
    A_uc_ess(temp.d+1 : temp.d+ess.duration(vdx), temp.v+1 : temp.v+ess.duration(vdx)*(N.market+N.sc)) = repmat(eye(ess.duration(vdx)),1,(N.market+N.sc));
    temp.d = temp.d + ess.duration(vdx);
    temp.v = temp.v + ess.duration(vdx)*(N.market+N.sc);
end
B_uc_ess = ones(N.dur_ess,1);

%%% UC-P
% p-ev
A_p_uc_ev = zeros(N.p_ev,lenN);
temp.v = 0;
A_p_uc_ev(1 : N.p_ev, temp.v+1:temp.v+N.p_ev) = eye(N.p_ev);
temp.d = 0;
temp.v = N.p_ev + N.p_ess;
for vdx = 1:N.ev
    A_p_uc_ev(temp.d+1 : temp.d+ev.duration(vdx)*(N.market+N.sc), temp.v+1:temp.v+ev.duration(vdx)*(N.market+N.sc)) = -eye(ev.duration(vdx)*(N.market+N.sc))*ev.pcs(vdx);
    temp.d = temp.d + ev.duration(vdx)*(N.market+N.sc);
    temp.v = temp.v + ev.duration(vdx)*(N.market+N.sc);
end
B_p_uc_ev = zeros(N.p_ev,1);

% p-ess
A_p_uc_ess = zeros(N.p_ess,lenN);
temp.v = N.p_ev;
A_p_uc_ess(1 : N.p_ess, temp.v+1:temp.v+N.p_ess) = eye(N.p_ess);
temp.d = 0;
temp.v = N.p_ess + N.p_ess + N.uc_ev;
for vdx = 1:N.ess
    A_p_uc_ess(temp.d+1 : temp.d+ess.duration(vdx)*(N.market+N.sc), temp.v+1:temp.v+ess.duration(vdx)*(N.market+N.sc)) = -eye(ess.duration(vdx)*(N.market+N.sc))*ess.pcs(vdx);
    temp.d = temp.d + ess.duration(vdx)*(N.market+N.sc);
    temp.v = temp.v + ess.duration(vdx)*(N.market+N.sc);
end
B_p_uc_ess = zeros(N.p_ess,1);


%%% SoC boundary #### <= stochastic#######################################################
% SoC-ev
A_soc_ev_ub = zeros(N.dur_ev,lenN);
B_soc_ev_ub = repelem((ev.maxSOC-ev.initialSOC).*ev.capacity,ev.duration,1);
B_soc_ev_lb = repelem((ev.initialSOC-ev.minSOC).*ev.capacity,ev.duration,1);
temp.d = 0; temp.v = 0;
for vdx = 1:N.ev
    dur = ev.in(vdx):ev.out(vdx);
    mat = tril(ones(ev.duration(vdx)));
    dev = [mat.*market.S_FCR(dur) mat.*market.S_aFRR_p(dur) mat.*market.S_aFRR_n(dur) mat];
    A_soc_ev_ub(temp.d+1 : temp.d+ev.duration(vdx), temp.v+1 : temp.v+ev.duration(vdx)*(N.market+N.sc))  = dev;
    temp.d = temp.d + ev.duration(vdx);
    temp.v = temp.v + ev.duration(vdx)*(N.market+N.sc);
end
A_soc_ev_lb = -A_soc_ev_ub;

% SoC-ess
A_soc_ess_ub = zeros(N.dur_ess,lenN);
B_soc_ess_ub = repelem((ess.maxSOC-ess.initialSOC).*ess.capacity,ess.duration,1);
B_soc_ess_lb = repelem((ess.initialSOC-ess.minSOC).*ess.capacity,ess.duration,1);
temp.d = 0; temp.v = N.uc_ev;
for vdx = 1:N.ess
    dur = 1:T;
    mat = tril(ones(ess.duration(vdx)));
    dev = [mat.*market.S_FCR(dur) mat.*market.S_aFRR_p(dur) mat.*market.S_aFRR_n(dur) mat];
    A_soc_ess_ub(temp.d+1 : temp.d+ess.duration(vdx), temp.v+1 : temp.v+ess.duration(vdx)*(N.market+N.sc))  = dev;
    temp.d = temp.d + ess.duration(vdx);
    temp.v = temp.v + ess.duration(vdx)*(N.market+N.sc);
end
A_soc_ess_lb = -A_soc_ess_ub;

%%% Goal SoC
% goalsoc-ev
A_goalsoc_ev = A_soc_ev_lb(cumsum(ev.duration),:);
B_goalsoc_ev = (ev.initialSOC - ev.goalSOC).*ev.capacity;
% goalsoc-ess
A_goalsoc_ess = A_soc_ess_lb(cumsum(ess.duration),:);
B_goalsoc_ess = (ess.initialSOC - ess.goalSOC).*ess.capacity;

%%% CBL
% def start-up, shut-down - ev
cbl_slack_ev = zeros(N.dur_ev - N.ev, lenN);
temp.d1 = 0; temp.d2 = N.p_ev + N.p_ess;
for vdx = 1:N.ev
    cbl_slack_ev(temp.d1+1: temp.d1+ev.duration(vdx)-1,temp.d2+1:temp.d2+ev.duration(vdx)*4) = repmat(eye(ev.duration(vdx)-1,ev.duration(vdx)) - [zeros(ev.duration(vdx)-1,1) eye(ev.duration(vdx)-1)],1,4);
    temp.d1 = temp.d1 + ev.duration(vdx) - 1;
    temp.d2 = temp.d2 + ev.duration(vdx)*4;
end

% def start-up, shut-down - ess
cbl_slack_ess = zeros(N.dur_ess - N.ess, lenN);
temp.d1 = 0; temp.d2 = N.p_ess + N.p_ess;
for vdx = 1:N.ess
    cbl_slack_ess(temp.d1+1: temp.d1+ess.duration(vdx)-1,temp.d2+1:temp.d2+ess.duration(vdx)*4) = repmat(eye(ess.duration(vdx)-1,ess.duration(vdx)) - [zeros(ess.duration(vdx)-1,1) eye(ess.duration(vdx)-1)],1,4);
    temp.d1 = temp.d1 + ess.duration(vdx) - 1;
    temp.d2 = temp.d2 + ess.duration(vdx)*4;
end

% cbl-ev
A_cbl_ev = zeros((N.dur_ev - N.ev)*16, lenN);
B_cbl_ev = ones((N.dur_ev - N.ev)*16, 1);
for R = 1:4    
    MAT = zeros((N.dur_ev - N.ev)*4, lenN);    
    notM = 1:4;
    notM = notM(notM~=R);
    for r = 1:3        
        temp.d1 = 0; temp.d2 = N.p_ev + N.p_ess;
        mat = zeros(N.dur_ev - N.ev, lenN);
        for vdx = 1:N.ev
            temp.r = ev.duration(vdx) * (R-1); % start-up
            mat(temp.d1+1: temp.d1+ev.duration(vdx)-1, temp.d2+temp.r+1 : temp.d2+temp.r+ev.duration(vdx)) ...
                = -cbl_slack_ev(temp.d1+1: temp.d1+ev.duration(vdx)-1, temp.d2+temp.r+1 : temp.d2+temp.r+ev.duration(vdx));
            temp.r = ev.duration(vdx) * (notM(r)-1); % shut-down
            mat(temp.d1+1: temp.d1+ev.duration(vdx)-1, temp.d2+temp.r+1 : temp.d2+temp.r+ev.duration(vdx)) ...
                = cbl_slack_ev(temp.d1+1: temp.d1+ev.duration(vdx)-1, temp.d2+temp.r+1 : temp.d2+temp.r+ev.duration(vdx));
            temp.d1 = temp.d1 + ev.duration(vdx) - 1;
            temp.d2 = temp.d2 + ev.duration(vdx)*4;
        end
        MAT((r-1)*(N.dur_ev - N.ev)+1:r*(N.dur_ev - N.ev), :) = mat;
    end

    A_cbl_ev((R-1)*(N.dur_ev - N.ev)*4+1:R*(N.dur_ev - N.ev)*4, :) = MAT;
end

% cbl-ess
A_cbl_ess = zeros((N.dur_ess - N.ess)*16, lenN);
B_cbl_ess = ones((N.dur_ess - N.ess)*16, 1);
for R = 1:4    
    MAT = zeros((N.dur_ess - N.ess)*4, lenN);    
    notM = 1:4;
    notM = notM(notM~=R);
    for r = 1:3        
        temp.d1 = 0; temp.d2 = N.p_ess + N.p_ess;
        mat = zeros(N.dur_ess - N.ess, lenN);
        for vdx = 1:N.ess
            temp.r = ess.duration(vdx) * (R-1); % start-up
            mat(temp.d1+1: temp.d1+ess.duration(vdx)-1, temp.d2+temp.r+1 : temp.d2+temp.r+ess.duration(vdx)) ...
                = -cbl_slack_ess(temp.d1+1: temp.d1+ess.duration(vdx)-1, temp.d2+temp.r+1 : temp.d2+temp.r+ess.duration(vdx));
            temp.r = ess.duration(vdx) * (notM(r)-1); % shut-down
            mat(temp.d1+1: temp.d1+ess.duration(vdx)-1, temp.d2+temp.r+1 : temp.d2+temp.r+ess.duration(vdx)) ...
                = cbl_slack_ess(temp.d1+1: temp.d1+ess.duration(vdx)-1, temp.d2+temp.r+1 : temp.d2+temp.r+ess.duration(vdx));
            temp.d1 = temp.d1 + ess.duration(vdx) - 1;
            temp.d2 = temp.d2 + ess.duration(vdx)*4;
        end
        MAT((r-1)*(N.dur_ess - N.ess)+1:r*(N.dur_ess - N.ess), :) = mat;
    end

    A_cbl_ess((R-1)*(N.dur_ess - N.ess)*4+1:R*(N.dur_ess - N.ess)*4, :) = MAT;
end


%%% 4 hour bidding slot
% def hourly participation
Pev{24,N.market*N.ev} = [];
cumdur = [0; cumsum(ev.duration)];
for vdx = 1:N.ev
    in = ev.in(vdx);
    out = ev.out(vdx);
    for t = in:out
        temp.m = 0;
        for m = 1:3 % FCR/aFRRp/aFRRn
            mat = zeros(1,lenN);
            mat(cumdur(vdx)*4+t-in+1+temp.m) = 1;
            Pev{t,(vdx-1)*3+m} = mat;
            temp.m = temp.m + ev.duration(vdx);
        end
    end
end

% plug-map
plug_map = zeros(24,1);
for vdx = 1:N.ev
    plug_map(ev.in(vdx):ev.out(vdx)) = plug_map(ev.in(vdx):ev.out(vdx)) + 1;
end
ava_market = zeros(6,1);
for t = 1:6
    if sum(plug_map((t-1)*4+1:t*4) == 0) > 0
        continue;
    end
    ava_market(t) = 1;
end

% 4hour bid-ev
Aeq_4hour_bid_ev = [];
bid_FCR_ev = [];
bid_aFRRp_ev = [];
bid_aFRRn_ev = [];
for t = 1 : 24/4
    in = (t-1)*4 + 1 - min(ev.in) + 1;
    if ava_market(t) == 0
        continue;
    end
    mat = zeros(N.market*4,lenN);
    for m = 1:N.market
        for vdx = 1:N.ev
            for tdx = 1:4
                if isempty(Pev{(t-1)*4+tdx,(vdx-1)*3+m}) == 1
                    continue
                end
                mat(tdx,:) = mat(tdx,:) + Pev{(t-1)*4+tdx,(vdx-1)*3+m};
            end
        end
    end
    Aeq_4hour_bid_ev = [Aeq_4hour_bid_ev; mat(1,:) - mat(2,:); mat(2,:) - mat(3,:); mat(3,:) - mat(4,:)];
    bid_FCR_ev = [bid_FCR_ev; ((t-1)*4+1:t*4)' mat(1:4,:)];
    bid_aFRRp_ev = [bid_aFRRp_ev; ((t-1)*4+1:t*4)' mat(5:8,:)];
    bid_aFRRn_ev = [bid_aFRRn_ev; ((t-1)*4+1:t*4)' mat(9:12,:)];
end
Beq_4hour_bid_ev = zeros(height(Aeq_4hour_bid_ev),1);

%%% Back-up resource #######################################################################################
% ESS = back-up resource??


%%% def matrix A,B,Aeq,Beq
constA_VPP = [A_cbl_ev; A_cbl_ess];
constA_EV = [A_uc_ev; A_p_uc_ev; A_soc_ev_ub; A_soc_ev_lb; A_goalsoc_ev];
constA_ESS = [A_uc_ess; A_p_uc_ess; A_soc_ess_ub; A_soc_ess_lb; A_goalsoc_ess];


constB_VPP = [B_cbl_ev; B_cbl_ess];
constB_EV = [B_uc_ev; B_p_uc_ev; B_soc_ev_ub; B_soc_ev_lb; B_goalsoc_ev];
constB_ESS = [B_uc_ess; B_p_uc_ess; B_soc_ess_ub; B_soc_ess_lb; B_goalsoc_ess];

A = [constA_VPP; constA_ESS; constA_EV]; 
B = [constB_VPP; constB_ESS; constB_EV];
Aeq=[Aeq_4hour_bid_ev]; Beq=[Beq_4hour_bid_ev];

%% test
% A = [constA_EV]; B=[constB_EV];
% Aeq = []; Beq = [];


%% Optimization - objective functions
%%% Objective: energy bid, 1~24
bid_energy = zeros(24,lenN);
% add ev
temp.v = 0;
for vdx = 1:N.ev
    dur = ev.in(vdx):ev.out(vdx);
    mat = eye(ev.duration(vdx));
    dev = [mat.*abs(market.S_FCR(dur)).*market.E_FCR(dur) mat.*abs(market.S_aFRR_p(dur)).*market.E_aFRR_p(dur) mat.*market.S_aFRR_n(dur).*market.E_aFRR_n(dur) zeros(ev.duration(vdx))];
    bid_energy(ev.in(vdx):ev.out(vdx), temp.v+1 : temp.v+ev.duration(vdx)*(N.market+N.sc)) = dev;
    temp.v = temp.v + ev.duration(vdx)*(N.market+N.sc);
end

% add ess
temp.v = N.uc_ev;
for vdx = 1:N.ess
    dur = 1:T;
    mat = eye(ess.duration(vdx));
    dev = [mat.*abs(market.S_FCR(dur)).*market.E_FCR(dur) mat.*abs(market.S_aFRR_p(dur)).*market.E_aFRR_p(dur) mat.*market.S_aFRR_n(dur).*market.E_aFRR_n(dur) zeros(ess.duration(vdx))];
    bid_energy(:,temp.v+1 : temp.v+ess.duration(vdx)*(N.market+N.sc)) = dev;
    temp.v = temp.v + ess.duration(vdx)*(N.market+N.sc);
end

%%% Objective: capacity bid, 1~24
bid_capacity = zeros(24,lenN);
temp.v = 0;
for vdx = 1:N.ev
    dur = ev.in(vdx):ev.out(vdx);
    mat = eye(ev.duration(vdx));
    dev = [mat.*market.C_FCR(dur) mat.*market.C_aFRR_p(dur) mat.*market.C_aFRR_n(dur) zeros(ev.duration(vdx))];
    bid_capacity(ev.in(vdx):ev.out(vdx), temp.v+1 : temp.v+ev.duration(vdx)*(N.market+N.sc)) = dev;
    temp.v = temp.v + ev.duration(vdx)*(N.market+N.sc);
end

% add ess
temp.v = N.uc_ev;
for vdx = 1:N.ess
    dur = 1:T;
    mat = eye(ess.duration(vdx));
    dev = [mat.*market.C_FCR(dur) mat.*market.C_aFRR_p(dur) mat.*market.C_aFRR_n(dur) zeros(ess.duration(vdx))];
    bid_capacity(:,temp.v+1 : temp.v+ess.duration(vdx)*(N.market+N.sc)) = dev;
    temp.v = temp.v + ess.duration(vdx)*(N.market+N.sc);
end

%%% Objective: charging cost, 1~24
cost_charge = zeros(24,lenN);
temp.d = 0; temp.v = 0;
for vdx = 1:N.ev
    dur = ev.in(vdx):ev.out(vdx);
    cost_charge(ev.in(vdx):ev.out(vdx), temp.v+ev.duration(vdx)*N.market+1 : temp.v+ev.duration(vdx)*(N.market+N.sc)) = eye(ev.duration(vdx)) .* market.sc(dur);
    temp.v = temp.v + ev.duration(vdx)*(N.market+N.sc);
end

% add ess
temp.v = N.uc_ev;
for vdx = 1:N.ess
    dur = 1:T;
    cost_charge(:,temp.v+ess.duration(vdx)*N.market+1 : temp.v+ess.duration(vdx)*(N.market+N.sc)) = eye(ess.duration(vdx)).*market.sc(dur);
    temp.v = temp.v + ess.duration(vdx)*(N.market+N.sc);
end

f = -sum(bid_capacity) - sum(bid_energy) + sum(cost_charge); % f: 목적함수
% f = zeros(lenN,1);
% f = -ones(1,lenN);
% f = sum(cost_charge);
x0 = zeros(lenN,1);
x0 = [];
%% Solver: Mixed-integer Linear Programming
options = optimoptions('intlinprog','MaxTime',300);%,'Display','off');
[x,fval,exitflag,output] = intlinprog(f,intcon,A,B,Aeq,Beq,lb,ub,x0,options);
% save('Bid.mat','bid_FCR_ev','bid_aFRRp_ev','bid_aFRRn_ev');

end