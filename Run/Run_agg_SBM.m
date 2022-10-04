clc, clear, close all, clear global
tic
folder = ['f09','\autoever_0531'];

pipeline % csv read
time_reduction % 48 hours => 12 slot
[UC, Bid, VPP, ESS, EV, SoC, exitflag] = CoSchedule(ESS, EV, ev, market, passive, Bid); % MILP scheduling algorithm
% "UC": 성능 검증을 위한 알고리즘 결과
% "Bid": 알고리즘이 도출한 최종 output
% "VPP": TE 출력의 총합
toc
%% excel format
netload = passive.D - passive.PV;
mat_EV = zeros(T*max(EV.class),5);
for i = 1:max(EV.class)
    mat_EV(T*(i-1)+1 : T*i, :) = [sum(UC.EV{i},2) SoC.EV{i}(:,1) SoC.EV{i}(:,2)-SoC.EV{i}(:,1) SoC.EV{i}(:,2) SoC.EV_per{i}];
end

str_EV = strings(T*max(EV.class),1);
list = ["FCR", "aFRRp" "aFRRn" "mFRRp" "mFRRn" "SCpos" "SCneg"];
for i = 1:max(EV.class)
    for t = 1:T
        if mat_EV(T*(i-1)+t,1) ~= 0
            idx = find(mat_EV((T*(i-1)+t),1) == UC.EV{i}(t,:));
            if idx == 6
                if UC.EV{i}(t,idx) >= 0
                    str_EV(T*(i-1)+t) = list(6);
                else
                    str_EV(T*(i-1)+t) = list(7);
                end
            else
                str_EV(T*(i-1)+t) = list(idx);
            end
        end
    end
end

mat_ESS = zeros(T*ESS.num,5);
str_ESS = strings(T*ESS.num,1);
for i = 1:ESS.num
    mat_ESS(T*(i-1)+1 : T*i, :) = [sum(UC.ESS{i},2) SoC.ESS{i}(:,1) SoC.ESS{i}(:,2)-SoC.ESS{i}(:,1) SoC.ESS{i}(:,2) SoC.ESS_per{i}];
end
for i = 1:ESS.num
    for t = 1:T
        if mat_ESS(T*(i-1)+t,1) ~= 0
            idx = find(mat_ESS((T*(i-1)+t),1) == UC.ESS{i}(t,:));
            if idx == 6
                if UC.ESS{i}(t,idx) >= 0
                    str_ESS(T*(i-1)+t) = list(6);
                else
                    str_ESS(T*(i-1)+t) = list(7);
                end
            else
                str_ESS(T*(i-1)+t) = list(idx);
            end
        end
    end
end

strEV_rows = reshape(str_EV,T,max(EV.class));



%% cost analysis for paper
% cost = zeros(1,max(EV.class));
% for i = 1:max(EV.class)
%     cost(i) = sum(sum(market.C{1}.*UC.EV{i}(:,1:5))) - sum(market.SC(:,1).*UC.EV{i}(:,6));    
% end
% 
% cost = [sum(cost) cost];
% save('mFRR.mat','cost')

%% soc check
SoC_box = [];
for i = 1:max(EV.class)
    SoC_box = [SoC_box; SoC.EV_per{i}];
end

save('SoC.mat','SoC_box')