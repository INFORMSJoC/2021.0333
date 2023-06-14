function [estmean] = evacuation(Ntpmu,Ntpsigma,Ncpmu,Ncpsigma,Route)

tt = zeros(22,22);
tc = zeros(22,22);

[row,column] = find(Ntpmu~=0);
for i=1:length(row)
    tt (row(i),column(i)) = lognrnd(Ntpmu(row(i),column(i)),Ntpsigma(row(i),column(i)));
    tc (row(i),column(i)) = ceil(lognrnd(Ncpmu(row(i),column(i)),Ncpsigma(row(i),column(i))));
end

flow = zeros(8,1);
time = zeros(8,1);
for i = 1:8
    flow(i,:) = tc(Route(i,1),Route(i,2));
    len = nnz(Route(i,:));
    for j = 1:5
        if j == len
            break
        else
            time(i,:) = time(i,:) + tt(Route(i,j),Route(i,j+1));
            flow(i,:) = min(flow(i,:),tc(Route(i,j),Route(i,j+1)));
        end
    end
end


ICS = [250 350 305 180];
TA = zeros(4,1);
TA(1) = (ICS(1)+sum((time(1:2)-1).*flow(1:2)))/sum(flow(1:2));
TA(2) = (ICS(2)+sum((time(3:4)-1).*flow(3:4)))/sum(flow(3:4));
TA(3) = (ICS(3)+sum((time(5:6)-1).*flow(5:6)))/sum(flow(5:6));
TA(4) = (ICS(4)+sum((time(7:8)-1).*flow(7:8)))/sum(flow(7:8));
[~,PRI] = sort(TA,'descend');

INA = zeros(8,1);
IPA = zeros(132,14);
Term = 1;
for i = 1:4
    for j = 1:2
        for h = 1:4
            for ell = 1:2
                if (PRI(i)==PRI(h)) && (j==ell)
                    ell = ell+1;
                else
                    idx=1;
                    for m = 1: (nnz(Route((2*PRI(i)-1+j-1),:))-1)
                        for n = 1:(nnz(Route((2*PRI(h)-1+ell-1),:))-1)
                            if Route((2*PRI(i)-1+j-1),m) == Route((2*PRI(h)-1+ell-1),n) && Route((2*PRI(i)-1+j-1),(m+1)) == Route((2*PRI(h)-1+ell-1),(n+1))
                                IPA(Term,(idx:idx+1)) = [Route((2*PRI(i)-1+j-1),m),Route((2*PRI(i)-1+j-1),(m+1))];
                                idx = idx+2;
                            end
                        end
                    end
                    if sum(IPA(Term,:))~=0
                        Term = Term+1;
                        INA(2*(i-1)+j) = INA(2*(i-1)+j)+1;
                    end
                end
            end
        end
    end
end

IPA(((sum(INA)+1):end),:) = [];
IPA(:,find(sum(IPA)==0)) = [];

Phi = zeros(size(IPA,1),size(IPA,2));
numPhi = zeros(8,1);
kk = unique(IPA((1:INA(1)),:),'rows');
idxx = size(kk,1);
numPhi(1) = idxx;
Phi((1:idxx),:) = kk;
for i=2:8
    st = sum(INA(1:(i-1)))+1;
    en = sum(INA(1:i));
    kk = unique(IPA((st:en),:),'rows');
    Phi(((idxx+1):(idxx+size(kk,1))),:) = kk;
    numPhi(i) = size(kk,1);
    idxx = idxx + size(kk,1);
end

Phi(((sum(numPhi)+1):end),:) = [];

flowchongfu = zeros(sum(numPhi),1);
timechongfu = zeros(sum(numPhi),1);
for i = 1:sum(numPhi)
    flowchongfu(i,:) = tc(Phi(i,1),Phi(i,2));
    lenchongfu=nnz(Phi(i,:));
    for j = 1:14
        if j == lenchongfu
            break
        else
            timechongfu(i,:) = timechongfu(i,:) + tt(Phi(i,j),Phi(i,j+1));
            if tc(Phi(i,j),Phi(i,j+1)) == 0
            else
                flowchongfu(i,:) = min(flowchongfu(i,:),tc(Phi(i,j),Phi(i,j+1)));
            end
        end
    end
end
timechongfu = ceil(timechongfu);

lambda = ceil(sum(TA));
timetime = ceil(time);
psiflow = zeros(8,(lambda+1));
gamma = flowchongfu.*ones(sum(numPhi),(lambda+1));
for i = 1:4
    for t = timetime(2*PRI(i)-1):lambda
        for j = 1:2
            if INA(2*(i-1)+j)>0
                if t >=timetime(2*PRI(i)-1+j-1)
                    psiflow(2*(i-1)+j,(t+1)) = flow(2*PRI(i)-1+j-1);
                    for r = 1:numPhi(2*(i-1)+j)
                        tplus = t-timechongfu(sum(numPhi(1:(2*(i-1)+j-1)))+r);
                        psiflow(2*(i-1)+j,(t+1)) = min(psiflow(2*(i-1)+j,(t+1)),gamma((sum(numPhi(1:(2*(i-1)+j-1)))+r),(tplus+1)));
                    end
                    if sum(sum(psiflow(((2*(i-1)+1):(2*i)),(1:(t+1))))) >= ICS(PRI(i))
                        rou = ICS(PRI(i)) - sum(sum(psiflow(((2*(i-1)+1):(2*i)),(1:t))));
                        psiflow(2*(i-1)+j,(t+1)) = min(psiflow(2*(i-1)+j,(t+1)),rou-sum(psiflow(((2*(i-1)+1):(2*(i-1)+j)-1),(t+1))));
                        for r = 1:numPhi(2*(i-1)+j)
                            tplus = t-timechongfu(sum(numPhi(1:(2*(i-1)+j-1)))+r);
                            gamma((sum(numPhi(1:(2*(i-1)+j-1)))+r),(tplus+1)) = gamma((sum(numPhi(1:(2*(i-1)+j-1)))+r),(tplus+1))-psiflow(2*(i-1)+j,(t+1));
                        end
                        break
                    else
                        for r = 1:numPhi(2*(i-1)+j)
                            tplus = t-timechongfu(sum(numPhi(1:(2*(i-1)+j-1)))+r);
                            gamma((sum(numPhi(1:(2*(i-1)+j-1)))+r),(tplus+1)) = gamma((sum(numPhi(1:(2*(i-1)+j-1)))+r),(tplus+1))-psiflow(2*(i-1)+j,(t+1));
                        end
                    end
                else
                    psiflow(2*(i-1)+j,(t+1)) = 0;
                end
            else
                if t>=timetime(2*PRI(i)-1+j-1)
                    psiflow(2*(i-1)+j,(t+1)) = flow(2*PRI(i)-1+j-1);
                    if sum(sum(psiflow(((2*(i-1)+1):(2*i)),(1:(t+1))))) >= ICS(PRI(i))
                        break
                    end
                else
                    psiflow(2*(i-1)+j,(t+1)) = 0;
                end
            end
        end
    end
end

psiflow(:,find(sum(psiflow)==0)) = [];
Tstar = size(psiflow,2)-1;

estmean = Tstar;

% idxxx = find(psiflow(:,(Tstar+1))~=0);
% gap = zeros(size(idxxx,2),1);
% for i=1:size(idxxx,1)
%     if rem(idxxx(i),3)==0
%         gap(i) = 3*PRI(fix(idxxx(i)/3));
%     else
%         gap(i) = 3*(PRI(fix(idxxx(i)/3)+1)-1) + rem(idxxx(i),3);
%     end
% end
% gaptime = timetime(gap)-time(gap);
% estmean = Tstar-min(gaptime);

end