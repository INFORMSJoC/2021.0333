function [PCS,EOC] = OCBAms(k,n0,T,v,num,m,mu)

PCS = zeros(1,T);
EOC = zeros(1,T);
X0 = zeros(n0,k);

tic

for t = 1:num
    
    for i = 1:n0
        X0(i,:) = normrnd(mu,v.^(1/2));
    end
    
    estmean = mean(X0);
    estvar = var(X0);
    %estvar = var(X0,1);
    
    N = n0*ones(1,k);
    Nrv = normrnd(0,1,1,T);
    [~,rb] = sort(mu,'descend');
    rb = sort(rb(1:m));
    wt = ones(1,k)/k;
    
    for budget = 1:T
        [~,id1] = sort(estmean,'descend');
        id1 = sort(id1(1:m));
        
        if id1==rb
            PCS(budget) = PCS(budget)+1/num;
        end
        EOC(budget) = EOC(budget)+(sum(mu(rb))-sum(mu(id1)))/num;
        
        gujisigma = sqrt(estvar./N);
        [pmm,id4] = sort(estmean,'descend');
        c = (gujisigma(id4(m+1))*pmm(m)+gujisigma(id4(m))*pmm(m+1))/(gujisigma(id4(m))+gujisigma(id4(m+1)));
        
        w  =  zeros(1,k);
        delta  =  estmean-c;
        Omega = 1:k;
        AW1 = (estvar(Omega).*delta(k).^2)./(estvar(k).*delta(Omega).^2);
        w(k) = 1/sum(AW1);
        w(Omega) = AW1*w(k);
        
        df = (k*n0+budget)*w-(k*n0+budget-1)*wt;
        
        [~,id2] = max(df);
        
        mm = estmean(id2);
        x = mu(id2)+(v(id2)).^(1/2).*Nrv(budget);
        estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        %estvar(id2) = (N(id2)./(N(id2)+1)).*(estvar(id2)+(mm-x).^2./(N(id2)+1));
        N(id2) = N(id2)+1;
        wt = N./sum(N);
    end
end
toc
end