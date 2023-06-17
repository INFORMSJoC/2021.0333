function [PCS,EOC]=OCBAmjia(k,n0,T,mu0,sigma0,v,num,m)
%function [PCS,EOC] = OCBAmjia(k,n0,T,sigma,num,m,truemu)
%function [PCS,EOC] = OCBAmjia(k,n0,T,num,m)

PCS = zeros(num,T);
EOC = zeros(num,T);

tic

parfor t=1:num
    
    mu = normrnd(mu0.* ones(n0,k),sqrt(sigma0).* ones(n0,k));
    
    X0 = normrnd(mu.* ones(n0,k), sqrt(v).* ones(n0,k));
    %X0 = normrnd(truemu .* ones(n0,k), sqrt(sigma) .* ones(n0,k));
    %alpha0 = unifrnd(1,20,1,k);
    %beta0 = unifrnd(1,20,1,k);
    %theta  =  betarnd(alpha0,beta0);
    %X0 = binornd(ones(n0,k),theta.*ones(n0,k));
    %alpha0 = unifrnd(2,10,1,k);
    %beta0 = unifrnd(1,2,1,k);
    %lambda = gamrnd(alpha0,beta0);
    %truemu = 1./lambda;
    
    X0 = exprnd(repmat(truemu,n0,1),[n0 k]);
    
    estmean = mean(X0);
    estvar = var(X0);
    %estvar = var(X0,1);
    
    N = n0*ones(1,k);
    pv = (1./sigma0+N./estvar).^(-1);
    %pv = (N ./ estvar) .^ (-1);
    pm = pv.*(mu0./sigma0+N.*estmean./estvar);
    %pm = pv .* (N .* estmean ./ estvar);
    %alpha = alpha0+sum(X0);
    %beta = beta0+n0-sum(X0);
    %pm = alpha./(alpha+beta);
    %alpha = alpha0+n0;
    %beta = beta0./(1+beta0.*n0.*estmean);
    %pm = alpha.*beta;
    
    Nrv = normrnd(0,1,1,T);
    
    [~,rb] = sort(mu,'descend');
    %[~,rb]=sort(mu);
    %[~,rb] = sort(truemu);
    %[~,rb] = sort(theta,'descend');
    %[~,rb] = sort(truemu,'descend');
    rb = sort(rb(1:m));
    wt = ones(1,k)/k;
    
    for budget = 1:T
        [~,id1] = sort(pm,'descend');
        %[~,id1]=sort(estmean,'descend');
        %[~,id1]=sort(estmean);
        %[~,id1]=sort(pm);
        %[~,id1] = sort(1./pm,'descend');
        id1 = sort(id1(1:m));

        PCS(t,budget) = prod(id1==rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(mu(rb))-sum(mu(id1)))/num;
        %EOC(t,budget) = EOC(t,budget)+(sum(truemu(id1))-sum(truemu(rb)))/num;
        %EOC(t,budget) = EOC(t,budget)+(sum(theta(rb))-sum(theta(id1)))/num;
        %EOC(t,budget) = EOC(t,budget)+(sum(truemu(rb))-sum(truemu(id1)))/num;
        
        [~,s] = min(estmean,[],2);
        [B,id4] = sort(estmean,'descend');
        %[B,id4] = sort(estmean);
        
        a = min(B(m)-B(m-1),B(m+1)-B(m));
        b = min(B(m+1)-B(m),B(m+2)-B(m+1));
        if a >= b
            w1 = zeros(1,k);
            Omega1 = setdiff((1:k),id4(m));
            alpha1 = zeros(1,k);
            alpha1(Omega1) = (estmean(Omega1)-B(m)).^2;
            AW11 = (estvar(Omega1)*alpha1(s))./(estvar(s)*alpha1(Omega1));
            AW21 = (estvar(id4(m))./estvar(Omega1)).*AW11.^2;
            w1(s) = 1/((sum(AW21))^(1/2)+sum(AW11));
            w1(Omega1) = AW11*w1(s);
            w1(id4(m)) = w1(s)*(sum(AW21))^(1/2);
            w = w1;
        else
            w2 = zeros(1,k);
            Omega2 = setdiff((1:k),id4(m+1));
            alpha2 = zeros(1,k);
            alpha2(Omega2) = (estmean(Omega2)-B(m+1)).^2;
            AW12 = (estvar(Omega2)*alpha2(s))./(estvar(s)*alpha2(Omega2));
            AW22 = (estvar(id4(m+1))./estvar(Omega2)).*AW12.^2;
            w2(s) = 1/((sum(AW22))^(1/2)+sum(AW12));
            w2(Omega2) = AW12*w2(s);
            w2(id4(m+1)) = w2(s)*(sum(AW22))^(1/2);
            w = w2;
        end
        
        df = (k*n0+budget)*w-(k*n0+budget-1)*wt;
        % qujian = zeros(1,k);
        %         for i=1:k
        %             qujian (i) = sum (w(1:i));
        %         end
        %         suijishu = rand(1);
        %         id2 = length(qujian(suijishu > qujian)) + 1;
        
        [~,id2] = max(df);

        mm = estmean(id2);
        x = mu(id2)+(v(id2)).^(1/2).*Nrv(budget);
        %x = truemu(id2) + (sigma(id2)) .^ (1/2) .* Nrv(budget);
        %x = binornd(1,theta(id2));
        %x = exprnd(truemu(id2));
        estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        %estvar(id2) = (N(id2)./(N(id2)+1)).*(estvar(id2)+(mm-x).^2./(N(id2)+1));
        N(id2) = N(id2)+1;
        wt = N./sum(N);
        
        pv(id2)=(1./sigma0(id2)+N(id2)./estvar(id2)).^(-1);
        %pv(id2) = (N(id2)./estvar(id2)).^(-1);
        pm(id2)=pv(id2).*(mu0(id2)./sigma0(id2)+N(id2).*estmean(id2)./estvar(id2));
        %pm(id2) = pv(id2).*(N(id2).*estmean(id2)./estvar(id2));
        %alpha(id2) = alpha(id2)+x;
        %beta(id2) = beta(id2)+1-x;
        %pm(id2) = alpha(id2)/(alpha(id2)+beta(id2));
        %alpha(id2) = alpha(id2)+1;
        %beta(id2) = beta0(id2)./(1+beta0(id2).*N(id2).*estmean(id2));
        %pm(id2) = alpha(id2)*beta(id2);
    end
end
toc

PCS = sum(PCS);
EOC = sum(EOC);
end