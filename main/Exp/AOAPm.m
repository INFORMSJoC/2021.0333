function [PCS,EOC]=AOAPm(k,n0,T,mu0,sigma0,v,num,m)
%function [PCS,EOC] = AOAPm(k,n0,T,sigma,num,m,truemu)
%function [PCS,EOC] = AOAPm(k,n0,T,num,m)

PCS = zeros(num,T);
EOC = zeros(num,T);

tic

parfor t = 1:num
    
    mu = normrnd(mu0.* ones(n0,k),sqrt(sigma0).* ones(n0,k));
    
    
    X0 = normrnd(mu.* ones(n0,k), sqrt(v).* ones(n0,k));
    %X0 = normrnd(truemu .* ones(n0,k), sqrt(sigma) .* ones(n0,k));
    %alpha0 = unifrnd(1,20,1,k);
    %beta0 = unifrnd(1,20,1,k);
    %theta = betarnd(alpha0,beta0);
    %X0 = binornd(ones(n0,k),theta.*ones(n0,k));
    %alpha0 = unifrnd(2,10,1,k);
    %beta0 = unifrnd(1,2,1,k);
    %lambda = gamrnd(alpha0,beta0);
    %truemu = 1./lambda; 
    %X0 = exprnd(repmat(truemu,n0,1),[n0 k]);
        
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

    %pv = pm.*(1-pm)./(alpha+beta+1);
    %pv = alpha.*beta.^2;
    
    
    Nrv = normrnd(0,1,1,T);
    
    [~,rb] = sort(mu,'descend');
    %[~,rb] = sort(mu);
    %[~,rb] = sort(truemu);
    %[~,rb] = sort(theta,'descend');
    %[~,rb] = sort(truemu,'descend');
    rb = sort(rb(1:m));
    
    for budget=1:T
        
        [~,id1] = sort(pm,'descend');
        %[~,id1] = sort(pm);
        %[~,id1] = sort(1./pm,'descend');
        id1 = sort(id1(1:m));
        
        PCS(t,budget) = prod(id1==rb) / num;
        EOC(t,budget) = EOC(t,budget)+(sum(mu(rb))-sum(mu(id1)))/num;
        %EOC(t,budget) = EOC(t,budget)+(sum(truemu(id1))-sum(truemu(rb)))/num;
        %EOC(t,budget) = EOC(t,budget)+(sum(theta(rb))-sum(theta(id1)))/num;
        %EOC(t,budget) = EOC(t,budget)+(sum(truemu(rb))-sum(truemu(id1)))/num;

        cm = setdiff((1:k),id1);
        
        V = zeros(1,k);

        for i = 1:k
            nv = pv;
            M = N;
            M(i) = N(i)+1;
            nv(i) = (1/sigma0(i)+M(i)/estvar(i))^(-1);
            %nv(i) = (M(i)/estvar(i))^(-1);
            %nv(i) = pm(i)*(1-pm(i))/(alpha(i)+beta(i)+2);
            %nv(i) = pm(i)*beta0(i)/(1+beta0(i)*(N(i)*estmean(i)+truemu(i)));
            a = 1:m;
            b = 1:(k-m);
            V1 = (pm(id1(a))'-pm(cm(b))).^2./(nv(id1(a))'+nv(cm(b)));
            %V1  =  (pm(id1(a))'-pm(cm(b))+10^(-5)).^2./(nv(id1(a))'+nv(cm(b)));
            V(i) = min(min(V1));
        end
        [~,id2] = max(V);
        
        mm = estmean(id2);
        x = mu(id2)+(v(id2)).^(1/2).*Nrv(budget);
        %x = truemu(id2) + (sigma(id2)) .^ (1/2) .* Nrv(budget);
        %x = binornd(1,theta(id2));
        %x = exprnd(truemu(id2));
        estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        %estvar(id2) = (N(id2)./(N(id2)+1)).*(estvar(id2)+(mm-x).^2./(N(id2)+1));
        N(id2) = N(id2)+1;
        
        pv(id2) = (1./sigma0(id2)+N(id2)./estvar(id2)).^(-1);
        %pv(id2) = (N(id2)./estvar(id2)).^(-1);
        
        pm(id2) = pv(id2).*(mu0(id2)./sigma0(id2)+N(id2).*estmean(id2)./estvar(id2));
        %pm(id2) = pv(id2).*(N(id2).*estmean(id2)./estvar(id2));
        %alpha(id2) = alpha(id2)+x;
        %beta(id2) = beta(id2)+1-x;
        %pm(id2) = alpha(id2)./(alpha(id2)+beta(id2));
        %alpha(id2) = alpha(id2)+1;
        %beta(id2) = beta0(id2)./(1+beta0(id2).*N(id2).*estmean(id2));
        %pm(id2) = alpha(id2).*beta(id2);

        %pv(id2) = pm(id2).*(1-pm(id2))./(alpha(id2)+beta(id2)+1);
    end
end
toc

PCS = sum(PCS);
EOC = sum(EOC);
end