function [R] = AOAPm(T,num)

R = zeros(T, 4, num);

parpool('local',24);

parfor_progress(num);

tic

parfor t = 1:num
    
    truemu = [10, 8.5, 5, 4];
    sigma = [3, 1.1, 3.5, 3.8];
    
    X0 = normrnd(truemu .* ones(10,4), sqrt(sigma) .* ones(10,4));
    
    estmean = mean(X0);
    estvar = var(X0);
    
    N = 10 * ones(1,4);
    
    pv = (N ./ estvar) .^ (-1);
    pm = pv .* (N .* estmean ./ estvar);
    
    Nrv = normrnd(0, 1, 1, T);
    
    for budget=1:T
        
        [~,id1] = sort(pm,'descend');
        id1 = sort(id1(1:2));
        
        R(budget, :, t) = N ./ sum(N);
        
        cm=setdiff((1:4),id1);
        
        V = zeros(1,4);
        
        for i=1:4
            nv=pv;
            M=N;
            M(i)=N(i)+1;
            nv(i)=(M(i)/estvar(i))^(-1);
            a=1:2;
            b=1:2;
            V1 = (pm(id1(a))'-pm(cm(b))).^2./(nv(id1(a))'+nv(cm(b)));
            V(i) = min(min(V1));
        end
        [~,id2]=max(V);
        
        mm = estmean(id2);
        x = truemu(id2) + (sigma(id2)) .^ (1/2) .* Nrv(budget);
        estmean(id2) = (estmean(id2) .* N(id2) + x) ./ (N(id2) + 1);
        estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
        N(id2) = N(id2) + 1;
        
        pv(id2)=(N(id2)./estvar(id2)).^(-1);
        pm(id2)= pv(id2).*(N(id2).*estmean(id2)./estvar(id2));
    end
    parfor_progress;
end
parfor_progress(0);

poolobj = gcp('nocreate');
delete(poolobj);

toc
R = sum(R, 3)./num;
R = R';
end