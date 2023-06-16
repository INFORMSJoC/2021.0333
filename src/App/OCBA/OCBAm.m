function [PCS,EOC] = OCBAm(k,n0,T,v,num,m,Delta,mu)

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
    
    ellell  =  T/Delta;
    
    for ell = 1:ellell
        
        gujisigma = sqrt(estvar./N);
        [pmm,id4] = sort(estmean,'descend');
        c = (gujisigma(id4(m+1))*pmm(m)+gujisigma(id4(m))*pmm(m+1))/(gujisigma(id4(m))+gujisigma(id4(m+1)));
        
        w  =  zeros(1,k);
        delta  =  estmean-c;
        Omega = 1:k;
        AW1 = (estvar(Omega).*delta(k).^2)./(estvar(k).*delta(Omega).^2);
        w(k) = 1/sum(AW1);
        w(Omega) = AW1*w(k);
        
        NN  =  Delta*w;
        NN = round(NN);
        [~,id10] = max(NN);
        if sum(NN)>Delta
            NN(id10) = NN(id10)-(sum(NN)-Delta);
        end
        [~,id11] = min(NN);
        if sum(NN)<Delta
            NN(id11) = NN(id11)+(Delta-sum(NN));
        end
        
        [A,B] = sort(NN);
        NN1 = diff(A);
        C(1) = A(1)*k;
        C(2:k) = [(k-1):-1:1].*NN1(1:(k-1));
        D(1) = 0;D(2) = C(1);
        for i = 3:(k+1)
            D(i) = sum(C(1:(i-1)));
        end
        
        for budget  =  1:Delta
            
            [~,id1] = sort(estmean,'descend');
            id1 = sort(id1(1:m));
            
            if id1 == rb
                PCS((ell-1)*Delta+budget) = PCS((ell-1)*Delta+budget)+1/num;
            end
            EOC((ell-1)*Delta+budget) = EOC((ell-1)*Delta+budget)+(sum(mu(rb))-sum(mu(id1)))/num;
            
            for i = 1:(k+1)
                if budget>D(i) && budget < = D(i+1)
                    qu = i;
                    break;
                end
            end
            
            if qu> = 2
                smm = setdiff((1:k),B(1:(qu-1)));
            else
                smm = [1:1:k];
            end
            
            budget1 = budget-sum(C(1:(qu-1)));
            
            q = fix(budget1/length(smm));
            res = budget1-length(smm)*q;
            if res>0
                id2 = smm(res);
            else
                id2 = smm(length(smm));
            end
            
            mm = estmean(id2);
            x = mu(id2)+(v(id2)).^(1/2).*Nrv(budget);
            estmean(id2) = (estmean(id2).*N(id2)+x)./(N(id2)+1);
            estvar(id2) = ((N(id2)-1).*estvar(id2)+(x-mm).*(x-estmean(id2)))./N(id2);
            %estvar(id2) = (N(id2)./(N(id2)+1)).*(estvar(id2)+(mm-x).^2./(N(id2)+1));
            N(id2) = N(id2)+1;
        end
    end
end
toc
end