function [PCS,EOC] = pSAR(k,n0,T,mu0,sigma0,sigma,num,m)
%function [PCS,EOC] = pSAR(k,n0,T,num,m)

PCS = zeros(num,T);
EOC=zeros(num,T);

parfor_progress(num);

tic

parfor t=1:num
    
    truemu = normrnd(mu0, sqrt(sigma0));
    
    [~,rb]=sort(truemu,'descend');
    %alpha0=unifrnd(1,10,1,k);
    %beta0=unifrnd(1,10,1,k);
    %theta = betarnd(alpha0,beta0);
    %[~,rb]=sort(theta,'descend');
    rb=sort(rb(1:m));
    
    LK = 1/2 + sum(1./[2:1:k]);
    
    for budget=1:T
        
        NK = ceil((n0*k+budget-1-k)/LK./(k+1-[1:1:(k-1)]));
        
        NS = diff(NK);
        A = [1:1:k];
        J=[];
        mm = m;
        
        X0 = normrnd(truemu.* ones(NK(1),k), sqrt(sigma) .* ones(NK(1),k));
        %X0 = binornd(ones(NK(1),k),theta.*ones(NK(1),k));
        
        estmean=mean(X0);
        N = NK(1).*ones(1,k);
        
        [~,id1]=sort(estmean,'descend');
        
        Delta1 = estmean(id1(1:mm))-estmean(id1(mm+1));
        Delta2 = estmean(id1(mm))-estmean(id1(mm+1:k));
        Delta = [Delta1 Delta2];
        
        for i = 1:(k-1)
            
            [~,id4]=max(Delta);
            id4 = A(id1(id4));
            
            if estmean(id4)>estmean(A(id1(mm+1)))
                mm = mm - 1;
                J = [J id4];
            end
            
            if length(J)==m
                break
            else
                
                A = setdiff(A,id4);
                
                X0 = normrnd(truemu(A).* ones(NS(i),length(A)), sqrt(sigma(A)) .* ones(NS(i),length(A)));
                %X0 = binornd(ones(NS(i),length(A)),theta(A).*ones(NS(i),length(A)));

                estmean(A) = (estmean(A).*N(A)+sum(X0))./(N(A)+NS(i));
                N(A) = N(A) + NS(i);
                
                [~,id1]=sort(estmean(A),'descend');
                
                Delta1 = estmean(A(id1(1:mm)))-estmean(A(id1(mm+1)));
                Delta2 = estmean(A(id1(mm)))-estmean(A(id1(mm+1:(k-i))));
                Delta = [Delta1 Delta2];
                
            end
            
        end
        
        J = sort(J);
        
        PCS(t,budget)=prod(J==rb) / num;
        EOC(t,budget)=EOC(t,budget)+(sum(truemu(rb))-sum(truemu(J)))/num;
        
    end
    
    parfor_progress;
    
end
parfor_progress(0);
toc
PCS = sum(PCS);
EOC = sum(EOC);
end