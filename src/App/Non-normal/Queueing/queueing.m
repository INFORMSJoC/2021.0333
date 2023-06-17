function [estmean] = queueing(x)

len = 100;

Serv1Num = x(1);
Serv2Num = x(2);

Serv1Tim = unifrnd(1,39,1,len);
Serv2Tim = unifrnd(5,45,1,len);

Serv1queue = zeros(2,1);
Serv2queue = zeros(2,1);

Wait1Time = zeros(2,1);
Wait2Time = zeros(2,1);

arrTim = exprnd(5,1,len);

AR = zeros(1,len);

DE = zeros(1,len);

nArrivals = 1;
nExits = 0;

AR(1) = arrTim(1);

Events = zeros(3,1);
Events(1,:) = AR(1);
Events(2,:) = 1;

[~,nextevent] = min(Events(1,:));
Time = Events(1,nextevent);
Customer = Events(2,nextevent);
Type = Events(3,nextevent);

while (nExits~=len)
    
    if Type == 0
        if (Serv1Num~=0)
            sTime = Serv1Tim(Customer);
            Events = [Events [(Time + sTime);Customer;1]];
            Serv1Num = Serv1Num - 1;
        else
            Serv1queue = [Serv1queue [Customer;Time]];
        end
        
        if nArrivals~=len
            nArrivals = nArrivals +1;
            AR(nArrivals) = Time + arrTim(nArrivals);
            Events = [Events [AR(nArrivals);nArrivals;0]];
        end
        
    elseif Type == 1
        if (length(Serv1queue(1,:)) > 1)
            Wait1Time = [Wait1Time [Time - Serv1queue(2,2);Serv1queue(1,2)]];
            sTime = Serv1Tim(Serv1queue(1,2));
            Events = [Events [(Time + sTime);Serv1queue(1,2);1]];
            Serv1queue(:,2) = [];
        else
            Serv1Num = Serv1Num + 1;
        end
        
        if (Serv2Num~=0)
            sTime = Serv2Tim(Customer);
            Events = [Events [(Time + sTime);Customer;2]];
            Serv2Num = Serv2Num - 1;
        else
            Serv2queue = [Serv2queue [Customer;Time]];
        end
        
    elseif Type == 2
        nExits = nExits + 1;
        DE(Customer) = Time;
        
        if (length(Serv2queue(1,:)) > 1)
            Wait2Time = [Wait2Time [Time - Serv2queue(2,2);Serv2queue(1,2)]];
            sTime = Serv2Tim(Serv2queue(1,2));
            Events = [Events [(Time + sTime);Serv2queue(1,2);2]];
            Serv2queue(:,2) = [];
        else
            Serv2Num = Serv2Num + 1;
        end
    end
    
    Events(:,nextevent) = [];
    [~,nextevent] = min(Events(1,:));
    Time = Events(1,nextevent);
    Customer = Events(2,nextevent);
    Type = Events(3,nextevent);
end

estmean = sum(DE - AR)./len;

end