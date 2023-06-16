num = 1000000;

Ncp = zeros(22,22);
Ntp = zeros(22,22);

Ncp(1,5)=50;Ntp(1,5)=2;Ncp(1,6)=30;Ntp(1,6)=2;Ncp(1,7)=20;Ntp(1,7)=1;
Ncp(5,9)=60;Ntp(5,9)=3;Ncp(6,2)=20;Ntp(6,2)=3;Ncp(6,7)=10;Ntp(6,7)=2;
Ncp(7,12)=40;Ntp(7,12)=2;Ncp(8,12)=20;Ntp(8,12)=4;Ncp(8,11)=30;Ntp(8,11)=5;
Ncp(2,10)=60;Ntp(2,10)=2;Ncp(10,9)=10;Ntp(10,9)=2;Ncp(12,20)=50;Ntp(12,20)=3;
Ncp(10,11)=30;Ntp(10,11)=3;Ncp(9,16)=20;Ntp(9,16)=4;Ncp(16,22)=50;Ntp(16,22)=3;
Ncp(10,4)=40;Ntp(10,4)=5;Ncp(4,17)=50;Ntp(4,17)=3;Ncp(11,3)=20;Ntp(11,3)=1;
Ncp(3,4)=20;Ntp(3,4)=6;Ncp(3,14)=50;Ntp(3,14)=2;Ncp(3,13)=10;Ntp(3,13)=4;
Ncp(13,12)=60;Ntp(13,12)=3;Ncp(13,15)=20;Ntp(13,15)=3;Ncp(14,15)=60;Ntp(14,15)=3;
Ncp(4,14)=30;Ntp(4,14)=2;Ncp(17,18)=20;Ntp(17,18)=2;Ncp(14,18)=20;Ntp(14,18)=1;
Ncp(18,19)=40;Ntp(18,19)=2;Ncp(19,21)=60;Ntp(19,21)=3;Ncp(17,16)=10;Ntp(17,16)=3;
Ncp(2,8)=40;Ntp(2,8)=2;Ncp(11,12)=20;Ntp(11,12)=4;Ncp(15,19)=20;Ntp(15,19)=3;

Ncpmu = zeros(22,22);
Ncpsigma = zeros(22,22);
Ntpmu = zeros(22,22);
Ntpsigma = zeros(22,22);

[row,column] = find(Ncp~=0);
for i=1:length(row)
    Ncpmu(row(i),column(i)) = log((Ncp(row(i),column(i))^2)./sqrt(1+Ncp(row(i),column(i))^2));
    Ntpmu(row(i),column(i)) = log((Ntp(row(i),column(i))^2)./sqrt(0.01+Ntp(row(i),column(i))^2));
    Ncpsigma(row(i),column(i)) = sqrt(log(1/(Ncp(row(i),column(i))^2)+1));
    Ntpsigma(row(i),column(i)) = sqrt(log(0.01/(Ntp(row(i),column(i))^2)+1));
end

PlanN1 = [1 7 12 20 0; 1 6 7 12 20; 1 5 9 16 22];
PlanN2 = [2 8 12 20 0;2 10 9 16 22;2 10 11 12 20];
PlanN3 = [3 14 18 19 21;3 13 12 20 0;3 14 15 19 21];
PlanN4 = [4 14 18 19 21;4 17 16 22 0;4 17 18 19 21];

TotalPlan = zeros(81,8);
Total = 1;
a = nchoosek([1,2,3],2);
b = nchoosek([1,2,3],2);
c = nchoosek([1,2,3],2);
d = nchoosek([1,2,3],2);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for ell = 1:3
                TotalPlan(Total,:) = [a(i,:) b(j,:) c(k,:) d(ell,:)];
                Total = Total +1;
            end
        end
    end
end

for i = 1:81
    Route = zeros(8,5);
    Route((1:2),:) = PlanN1(TotalPlan(i,(1:2)),:);
    Route((3:4),:) = PlanN2(TotalPlan(i,(3:4)),:);
    Route((5:6),:) = PlanN3(TotalPlan(i,(5:6)),:);
    Route((7:8),:) = PlanN4(TotalPlan(i,(7:8)),:);
    for j = 1:num
        Q(j,i) = evacuation(Ntpmu,Ntpsigma,Ncpmu,Ncpsigma,Route);
    end
end

truemu = mean(Q);

[mu,rb] = sort(truemu);