X = zeros(10,2);
X(:,1) = 11:1:20;
X(:,2) = 20:-1:11;

num = 1000000;
QTim = zeros(num,10);

for i = 1 : 10
    for j = 1 : num
        QTim(j,i) = queueing(X(i,:));
    end
end

truemu = mean(QTim);
sigma = var(QTim);