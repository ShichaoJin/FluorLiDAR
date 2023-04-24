% poissonf(750,10,350)
% 泊松分布函数

function Pnn=poissonf(d,n,NN)
% d是像元中的株数；n是像元中的样地数量；NN是样地中最大允许有的株数，规定为350

Px_tot=0; % 累计概率初值
m=d/n;    % 计算样地中的平均株数m

Pnn(1,1)=exp(-m); % 计算x==0时的概率
Px_tot=Px_tot+Pnn; % 计算累计概率

for ak=1:NN % 循环NN次，即从样地中有1株树开始，直到样地中有NN株
    Pnn(1,ak+1)=Pnn(1,ak)*m/ak; % 利用有ak-1株树时的概率计算ak株树时的概率，注意：Pnn中的序号比实际的株数大1
    Px_tot=Px_tot+Pnn(1,ak+1); % 计算累计概率
end

if Px_tot<0.95 % 如果概率全部相加后小于0.95，则提示错误信息
    msgbox('Possible error in Neyman distribution');
end

for k=1:(NN+1) % 重新划分概率，使之分布在0-1之间
    Pnn(1,k)=Pnn(1,k)/Px_tot;
end