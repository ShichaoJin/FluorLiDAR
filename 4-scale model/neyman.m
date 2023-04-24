% neyman(3000,40,20,350);
% 纽曼A分布
function Px=neyman(d,n,m2,NN)
% d是像元中的株数；n是像元中的样地数量；m2是一个群落的平均株数；NN是样地中最大允许有的株数，规定为350

m1=d/(n*m2); % 计算m1，一个样地平均有几个群落

Px(1,1)=exp(-m1*(1-exp(-m2))); % 样地里有0棵树的概率,株数比Px数组中的序号少1
Px_tot=Px(1,1); % 累计概率
for k=2:(NN+1) % 样地中有1棵树到NN棵树时的概率
    Px(1,k)=0; % 赋初始值
    if (k-1)>d/n && Px(1,k-1)<0.00001 % 为加快速度,在株数超过样地的平均值且概率小于0.00001的不计算
		Px(1,k)=0; 
    else
        for t=1:(k-1)
            PXX=m1*m2*exp(-m2)/(k-1)*Px(1,k-t);
            for i=2:t
                PXX=PXX*m2/(i-1);
            end
            Px(1,k)=Px(1,k)+PXX;
        end
    end
    Px_tot=Px_tot+Px(1,k); %  making sure that the sum of Px is unity
end

if Px_tot<0.98 % 如果概率全部相加后小于0.98，则提示错误信息
    msgbox('Possible error in Neyman distribution');
end

for k=1:(NN+1) % 重新划分概率，使之分布在0-1之间
    Px(1,k)=Px(1,k)/Px_tot;
end