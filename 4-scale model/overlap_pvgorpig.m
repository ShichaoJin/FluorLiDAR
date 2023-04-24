% [pig]=overlap_pvgorpig(10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,szap,azimuth_dem,azimuth_view_or_sun,theta_g)
% [pvg]=overlap_pvgorpig(10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,vzap,azimuth_dem,azimuth_view_or_sun,theta_g)
% Subroutine that calculates the gap fraction (Pvg and Pig)

function [pvgorpig pvg_c]=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,zap,azimuth_dem,azimuth_view_or_sun,theta_g) % 输入角度
% b是像元面积,n是像元中的样地数量,norp,cors只能等于1或2,等于3就是随机分布,d是像元中的株数,m2是一个群落的平均株数,
% k是样地中最大允许有的株数(规定为350),r筒半径或椭球的水平半径,za太阳天顶角或观测天顶角(角度),
% hb筒高或椭球垂直长,alpha锥顶角的一半,a和c是计算平均叶倾角的参数(其中c等于论文中的b),
% omegae嫩枝聚集度指数,gamae针对枝的面积比例,l叶面积指数,theta_g坡度,azimuth_view观测方位角

za=acos(cos(theta_g*pi/180)*cos(zap*pi/180)+sin(theta_g*pi/180)*sin(zap*pi/180)*cos(azimuth_view_or_sun*pi/180-azimuth_dem*pi/180));

% 选择样地内树的某种分布是纽曼分布(1)还是泊松分布(2)
if norp==1
    px=neyman(d,n,m2,k);
end
if norp==2
    px=poissonf(d,n,k);
end

% 选择样地内植被的形状是筒形(1)还是椭球形(2)，计算crown projection on the ground,vg或sg
if cors==1
    [tab,tac]=cone_taborc(r,zap,hb,alpha);
    ta=tab+tac;
end
if cors==2
    ta=spheroid_ta(zap,hb,r);
end

pgap=pgap_ax_b(a,c,omegae,gamae,l,d,b,zap,r,hb,alpha,cors); % 计算一株树的孔隙率
A=b/n*cos(za)/cos(theta_g*pi/180); % 斜面面积
i_m=d/n; % 样地平均株数

% ------------------------------------- 4-scale gap fraction -------------------------------------
% 计算pvg_c,样地内有0株树挡住视线的情况
pvg_c=0;
tat=ta; % 临时变量
for i=1:k % 样地内有0棵树挡住视线时有1-350株树的情况
    if i>i_m % 株数多时缩小ta % 在增加锥筒和椭球投影面积的控制条件后，这里的控制可以使看到地面的概率变的比不用这个控制大一点
        tat=ta*i_m/i;
    end
    if tat>A
        tat=A;
    end
    pvg_c=pvg_c+px(i+1)*(1-tat/A)^i;
end

aaa=ta; % 整个树的Vg or Sg
P=0;
for j=1:k % 循环 Pvg,j==0时只计算pvg_c
    Ptj=0;
    for i=j:k % 循环 Ptj
        Ptj_i(i)=0;
        if i>i_m && Ptj_i(i-1)<0.000000001 && i~=j % 如果i-1棵树中有遮住光线或视线的概率很小，株数大于平均水平，存在Ptj(i-1)，则Ptj(i)不用计算直接赋值
            Ptj_i(i)=0;
        else
            aa=aaa; % 临时变量 % Vg or Sg
            if i>i_m % 在增加锥筒和椭球投影面积的控制条件后，这里的控制可以使看到地面的概率变的比不用这个控制大一点
                aa=aaa*i_m/i; % 如果株数大于平均水平，则sg or vg按比例缩放
            end
            if aaa>A % 在增加锥筒和椭球投影面积的控制条件后，这里的控制不起作用，但没有前面的控制时这个作用很大
                aa=A;  % 如果这个角度下的阴影面积超过了样地的面积，则令样地面积等于阴影面积
            end
            PNN=px(i+1); % 临时变量，样地中有i株数的概率
            if i~=j
                for v=1:(i-j) % 循环计算二项分布
                    PNN=PNN*(j+v)/v*(1-aa/A);
                end
            end
            Ptj_i(i)=PNN*(aa/A)^j;
        end
        Ptj=Ptj_i(i)+Ptj;
    end
    %往下开始计算Ptj之外的
    P=P+Ptj*(pgap^j);
end

pvgorpig=pvg_c+P;
if pvgorpig<0
    pvgorpig=0;
end
if pvgorpig>1
    pvgorpig=1;
end
