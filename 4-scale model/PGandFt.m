% [PG,Ft]=PGandFt(0.5,10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g)
% Subroutine that calculates the gap fraction (Pvg and Pig)

function [PG,Ft]=PGandFt(ha,b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g) % 输入角度
% b是像元面积,n是像元中的样地数量,norp,cors只能等于1或2,d是像元中的株数,m2是一个群落的平均株数,
% k是样地中最大允许有的株数(规定为350),r筒半径或椭球的水平半径,za太阳天顶角或观测天顶角(角度),
% hb筒高或椭球垂直长,alpha锥顶角的一半,a和c是计算平均叶倾角的参数(其中c等于论文中的b),
% omegae嫩枝聚集度指数,gamae针对枝的面积比例,l叶面积指数,theta_g坡度,azimuth_view观测方位角

sza=acos(cos(theta_g*pi/180)*cos(szap*pi/180)+sin(theta_g*pi/180)*sin(szap*pi/180)*cos(azimuth_sun*pi/180-azimuth_dem*pi/180));

if cors==1
    hc=r/tan(alpha*pi/180); % 圆锥高
    H=(hc/3+hb+ha); % 三部分的高之和为树高，hc为什么/3？？？？？H 此处先不随cos变化，在后面变
    [tab,tac]=cone_taborc(r,szap,hb,alpha);
    ta=tab+tac;
end
if cors==2
    H=hb+ha; % 对于椭球形树冠hb就是树高
    ta=spheroid_ta(szap,hb,r);
end

[pigneyman pvg_cn]=overlap_pvgorpig(b,n,1,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g); % neyman时的pig
[pigpoisson pvg_cp]=overlap_pvgorpig(b,n,2,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g);% poisson时的pig
OmegaT = log(pvg_cn)/log(pvg_cp); % 以e为底----------------------------------------------------------------------------------------!!!!!!!!!!??????????
% 为什么是pig而不是pvg，聚集度随着太阳角度改变吗？不随观测角度改变吗???????????????????????????????????????????????????????????

Lt=OmegaT*ta*d/(b*cos(sza)/cos(theta_g*pi/180)); % 在斜面上的投影
Wt=ta^0.5; % 在斜面上的投影

H=H*cos(theta_g*pi/180)/cos(sza);

xi=acos(cos(szap*pi/180)*cos(vzap*pi/180)+sin(szap*pi/180)*sin(vzap*pi/180)*cos(azimuth_sun*pi/180-azimuth_view*pi/180));
if xi<pi
    lambda_m=H*tan(xi); % 计算最小lambda
end
if xi>=pi
    lambda_m=0;
end

pvg=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,vzap,azimuth_dem,azimuth_view,theta_g);
pig=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g);

Cold=pig*pvg; % ps.c中的PS
Hot=pig;

if xi<=pi/2
    IN1=0; % 分子积分
    IN2=0; % 分母积分
    FLAG_IN1=0;
    FLAG_IN2=0;
    MAX=200000; % 积分的最大次数
    INC=0.1; % 积分间隔
    for i=0:MAX % 求积分
        if (lambda_m+INC*i)==0 % lambda_m=0.
            IN2=1; % 随便赋值
            IN1=0;
            break;
        else % lambda_m~=0
            IN1=IN1+exp(-Lt*(1+(lambda_m+INC*i)/Wt))/atan((lambda_m+INC*i)/H); % 积分中已经上下同时约分约掉了(Lt/Wt)这一项
            IN2=IN2+exp(-Lt*(1+(lambda_m+INC*i)/Wt));
            if FLAG_IN1==IN1 && FLAG_IN2==IN2 % 如果上下积分不变了，就跳出循环
                i=MAX;
            end
            FLAG_IN1=IN1;
            FLAG_IN2=IN2;
        end
    end
    Ft=0;
    if IN2~=0
        Ft=1-IN1/IN2*xi; % 计算热点函数
    end
    if xi<0.0001
        Ft=1;
    end
    PG=Cold+(Hot-Cold)*Ft; % 计算可见光照地面
else
    Ft=0;
    PG=Cold;
end

if Ft<0
    Ft=0;
end