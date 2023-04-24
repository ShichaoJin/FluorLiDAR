% [PT,Fs]=PTandFs(0.035,10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g)
% Subroutine that calculates the gap fraction (Pvg and Pig)

function [PT,Fs]=PTandFs(Ws,b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,cp,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g) % 输入角度
% b是像元面积,n是像元中的样地数量,norp,cors只能等于1或2,d是像元中的株数,m2是一个群落的平均株数,
% k是样地中最大允许有的株数(规定为350),r筒半径或椭球的水平半径,za太阳天顶角或观测天顶角(角度),
% hb筒高或椭球垂直长,alpha锥顶角的一半,a和c是计算平均叶倾角的参数(其中c等于论文中的b),
% omegae嫩枝聚集度指数,gamae针对枝的面积比例,l叶面积指数,theta_g坡度,azimuth_view观测方位角



[Q1tot Q2tot]=Q(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,cp,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g);
pti=pti_adv(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g); % 可视光照面比可视树冠

Hot=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g); % Hot==Pig，光照地面面积比样方面积??????????????????????vzap????????szap???????????
ptf=pti*Q1tot+(1-pti)*Q2tot; % 可视光照叶片比样方面积
% 选择样地内植被的形状是筒形(1)还是椭球形(2)
if cors==1
    hc=r/tan(alpha*pi/180); % 计算锥的高度
    V=pi*r^2*hb+1/3*pi*r^2*hc; % 计算树冠的体积
    [tab,tac]=cone_taborc(r,szap,hb,alpha);
    ta=tab+tac;
end
if cors==2
    V=4/3*pi*r^2*(hb*0.5); % 计算椭球的体积
    ta=spheroid_ta(szap,hb,r);
end

xi=acos(cos(szap*pi/180)*cos(vzap*pi/180)+sin(szap*pi/180)*sin(vzap*pi/180)*cos((azimuth_sun-azimuth_view)*pi/180));

Ss=V/ta;% tas
Lo=l*b/(d*ta);% p_gap_ax_b.c中的
% Lo=l*pi*r^2/ta;  % QFeng,20190327
hs=r/(Ss*Lo);   % ls.c中的？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
lambda_ms=tan(xi)*hs; %ls.c中的

LAI_HOT_SPOT=1;% 为什么==1  和Lt一样，类似于覆盖度指数，不是叶面积指数
Ls=(a+c*szap*pi/180)*omegae*LAI_HOT_SPOT; %ls.c中的 shoot尺度不考虑gamae

if xi<=pi/2  %？？？？？？？？？？？？？？？？？？
    IN1=0; % 分子积分
    IN2=0; % 分母积分
    FLAG_IN1=0;
    FLAG_IN2=0;
    MAX=40000; % 积分的最大次数
    INC=0.001; % 积分间隔
    
    for i=0:MAX % 求积分
        if (lambda_ms+INC*i)==0 % lambda_ms==0
            IN2=1; % 随便赋值
            IN1=0;
            break;
        else % lambda_ms~=0
            IN1=IN1+exp(-Ls*(1+(lambda_ms+INC*i)/Ws))/atan((lambda_ms+INC*i)/hs); % 积分中已经上下同时约分约掉了(Ls/Ws)这一项
            IN2=IN2+exp(-Ls*(1+(lambda_ms+INC*i)/Ws));
            if FLAG_IN1==IN1 && FLAG_IN2==IN2
                i=MAX+1;
            end
            FLAG_IN1=IN1;
            FLAG_IN2=IN2;
        end
    end
    Fs=0;
    if IN2~=0
        Fs=1-IN1/IN2*xi; % 计算热点函数
    end
    if xi<=0.000001
        Fs=1;
    end
    PT=ptf+(1-Hot-ptf)*Fs; % 这个函数与.c一样，只是加一减一的位置不同，计算可见光照叶片
%     PT=ptf+(1-ptf)*Fs;
else
    Fs=0;  %     修改？？？？？？？？？？？？？？？？？？？？？
    PT=ptf;
end