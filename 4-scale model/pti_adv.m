% pti=pti_adv(10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g)
% Subroutine that calculates the gap fraction (Pvg and Pig)

function [pti]=pti_adv(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g) % 输入角度
% b是像元面积,n是像元中的样地数量,norp,cors只能等于1或2,d是像元中的株数,m2是一个群落的平均株数,
% k是样地中最大允许有的株数(规定为350),r筒半径或椭球的水平半径,za太阳天顶角或观测天顶角(角度),
% hb筒高或椭球垂直长,alpha锥顶角的一半,a和c是计算平均叶倾角的参数(其中c等于论文中的b),
% omegae嫩枝聚集度指数,gamae针对枝的面积比例,l叶面积指数,theta_g坡度,azimuth_view观测方位角

sza=acos(cos(theta_g*pi/180)*cos(szap*pi/180)+sin(theta_g*pi/180)*sin(szap*pi/180)*cos(azimuth_sun*pi/180-azimuth_dem*pi/180));
vza=acos(cos(theta_g*pi/180)*cos(vzap*pi/180)+sin(theta_g*pi/180)*sin(vzap*pi/180)*cos(azimuth_view*pi/180-azimuth_dem*pi/180));

% 选择样地内树的某种分布是纽曼分布(1)还是泊松分布(2)
if norp==1
    px=neyman(d,n,m2,k);
end
if norp==2
    px=poissonf(d,n,k);
end

% 选择样地内植被的形状是筒形(1)还是椭球形(2)，计算crown projection on the ground
if cors==1
    [tab_v,tac_v]=cone_taborc(r,vzap,hb,alpha);
    ta_v=tab_v+tac_v;
    [tab_s,tac_s]=cone_taborc(r,szap,hb,alpha);
    ta_s=tab_s+tac_s;
end
if cors==2
    ta_v=spheroid_ta(vzap,hb,r);
    ta_s=spheroid_ta(szap,hb,r);
end

pgapv=pgap_ax_b(a,c,omegae,gamae,l,d,b,vzap,r,hb,alpha,cors); % 计算一株树的孔隙率
pgaps=pgap_ax_b(a,c,omegae,gamae,l,d,b,szap,r,hb,alpha,cors); % 计算一株树的孔隙率
A=b/n/cos(theta_g*pi/180); % 斜面面积
i_m=d/n; % 样地平均株数

for choice=1:2
    if choice==1
        aaa=ta_s; % 整个树的tas
        za=sza;
        P_s_up=0;
        P_s_down=0;
    end
    if choice==2
        aaa=ta_v; % 整个树的tav
        za=vza;
        P_v_up=0;
        P_v_down=0;
    end
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
                        PNN=PNN*(j+v)/v*(1-aa/(A*cos(za)));
                    end
                end
                Ptj_i(i)=PNN*(aa/(A*cos(za)))^j;
            end
            Ptj=Ptj_i(i)+Ptj;
        end
        %往下开始计算Ptj之外的
        if choice==1
            P_s_up=P_s_up+Ptj*(1-pgaps^j);
            P_s_down=P_s_down+Ptj*(1-pgaps)*j;
        end
        if choice==2
            P_v_up=P_v_up+Ptj*(1-pgapv^j);
            P_v_down=P_v_down+Ptj*(1-pgapv)*j;
        end
    end
end

Qs=P_s_up/P_s_down; % 光线相对于斜面和相对于平面一样
Qv=P_v_up/P_v_down; % 视线相对于斜面和相对于平面一样
Qsv=min(Qs,Qv); % 相对于斜面和相对于平面一样

pigneyman=overlap_pvgorpig(b,n,1,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g); % neyman时的pig
pigpoisson=overlap_pvgorpig(b,n,2,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g);% poisson时的pig
OmegaT = log(pigneyman)/log(pigpoisson); % 以e为底
% 为什么是pig而不是pvg，聚集度随着太阳角度改变吗？不随观测角度改变吗？？？？？？？？？？？？？？？？？？？？？？？？？？？？？
Lt=OmegaT*pi*r*r*d/b; % 在垂直方向上的投影
Wt=(pi*r*r)^0.5; % 在垂直方向上的投影
E_r=Wt/Lt;
thetae=atan(2*r/E_r); % 文章中与程序中不同!!用斜面上的角度？？？？？？？？？？？？？？
if thetae<0 % 保证theta不小于0
    thetae=-thetae;
end

phi_r_sv=abs(azimuth_sun-azimuth_view);
if phi_r_sv>180
    phi_r_sv=360-phi_r_sv;
end
F=1-phi_r_sv*pi/180/thetae;%????????????????????????????????????????
if F<0
    F=0;
end
if F>1
    F=1;
end

Qcb=Qv*Qs+F*(Qsv-Qv*Qs);

if cors==1
    [tib,tic]=cone_ticorb(szap,vzap,alpha,phi_r_sv,r,hb);
    ti=tic+tib;
    pti=Qcb*ti/(Qv*ta_v);
end
if cors==2
    ti=spheroid_ti(vzap,szap,hb,r,phi_r_sv);
    pti=Qcb*ti/(Qv*ta_v);
end