% [PG,Ft]=PGandFt(0.5,10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g)
% Subroutine that calculates the gap fraction (Pvg and Pig)

function [PG,Ft]=PGandFt(ha,b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g) % ����Ƕ�
% b����Ԫ���,n����Ԫ�е���������,norp,corsֻ�ܵ���1��2,d����Ԫ�е�����,m2��һ��Ⱥ���ƽ������,
% k����������������е�����(�涨Ϊ350),rͲ�뾶�������ˮƽ�뾶,za̫���춥�ǻ�۲��춥��(�Ƕ�),
% hbͲ�߻�����ֱ��,alpha׶���ǵ�һ��,a��c�Ǽ���ƽ��Ҷ��ǵĲ���(����c���������е�b),
% omegae��֦�ۼ���ָ��,gamae���֦���������,lҶ���ָ��,theta_g�¶�,azimuth_view�۲ⷽλ��

sza=acos(cos(theta_g*pi/180)*cos(szap*pi/180)+sin(theta_g*pi/180)*sin(szap*pi/180)*cos(azimuth_sun*pi/180-azimuth_dem*pi/180));

if cors==1
    hc=r/tan(alpha*pi/180); % Բ׶��
    H=(hc/3+hb+ha); % �����ֵĸ�֮��Ϊ���ߣ�hcΪʲô/3����������H �˴��Ȳ���cos�仯���ں����
    [tab,tac]=cone_taborc(r,szap,hb,alpha);
    ta=tab+tac;
end
if cors==2
    H=hb+ha; % ��������������hb��������
    ta=spheroid_ta(szap,hb,r);
end

[pigneyman pvg_cn]=overlap_pvgorpig(b,n,1,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g); % neymanʱ��pig
[pigpoisson pvg_cp]=overlap_pvgorpig(b,n,2,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g);% poissonʱ��pig
OmegaT = log(pvg_cn)/log(pvg_cp); % ��eΪ��----------------------------------------------------------------------------------------!!!!!!!!!!??????????
% Ϊʲô��pig������pvg���ۼ�������̫���Ƕȸı��𣿲���۲�Ƕȸı���???????????????????????????????????????????????????????????

Lt=OmegaT*ta*d/(b*cos(sza)/cos(theta_g*pi/180)); % ��б���ϵ�ͶӰ
Wt=ta^0.5; % ��б���ϵ�ͶӰ

H=H*cos(theta_g*pi/180)/cos(sza);

xi=acos(cos(szap*pi/180)*cos(vzap*pi/180)+sin(szap*pi/180)*sin(vzap*pi/180)*cos(azimuth_sun*pi/180-azimuth_view*pi/180));
if xi<pi
    lambda_m=H*tan(xi); % ������Сlambda
end
if xi>=pi
    lambda_m=0;
end

pvg=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,vzap,azimuth_dem,azimuth_view,theta_g);
pig=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g);

Cold=pig*pvg; % ps.c�е�PS
Hot=pig;

if xi<=pi/2
    IN1=0; % ���ӻ���
    IN2=0; % ��ĸ����
    FLAG_IN1=0;
    FLAG_IN2=0;
    MAX=200000; % ���ֵ�������
    INC=0.1; % ���ּ��
    for i=0:MAX % �����
        if (lambda_m+INC*i)==0 % lambda_m=0.
            IN2=1; % ��㸳ֵ
            IN1=0;
            break;
        else % lambda_m~=0
            IN1=IN1+exp(-Lt*(1+(lambda_m+INC*i)/Wt))/atan((lambda_m+INC*i)/H); % �������Ѿ�����ͬʱԼ��Լ����(Lt/Wt)��һ��
            IN2=IN2+exp(-Lt*(1+(lambda_m+INC*i)/Wt));
            if FLAG_IN1==IN1 && FLAG_IN2==IN2 % ������»��ֲ����ˣ�������ѭ��
                i=MAX;
            end
            FLAG_IN1=IN1;
            FLAG_IN2=IN2;
        end
    end
    Ft=0;
    if IN2~=0
        Ft=1-IN1/IN2*xi; % �����ȵ㺯��
    end
    if xi<0.0001
        Ft=1;
    end
    PG=Cold+(Hot-Cold)*Ft; % ����ɼ����յ���
else
    Ft=0;
    PG=Cold;
end

if Ft<0
    Ft=0;
end