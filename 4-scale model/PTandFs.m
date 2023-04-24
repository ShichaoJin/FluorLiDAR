% [PT,Fs]=PTandFs(0.035,10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g)
% Subroutine that calculates the gap fraction (Pvg and Pig)

function [PT,Fs]=PTandFs(Ws,b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,cp,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g) % ����Ƕ�
% b����Ԫ���,n����Ԫ�е���������,norp,corsֻ�ܵ���1��2,d����Ԫ�е�����,m2��һ��Ⱥ���ƽ������,
% k����������������е�����(�涨Ϊ350),rͲ�뾶�������ˮƽ�뾶,za̫���춥�ǻ�۲��춥��(�Ƕ�),
% hbͲ�߻�����ֱ��,alpha׶���ǵ�һ��,a��c�Ǽ���ƽ��Ҷ��ǵĲ���(����c���������е�b),
% omegae��֦�ۼ���ָ��,gamae���֦���������,lҶ���ָ��,theta_g�¶�,azimuth_view�۲ⷽλ��



[Q1tot Q2tot]=Q(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,cp,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g);
pti=pti_adv(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g); % ���ӹ�����ȿ�������

Hot=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g); % Hot==Pig�����յ���������������??????????????????????vzap????????szap???????????
ptf=pti*Q1tot+(1-pti)*Q2tot; % ���ӹ���ҶƬ���������
% ѡ��������ֲ������״��Ͳ��(1)����������(2)
if cors==1
    hc=r/tan(alpha*pi/180); % ����׶�ĸ߶�
    V=pi*r^2*hb+1/3*pi*r^2*hc; % �������ڵ����
    [tab,tac]=cone_taborc(r,szap,hb,alpha);
    ta=tab+tac;
end
if cors==2
    V=4/3*pi*r^2*(hb*0.5); % ������������
    ta=spheroid_ta(szap,hb,r);
end

xi=acos(cos(szap*pi/180)*cos(vzap*pi/180)+sin(szap*pi/180)*sin(vzap*pi/180)*cos((azimuth_sun-azimuth_view)*pi/180));

Ss=V/ta;% tas
Lo=l*b/(d*ta);% p_gap_ax_b.c�е�
% Lo=l*pi*r^2/ta;  % QFeng,20190327
hs=r/(Ss*Lo);   % ls.c�еģ���������������������������������������������������������������
lambda_ms=tan(xi)*hs; %ls.c�е�

LAI_HOT_SPOT=1;% Ϊʲô==1  ��Ltһ���������ڸ��Ƕ�ָ��������Ҷ���ָ��
Ls=(a+c*szap*pi/180)*omegae*LAI_HOT_SPOT; %ls.c�е� shoot�߶Ȳ�����gamae

if xi<=pi/2  %������������������������������������
    IN1=0; % ���ӻ���
    IN2=0; % ��ĸ����
    FLAG_IN1=0;
    FLAG_IN2=0;
    MAX=40000; % ���ֵ�������
    INC=0.001; % ���ּ��
    
    for i=0:MAX % �����
        if (lambda_ms+INC*i)==0 % lambda_ms==0
            IN2=1; % ��㸳ֵ
            IN1=0;
            break;
        else % lambda_ms~=0
            IN1=IN1+exp(-Ls*(1+(lambda_ms+INC*i)/Ws))/atan((lambda_ms+INC*i)/hs); % �������Ѿ�����ͬʱԼ��Լ����(Ls/Ws)��һ��
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
        Fs=1-IN1/IN2*xi; % �����ȵ㺯��
    end
    if xi<=0.000001
        Fs=1;
    end
    PT=ptf+(1-Hot-ptf)*Fs; % ���������.cһ����ֻ�Ǽ�һ��һ��λ�ò�ͬ������ɼ�����ҶƬ
%     PT=ptf+(1-ptf)*Fs;
else
    Fs=0;  %     �޸ģ�����������������������������������������
    PT=ptf;
end