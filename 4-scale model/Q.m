% Ϊʲô���ֵ������ߣ�ΪʲôQ1��Q2����С����������������������������������������������������������������
% [QQ1 QQ2]=Q(10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g)
% Subroutine that calculates the gap fraction (Pvg and Pig)

function [QQ1 QQ2]=Q(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,cp,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g) % ����Ƕ�
% b����Ԫ���,n����Ԫ�е���������,norp,corsֻ�ܵ���1��2,d����Ԫ�е�����,m2��һ��Ⱥ���ƽ������,
% k����������������е�����(�涨Ϊ350),rͲ�뾶�������ˮƽ�뾶,za̫���춥�ǻ�۲��춥��(�Ƕ�),
% hbͲ�߻�����ֱ��,alpha׶���ǵ�һ��,a��c�Ǽ���ƽ��Ҷ��ǵĲ���(����c���������е�b),
% omegae��֦�ۼ���ָ��,gamae���֦���������,lҶ���ָ��,theta_g�¶�,azimuth_view�۲ⷽλ��

pgapv=pgap_ax_b(a,c,omegae,gamae,l,d,b,vzap,r,hb,alpha,cors); % ����һ�����Ŀ�϶��

if cors==1
    hc=r/tan(alpha*pi/180); % ����׶�ĸ߶�
    V=pi*r^2*hb+1/3*pi*r^2*hc; % �������ڵ����
    Lo_90=l*b/(V*d)*pi*r/2; % һ����ˮƽ����Ҷ���ָ�����ܶ�*·������
%     Lo_90=l.*pi.*r./(hb*2+hc); % QiuFeng,20190327
    [tab_v,tac_v]=cone_taborc(r,vzap,hb,alpha);
    ta_v=tab_v+tac_v;
    [tab_s,tac_s]=cone_taborc(r,szap,hb,alpha);
    ta_s=tab_s+tac_s;
    Sv=V/ta_v; % һ������б·������,����׼ȷ��ֵ���൱�ھ�ֵ
    Ss=V/ta_s;
    mu=l*b/(d*V);
%     mu=l/V; % QiuFeng,20190327
    Cs=(a+c*szap*pi/180)*omegae/gamae*Ss*mu/Lo_90; % �൱�� Cs=G(sza)*omegae/gamae*s(sza)/s(90)
    Cv=(a+c*vzap*pi/180)*omegae/gamae*Sv*mu/Lo_90;
end
if cors==2
    V=4/3*pi*r^2*(hb*0.5); % ������������
    Lo_90=l*b/(V*d)*pi*r/2;%
%     Lo_90=l.*2.*r./(hb); % QiuFeng,20190327
    ta_v=spheroid_ta(vzap,hb,r); %%% ???
    ta_s=spheroid_ta(szap,hb,r); %%% ???
    Sv=V./ta_v;
    Ss=V./ta_s;
    mu=l/V;
    Cs=(a+c*szap*pi/180)*omegae/gamae*Ss*mu/Lo_90; 
    Cv=(a+c*vzap*pi/180)*omegae/gamae*Sv*mu/Lo_90;
%     Cs=(a+c*szap*pi/180)*omegae/gamae/sin(szap*pi/180); 
%     Cv=(a+c*vzap*pi/180)*omegae/gamae/sin(vzap*pi/180);

end

QQ1=0;
QQ2=0;
Q1=(1-exp(-(Cs+Cv)*Lo_90))*Cs*Cv/(Cs+Cv);
if abs(Cs-Cv)<0.0000001
    Cs=Cs-0.0000001;
end
Q2=(exp(-Cs*Lo_90)-exp(-Cv*Lo_90))*Cs*Cv/(Cv-Cs);
sza=acos(cos(theta_g*pi/180)*cos(szap*pi/180)+sin(theta_g*pi/180)*sin(szap*pi/180)*cos(azimuth_sun*pi/180-azimuth_dem*pi/180));
vza=acos(cos(theta_g*pi/180)*cos(vzap*pi/180)+sin(theta_g*pi/180)*sin(vzap*pi/180)*cos(azimuth_view*pi/180-azimuth_dem*pi/180));

[sgap spvg_c]=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g);
[vgap vpvg_c]=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,vzap,azimuth_dem,azimuth_sun,theta_g);

aa=(a+c*vzap*pi/180)*omegae/gamae*l*b/(d*ta_v)*(1-sgap)*(1-vgap);%spvg_c++vpvg_c*cos(vza)/cos(sza)???????*(1-vgap)???????????????????????��manual�ϲ�ͬ???????????????????????????????vza and sza wrong
% aa2=(a+c*vzap*pi/180)*omegae/gamae*l*b/(d*ta_v);%*(1-sgap+spvg_c)*(1-vgap+vpvg_c)���Ƕ�Խ�󴩹������²�ĸ���ԽС���ӹ��������Ӱ�濴������ҶƬ�ĸ���ԽС

Ptreev=Q_sub(b,n,norp,cors,d,m2,k,r,hb,alpha,vzap,azimuth_dem,azimuth_view,theta_g);
NN=k;
for i=1:NN % ��i������ס�ĸ���
    Pav_i=0;
    if exp(-(i-1)*a) <0.00000000000000000000000001
        i=NN+1;
    end
    for j=i:NN % ������i������ס�ĸ�����ͣ��õ�Pat(i)
        Pav_i=Pav_i+Ptreev(j);
        if Ptreev(j)<0.00000000000000000000000001 && j>NN/2
            j=NN+1;
        end
    end
    QQ1=QQ1+Q1*Pav_i*exp(-(i-1)*aa);%������������������������������������������������������*pgapv^(i-1)
    QQ2=QQ2+Q2*Pav_i*exp(-(i-1)*aa);%!!!!!!!!!!!!!!!!!!!!!!*pgapv^(i-1)
end

xi=acos(cos(szap*pi/180)*cos(vzap*pi/180)+sin(szap*pi/180)*sin(vzap*pi/180)*cos((azimuth_sun-azimuth_view)*pi/180));
QQ1=QQ1*(1-cp.*xi/pi);% ��������������������������������������������������������������������������������
QQ2=QQ2*(1-cp.*xi/pi);% 1==cp ��������������������������������������������������������������������������������

if QQ1>1
    QQ1=1;
end
if QQ1<0
    QQ1=0;
end
if QQ2>1
    QQ2=1;
end
if QQ2<0
    QQ2=0;
end
