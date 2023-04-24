clear all;
% load('testdata.mat');

LAI_canopy=3.5; % LAI of a sletect tree;	һ������LAI
b=1000;    % b is the pixel area ;	��Ԫ���m2
d=186;     % d is the tree number;          ��Ԫ�е�����
n=10;       % n is the plot number; 	    ��Ԫ�е���������
m2=3;      % m2	is the average number of plants in a community; һ��Ⱥ���ƽ������
k=350;      % k is the Maximum number of plants in a plot (specified as 350);	��������������е�����(�涨Ϊ350)
norp=2;     % norp is Newman distribution (1) or Poisson distribution (2);	Ŧ���ֲ�(1)���ǲ��ɷֲ�(2)
shape=2;    % shape is Cylindrical(1)Ellipsoid(2);	Ͳ��(1)������(2)
r=1.09;      % r is canopy radius;	���ڰ뾶
ha=2;       % ha is the stump height;	��׮��
hb=2.32;     % hb is the vertical height of spheroid	Ͳ�߻�����ֱ��
alpha=20;  % alpha is top angle of a half cone;	��׶����
omegae=0.8;   % omegae is clumping index of Shoots	��֦�ۼ���ָ��
gamae=1.1; % gamae	is area ratio for branches; ���֦���������
cp=0.8; % cp phase function, 0.9 for NIR, 0.35 for red
Ws=0.025; % Element width(m)
a=0.5;      % a	is Average leaf inclination parameters; ƽ��Ҷ��ǲ�����������G(theta)�õģ���a=0.5,c=0,��G(theta)=0.5;
c=0;        % c	is Mean leaf inclination parameters, which is depicted by b in the paper. ƽ��Ҷ��ǲ����������е�b
% Rleaf ҶƬ�����ʣ����ļ�����
% Tleaf ҶƬ͸���ʣ����ļ�����
% Rground���汳�������ʣ����ļ�����
SdifSg=0.15;         % SdifSg is the proportion of scattered radiation; ɢ��������
sza_saa = load("\\Jinlab309\d\JinShichao\3DSIF\Malixia_14Layer30Points\ToMaliXia\XYZ_NPD0816\sza_saa_out.txt");
szap= sza_saa(:,9); % load('E:\NJU\3DSIF_Zch\���&����\SCOPE_USHa1\SCOPE_v1.61\���ɢ��\ZA.txt');            % szap	       Solar zenith angle�� ̫���춥�ǣ��Ƕ�
azimuth_sun=sza_saa(:,10); % load('E:\NJU\3DSIF_Zch\���&����\SCOPE_USHa1\SCOPE_v1.61\���ɢ��\AA.txt');     % azimuth_sun Solar azimuth; ̫����λ��
azimuth_dem=0;      % azimuth_dem	Slope azimuth; б�·�λ�ǣ�0��ʾ����
theta_g=0;          % theta_g	Slope inclination; б�����

% Rground0=load('E:\G-6-SHB\G-6-SHB2018\2-����������\BackgroundRefl_2018Jul_Za_400900_1nm.txt'); % 1 nm interval
% rG=Rground0(1:5:end,2); % 5 nm interval

% % PROSPECT-D
% N=3.5;
% Cab=40;
% Car=4;
% Canth=2;
% Cw=0.006;
% Cm=0.007;
% RT=prospect_DB(N,Cab,Car,Canth,0,Cw,Cm);
% Rleaf=RT(1:5:501,2)*0.8; % 400-900 5nm
% Tleaf=RT(1:5:501,3)*0.5; % 400-900 5nm
% wvl=(400:5:900)';
% 
% directly input Rleaf Tleaf rG
Rleaf=0.359745336498856; % 737nm
Tleaf=0.411655848061785;
rG=0;
if shape==1
    hc=r/tan(alpha*pi/180); % ����׶�ĸ߶�   
else
    hc=0;   
end
xtheta=180;
RT_MS_total=zeros(601,1);
RZT_MS_total=zeros(601,1);
try1=zeros(601,1);
try2=zeros(601,1);
for time=1:601
za=szap(time)/180*pi;
aa=azimuth_sun(time)/180*pi;
azimuth_view=aa+xtheta;
intensity=atan((hb+hc)/(2*r));
intensity=abs(intensity-za);
intensity=abs(cos(intensity));
try1(time)=intensity;
[RT_MS,RZT_MS,RG_MS,RZG_MS]=multi_scattering_4scale2001_QFeng(Rleaf,Tleaf,rG,za,LAI_canopy,SdifSg,b,d,n,m2,k,norp,shape,alpha,ha,hb,r,omegae,gamae,cp,Ws,a,c,azimuth_dem,azimuth_view,aa,theta_g);
try2(time)=RT_MS;
%RT_MS_total(time)=Rleaf./pi.*intensity+RT_MS;
RT_MS_total(time)=Rleaf+RT_MS;
RZT_MS_total(time)=RZT_MS;
end
% %VZA=(-85:5:85)';
% hold on;plot(VZA,Rcanopy(2,:));
% % figure;plot(wvl,Rcanopy(:,1),wvl,Rleaf,wvl,Tleaf,wvl,rG);
% x=Rcanopy';

close all;




