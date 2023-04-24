clear all;
% load('testdata.mat');

LAI_canopy=3.5; % LAI of a sletect tree;	一株树的LAI
b=1000;    % b is the pixel area ;	像元面积m2
d=186;     % d is the tree number;          像元中的株数
n=10;       % n is the plot number; 	    像元中的样地数量
m2=3;      % m2	is the average number of plants in a community; 一个群落的平均株数
k=350;      % k is the Maximum number of plants in a plot (specified as 350);	样地中最大允许有的株数(规定为350)
norp=2;     % norp is Newman distribution (1) or Poisson distribution (2);	纽曼分布(1)还是泊松分布(2)
shape=2;    % shape is Cylindrical(1)Ellipsoid(2);	筒形(1)椭球形(2)
r=1.09;      % r is canopy radius;	树冠半径
ha=2;       % ha is the stump height;	树桩高
hb=2.32;     % hb is the vertical height of spheroid	筒高或椭球垂直长
alpha=20;  % alpha is top angle of a half cone;	半锥顶角
omegae=0.8;   % omegae is clumping index of Shoots	嫩枝聚集度指数
gamae=1.1; % gamae	is area ratio for branches; 针对枝的面积比例
cp=0.8; % cp phase function, 0.9 for NIR, 0.35 for red
Ws=0.025; % Element width(m)
a=0.5;      % a	is Average leaf inclination parameters; 平均叶倾角参数，用于求G(theta)用的，若a=0.5,c=0,则G(theta)=0.5;
c=0;        % c	is Mean leaf inclination parameters, which is depicted by b in the paper. 平均叶倾角参数，论文中的b
% Rleaf 叶片反射率，由文件导入
% Tleaf 叶片透射率，由文件导入
% Rground地面背景反射率，由文件导入
SdifSg=0.15;         % SdifSg is the proportion of scattered radiation; 散射辐射比例
sza_saa = load("\\Jinlab309\d\JinShichao\3DSIF\Malixia_14Layer30Points\ToMaliXia\XYZ_NPD0816\sza_saa_out.txt");
szap= sza_saa(:,9); % load('E:\NJU\3DSIF_Zch\结果&代码\SCOPE_USHa1\SCOPE_v1.61\多次散射\ZA.txt');            % szap	       Solar zenith angle； 太阳天顶角，角度
azimuth_sun=sza_saa(:,10); % load('E:\NJU\3DSIF_Zch\结果&代码\SCOPE_USHa1\SCOPE_v1.61\多次散射\AA.txt');     % azimuth_sun Solar azimuth; 太阳方位角
azimuth_dem=0;      % azimuth_dem	Slope azimuth; 斜坡方位角，0表示北坡
theta_g=0;          % theta_g	Slope inclination; 斜坡倾角

% Rground0=load('E:\G-6-SHB\G-6-SHB2018\2-背景反射率\BackgroundRefl_2018Jul_Za_400900_1nm.txt'); % 1 nm interval
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
    hc=r/tan(alpha*pi/180); % 计算锥的高度   
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




