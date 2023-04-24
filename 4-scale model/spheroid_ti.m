% spheroid_ti(vzap,30,6.5,0.45,0) % 输入角度
% 视线看到的光照椭球型树冠面积(不考虑冠层中的空隙)
% 不会推导里面的公式，直接抄居老师给的程序？？？？？？？？？？？？？？？

function ti=spheroid_ti(vzap,szap,hb,r,phi_r_sv)
% vza观测天顶角,sza太阳天顶角,hb椭球垂直长,r椭球的水平半径,phi相对方位角

vza=vzap*pi/180; % 角度转弧度
sza=szap*pi/180;
phi=phi_r_sv*pi/180;

vza_prime=atan(hb/(2*r)*tan(vza));
sza_prime=atan(hb/(2*r)*tan(sza));
cs_prime=cos(sza_prime)*cos(vza_prime)+sin(sza_prime)*sin(vza_prime)*cos(phi);

vza=vza*180/pi; % spheroid_ta需要输入角度，因此需要这一步转换
ti=spheroid_ta(vza,hb,r)*0.5*(1+cs_prime);