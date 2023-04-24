% spheroid_ta(zap,6.5,0.45) % 输入角度
% 视线看到的椭球型树冠面积(不考虑冠层中的空隙)
% 不会推导里面的公式，直接抄的.c的程序？？？？？？？？？？？？？？？

function tab=spheroid_ta(zap,hb,r)
% vza观测天顶角,hb椭球垂直长,r椭球的水平半径

zap=zap*pi/180; % 角度转弧度
tab=pi*r*(hb/2*sin(zap)+r*cos(zap));

% .c中求体积的公式未写