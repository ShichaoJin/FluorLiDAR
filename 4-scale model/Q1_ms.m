% compute Q1 eq.(25) in Chen,2001
function Q1_theta=Q1_ms(theta,SZA,phi,LAIH,R,Hb,Hc,cp,shape)
G=0.5;
% cp=0.75;
% cs=G/cos(SZA);% wrong!!!
% cv=G/cos(theta);% wrong!!!

% % crown volume
% crownV=2/3*pi*R*R*Hb;  %shape=2;
% Hindex=B/(crownV*D)*pi*R/2;  % 4scale.c÷–shape=2 ±”–¥ÌŒÛ
% LAIH=LAI*Hindex;



if theta<0
    theta=-theta;
else
end
    
diff=XI(SZA,theta,phi);
if diff<0
    diff=-diff;
else
end

del=1-cp*diff/pi;

if del<0
    del=-del;
elseif del>1
    del=1;
end  

if theta>pi/2
    theta=pi-theta;
elseif theta==0
    theta=0.0000000001;
elseif sqrt((theta-pi/2)^2)<0.0000000001 
    theta=pi/2-0.0000000001;
else
end

if shape == 1 % needle£¨◊∂–Œ+‘≤÷˘
    V = pi*R*R*(Hb+Hc/3); 
else       % shape == 2 ¿´“∂£¨Õ÷«Ú–Œ
    V = 2/3*pi*R*R*Hb; 
end
% mu=LAI*B/(D*V);


cv=G/sin(theta);

if SZA==0
    cs=G/sin(SZA+0.1/180*pi);
else
    cs=G/sin(SZA);
end
    

% vg_theta=2*tan(theta)*R*Hb+pi*R*R;
% PgapV_theta=exp(-cv*LAI*B/(D*vg_theta*cos(theta)));
% PgapV_theta=exp(-G*LAI*B/(D*vg_theta*cos(theta)));
% PgapV_theta=exp(-G*LAI/cos(theta)); 
PgapV_theta=exp(-G*LAIH/sin(theta)); 

Q1_theta=del*(1-exp(-LAIH*(cs+cv)))*(cs*cv/(cs+cv))/(1-PgapV_theta);
% Q1_theta=(1-exp(-LAIH*(cs+cv)))*(cs*cv/(cs+cv))/(1-PgapV_theta);




