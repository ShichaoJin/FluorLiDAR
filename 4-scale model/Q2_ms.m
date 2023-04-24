% compute Q2 eq.(25) in Chen,2001
function Q2_theta=Q2_ms(theta,SZA,phi,LAIH,R,Hb,cp)
G=0.5;
% cp=0.75;

% diff=XI(SZA,theta,phi);   % phi=0还是pi?
% del=1-cp*diff/pi;
% if del<0
%     del=-del;
% elseif del>1
%     del=1;
% end
% 
% cv=G/sin(theta);
% cs=G/sin(SZA);
% 
% % cv_mean=G/sin(theta);
% cv_mean=G*Sv;
% TG=tan(theta+pi/2);
% if TG<0 
%     TG=-TG;
% else
% end



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

cv=G/sin(theta);
if SZA==0
    cs=G/sin(SZA+0.1/180*pi);
else
    cs=G/sin(SZA);
end

PgapV_theta=exp(-G*LAIH/sin(theta)); 
% PgapV_theta=exp(-G*LAI/cos(theta)); 
% vg_mean=2*tan(theta)*R*Hb+pi*R*R;
% PgapV_mean=exp(-cv*LAI*B/(D*vg_mean*cos(theta))); %原程序中用的是cv_mean

Q2_theta=del*(exp(-LAIH*cs)-exp(-LAIH*cv))*(cs*cv)/(cv-cs)/(1-PgapV_theta);

