% Chen et al. 2001 IETGARS, revised by QiuFeng 2015 & 2018
% function [RT_MS,RZT_MS,RG_MS]=multi_scattering_4scale2001(ref,rG,LAI,SZA,VZA,SdifSg,Ha,Hb,d)
function [RT_MS,RZT_MS,RG_MS,RZG_MS]=multi_scattering_4scale2001_QFeng(ref,trans,rG,szap,LAI,SdifSg,b,d,n,m2,k,norp,cors,alpha,Ha,Hb,r,omegae,gamae,cp,Ws,a,c,azimuth_dem,azimuth_view,azimuth_sun,theta_g)
                                                            
% ref/trans -- leaf reflectance/transmittance
% rG -- background reflectance
% SZA,VZA -- radian
% SdifSg -- 天空中散射光的比例，可直接给出，或在下面计算瑞利散射
% Ha -- stick height
% Hb -- cylinder or spheroid height
% R 树冠半径
% B
% d -- mean distance between trees


G=0.5;
SZA=szap/180*pi;  %%%%
deltaLAI=0.1;   %%%
fd=SdifSg;


if cors==1
    Hc=r/tan(alpha*pi/180); % 计算锥的高度   
else
    Hc=0;   
end

% % /* rayleigh Scattering */
% tauR = 0.008569*pow(lambda[i],-4)*(1-0.0113*pow(lambda[i],-2.)+ 0.00013*pow(lambda[i],-4));
% Rr = (3*tauR+(2-3*cos(in_p.SZA))*(1-exp(-tauR/cos(in_p.SZA))))/(4+3*tauR);
% fd[i] = ((1-Rr) -  exp(-tauR/cos(in_p.SZA)))/(1-Rr);


% % calculating d：mean distance between trees
% out_p->OmegaT = log(out_p->Pig)/log(out_p->Pig_poisson);
% Lt = out_p->OmegaT*PI*in_p.R*in_p.R*in_p.D/in_p.B;
% Wt = sqrt(PI*in_p.R*in_p.R);
% d=Wt/Lt;
OmegaT=omegae;
Lt=OmegaT.*pi.*r^2.*d./b;
Wt=sqrt(pi.*r^2);
distance_out= Wt/Lt;


% calculation of view factors

% sky view factor from ground 
cos_vza_mean=0.537+0.025*LAI;
vza_mean=acos(cos_vza_mean);
pvg_mean=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,Hb,alpha,a,c,omegae,gamae,LAI,vza_mean,azimuth_dem,azimuth_view,theta_g); 
pvgr_mean=exp(-G.*LAI/cos_vza_mean); % equation (21） in Chen and Leblanc, 2001
omega_Total=log(pvg_mean)./log(pvgr_mean);  % equation (20） in Chen and Leblanc, 2001
F_S_g=exp(-G.*omega_Total.*LAI/cos_vza_mean); % equation (19） in Chen and Leblanc, 2001

CV=G.*omega_Total./cos_vza_mean;

Pig=exp(-G*LAI/cos(SZA));

% calculating ft in equation (34)
[x,ft]=PTandFs(Ws,b,n,norp,cors,d,m2,k,r,Hb,alpha,a,c,omegae,gamae,cp,LAI,vza_mean,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g); % 直接用4-scale里面的计算热点函数的方法
% % 用equation(34)后提到的1-()的方法?
% ee=XI(SZA,vza_mean,0);
% Px=1-ee*cp/pi;% 1st-order scattering phase function
% ft2=1-Px;



% sky view factor from foliage
F_s_t=0;
F_s_tRest=0;
P_tot=0;
F_g_T=0;
F_T_ZT=0;
F_TT_ZT=0;
for Li=deltaLAI:deltaLAI:(LAI-deltaLAI)
    arg1=-G*Li.*omega_Total/cos(asin(0.537+0.025*Li));  % Fst(Li) eq.(29)[Chen,2001] %%%
    F_s_t=F_s_t+exp(-CV*Li)*0.5*exp(arg1)*deltaLAI; %sky view factor to all foliage
    F_s_tRest=F_s_tRest+0.5*exp(arg1)*deltaLAI;    %sky view factor to foliage not weighted with height ???
    
    P_tot=P_tot+exp(-CV*Li)*deltaLAI;
    
    arg2=-G*(LAI-Li).*omega_Total./cos(asin(0.537+0.025*(LAI-Li)));
    F_g_T=F_g_T+0.5*exp(arg2)*deltaLAI;  % F_g_t=F_g_T=F_g_ZT ???
    
    inc_theta=deltaLAI;
    thetah_max=Li/LAI*(Hb+Hc/3)/distance_out;
    thetah_min=(1-Li/LAI)*(Hb+Hc/3)/distance_out;
    theta_max=0.5*(pi/2+atan(thetah_max));
    theta_min=0.5*(pi/2+atan(thetah_min));   
    if theta_min>theta_max  % swap min and max
        t=theta_min;
        theta_min=theta_max;
        theta_max=t;
    end
    
    f_tt1=0;
    f_tt2=0;
    f_tt3=0;
%     LHi=Li*Hindex;
%     theta0=0.5*(inc_theta+pi/2);
    theta0=0.5.*(0+pi/2); % QFeng, 20191026
    for theta=theta0:inc_theta:theta_max
        theta_h=2*theta-pi/2;
        theta_h=pi/2-theta_h; % 计算Q1_avg的公式里面的theta是天顶角，而这里的theta不是，互余(Chen,2001,eq(30))
%         y=Q1_ms(theta_h,SZA,0,LAI,LAIH,B,D,R,Hb)*sin(theta)*cos(theta);
%         %???  /*此处计算Q1时用的LAI是总的LAI还是Li?*/
%         %%%sin()cos()里面用的是theta不是theta_h
        diffQ1=XI(SZA,theta_h,0);
        if diffQ1<0
            diffQ1=-diffQ1;
        else
        end
        delQ1=1-cp*diffQ1/pi;
        if delQ1<0
            delQ1=-delQ1;
        elseif delQ1>1
            delQ1=1;
        end  
        y=Q1_ms(theta_h,SZA,0,Li,r,Hb,Hc,cp,cors)*sin(theta)*cos(theta); %Q1里面已经计算了del，后面对p的积分计算的也是del,是否重复了？??故在此去掉del项(2015/1/24) %2015/2/1这里sin()cos()是theta非theta_h
%         f_tt1=f_tt1+y.*inc_theta.*10/180.*pi;  % 为何要*10/180*pi ???
        if theta<=theta_min % 从theta_min到0，再从0到theta_min积分两次 % revised 2015/1/29
            f_tt1=f_tt1+2*y.*inc_theta.*10/180*pi;  %%%    % 为何要*10/180*pi ???
        else
            f_tt1=f_tt1+y.*inc_theta.*10/180*pi;  %
        end
        
        del1=0;
        del2=0;   
        for p=0:10:180  % integration of solar azimuth angle
            phi=p/180*pi;
            inc_p=10/180*pi;
            diff_h=XI(SZA,theta_h,phi);
            del1=del1+(1-cp*diff_h/pi)*inc_p;
            del2=del2+cp*diff_h/pi*inc_p;
        end
        
        f_tt2=f_tt1;
        f_tt2=f_tt2*del1;  %%%
        f_tt3=f_tt1*del2;  %%%

    end
    
    F_T_ZT=F_T_ZT+f_tt2*deltaLAI;
    F_TT_ZT=F_TT_ZT+f_tt3*deltaLAI;
    
end

F_s_t=F_s_t/P_tot;
F_s_tRest=F_s_tRest/LAI;

F_g_T=F_g_T/LAI;
F_G_T=F_g_T*Pig;
F_ZG_T=F_g_T*(1-Pig);
F_G_ZT=F_G_T;

F_T_ZT=F_T_ZT/LAI;
F_TT_ZT=F_TT_ZT/LAI;
F_ZT_T=F_T_ZT;



% Q1 & Q2 (Chen,2001 eq(25)(26))%
theta_halfmean=atan((Ha+0.5*Hb)/(0.5*distance_out));
% pgap_theta=exp(-G*LAI/cos(pi/2-theta_halfmean));

theta1=theta_halfmean+pi/2;
diff1=XI(SZA,theta1,0);
del11=1-diff1*cp/pi;
if del11<0
    del11=-del11;
elseif del11>1
    del11=1;
else
end
Q1_mean=Q1_ms(theta1,SZA,0,LAI,r,Hb,Hc,cp,cors);%为何要+pi/2 ？若用theta_halfmean大了十倍
Q1B_mean=Q1_mean/del11*(1-del11);

if Q1B_mean>1
    Q1B_mean=1;
elseif Q1B_mean<0
    Q1B_mean=0;
else 
end

% theta2=3/2*pi-theta_halfmean;
% theta2=theta_halfmean+pi/2;
% del2=1-cp*XI(SZA,theta_halfmean,pi)/pi;  %是否要判断del<0???a或Q1Q1<0 ???  %%% %%%%

diff22=XI(SZA,theta1,pi);
del22=1-cp*diff22/pi;
if del22<0
    del22=-del22;
elseif del22>1
    del22=1;
else
end

Q2_mean=Q2_ms(theta1,SZA,pi,LAI,r,Hb,cp); %这里应该是用theta_halfmean,若用theta_halfmean+pi/2大了十倍 ？  %%% %phi=0
Q2B_mean=Q2_mean/del22*(1-del22);

if Q2B_mean>1
    Q2B_mean=1;
elseif Q2B_mean<0
    Q2B_mean=0;
else 
end

% sunlit trees view factor from ground
F_T_g=0.5*(Q1_mean+Q2_mean)*(1-F_S_g);
F_TT_g=0.5*(Q1B_mean+Q2B_mean)*(1-F_S_g);

F_T_T=1-F_s_t-F_g_T-F_ZT_T;

F_g_t=F_g_T;
F_T_t=F_T_T+F_T_ZT;



% calculation of the shaded reflectivities
% trans=ref; %！！！！！

% RT_1st=ref.*F_s_t;
% RG_1st=rG.*F_S_g;

% trans=trans/pi;
% ref=ref/pi;


RT_1st=ref;
RG_1st=rG;

% ref=RT_1st;
% ref=ref./pi;
% trans=trans./pi;

RT_2nd=ref.*((ref+trans*ft)*F_T_T+rG*F_G_T+fd*F_s_t);
RZT_2nd=ref.*(ref*F_T_ZT+trans*F_TT_ZT+rG*F_G_ZT+fd*F_s_t);
RG_2nd=rG.*(ref*F_T_g+trans*F_TT_g+fd*F_S_g);
RTG_2nd=(RG_2nd*F_g_t+RZT_2nd*(F_ZT_T)+RT_2nd*F_T_t)/(1-F_s_t);
RTG_3rd=0.5*(ref+trans).*(1-F_s_t).*RTG_2nd;  % ???

% RT_MS_total=RT_1st.*intensity+RT_2nd+RTG_3rd./(1-0.5*(ref+trans)*(1-F_s_t));  % ???
RT_MS=RT_2nd+RTG_3rd./(1-0.5*(ref+trans)*(1-F_s_t));  % ???
RZT_MS=RZT_2nd+RTG_3rd./(1-0.5*(ref+trans)*(1-F_s_t));  % ???
% RG_MS_total=RG_1st+RG_2nd+rG.*RTG_2nd*(1-F_S_g)./(1-0.5*(ref+trans)*(1-F_s_t)); % ???
RG_MS=RG_2nd+rG.*RTG_2nd*(1-F_S_g)./(1-0.5*(ref+trans)*(1-F_s_t)); % ???
RZG_MS=RG_2nd+rG.*RTG_2nd*(1-F_S_g)./(1-0.5*(ref+trans)*(1-F_s_t)); % ???

% % 去掉1/2
% trans=ref;
% RT_2nd=ref.*((ref+trans*ft)*F_T_T+rG*F_G_t+fd*F_s_t);
% RZT_2nd=ref.*(ref*F_T_ZT+trans*F_TT_ZT+rG*F_G_t+fd*F_s_t);
% RG_2nd=rG.*(ref*F_T_g+trans*F_TT_g+fd*F_S_g);
% RTG_2nd=(RG_2nd*F_g_t+RZT_2nd*(1-F_g_t-F_T_T-F_s_t)+RT_2nd*(F_T_T+F_T_ZT))/(1-F_s_t);
% RTG_3rd=(ref+trans).*(1-F_s_t).*RTG_2nd;  % ???
% RT_MS_total=ref+RT_2nd+RTG_3rd./(1-(ref+trans)*(1-F_s_t));  % ???
% RZT_MS=RZT_2nd+RTG_3rd./(1-(ref+trans)*(1-F_s_t));  % ???
% RG_MS_total=rG+RG_2nd+rG.*RTG_2nd*(1-F_S_g)./(1-(ref+trans)*(1-F_s_t)); % ???
% RZG_MS=RG_2nd+rG.*RTG_2nd*(1-F_S_g)./(1-(ref+trans)*(1-F_s_t)); % ???


% 直接求解方法太慢，改用fminsearch
% syms ref RT_2nd RZT_2nd RG_2nd RTG_2nd RTG_3rd eq
% ref_all=zeros(2151,3);
% for i=1:2151
%     rg=rG(i,1);
%     rpt=Rpt(i,1);
%     % if leaf reflectance=transmittance=ref
%     RT_2nd=ref*((ref+ref*ft)*F_T_T+rg*F_G_t+fd*F_s_t);
%     RZT_2nd=ref*(ref*F_T_ZT+ref*F_TT_ZT+rg*F_G_t+fd*F_s_t);
%     RG_2nd=rg*(ref*F_T_g+ref*F_TT_g+fd*F_S_g);
%     RTG_2nd=(RG_2nd*F_g_t+RZT_2nd*(1-F_g_t-F_T_T-F_s_tRest)+RT_2nd*(F_T_T+F_T_ZT))/(1-F_s_tRest);
%     RTG_3rd=0.5*(ref+ref)*(1-F_s_tRest)*RTG_2nd;  % ???
%     % RT_MS=RT_2nd+RTG_3rd./(1-0.5*(ref+ref)*(1-F_s_tRest));  % ???
% 
%     eq=RT_2nd+RTG_3rd./(1-0.5*(ref+ref)*(1-F_s_tRest))-rpt;
%     output=solve(eq,ref);
%     output=double(output);
%     output=output';
% %     output(find(imag(output)~=0))=[]; % remove complex results
%     ref_all(i,:)=output;
% end





