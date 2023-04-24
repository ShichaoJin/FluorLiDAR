% pti=pti_adv(10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g)
% Subroutine that calculates the gap fraction (Pvg and Pig)

function [pti]=pti_adv(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,vzap,szap,azimuth_dem,azimuth_view,azimuth_sun,theta_g) % ����Ƕ�
% b����Ԫ���,n����Ԫ�е���������,norp,corsֻ�ܵ���1��2,d����Ԫ�е�����,m2��һ��Ⱥ���ƽ������,
% k����������������е�����(�涨Ϊ350),rͲ�뾶�������ˮƽ�뾶,za̫���춥�ǻ�۲��춥��(�Ƕ�),
% hbͲ�߻�����ֱ��,alpha׶���ǵ�һ��,a��c�Ǽ���ƽ��Ҷ��ǵĲ���(����c���������е�b),
% omegae��֦�ۼ���ָ��,gamae���֦���������,lҶ���ָ��,theta_g�¶�,azimuth_view�۲ⷽλ��

sza=acos(cos(theta_g*pi/180)*cos(szap*pi/180)+sin(theta_g*pi/180)*sin(szap*pi/180)*cos(azimuth_sun*pi/180-azimuth_dem*pi/180));
vza=acos(cos(theta_g*pi/180)*cos(vzap*pi/180)+sin(theta_g*pi/180)*sin(vzap*pi/180)*cos(azimuth_view*pi/180-azimuth_dem*pi/180));

% ѡ������������ĳ�ֲַ���Ŧ���ֲ�(1)���ǲ��ɷֲ�(2)
if norp==1
    px=neyman(d,n,m2,k);
end
if norp==2
    px=poissonf(d,n,k);
end

% ѡ��������ֲ������״��Ͳ��(1)����������(2)������crown projection on the ground
if cors==1
    [tab_v,tac_v]=cone_taborc(r,vzap,hb,alpha);
    ta_v=tab_v+tac_v;
    [tab_s,tac_s]=cone_taborc(r,szap,hb,alpha);
    ta_s=tab_s+tac_s;
end
if cors==2
    ta_v=spheroid_ta(vzap,hb,r);
    ta_s=spheroid_ta(szap,hb,r);
end

pgapv=pgap_ax_b(a,c,omegae,gamae,l,d,b,vzap,r,hb,alpha,cors); % ����һ�����Ŀ�϶��
pgaps=pgap_ax_b(a,c,omegae,gamae,l,d,b,szap,r,hb,alpha,cors); % ����һ�����Ŀ�϶��
A=b/n/cos(theta_g*pi/180); % б�����
i_m=d/n; % ����ƽ������

for choice=1:2
    if choice==1
        aaa=ta_s; % ��������tas
        za=sza;
        P_s_up=0;
        P_s_down=0;
    end
    if choice==2
        aaa=ta_v; % ��������tav
        za=vza;
        P_v_up=0;
        P_v_down=0;
    end
    for j=1:k % ѭ�� Pvg,j==0ʱֻ����pvg_c
        Ptj=0;
        for i=j:k % ѭ�� Ptj
            Ptj_i(i)=0;
            if i>i_m && Ptj_i(i-1)<0.000000001 && i~=j % ���i-1����������ס���߻����ߵĸ��ʺ�С����������ƽ��ˮƽ������Ptj(i-1)����Ptj(i)���ü���ֱ�Ӹ�ֵ
                Ptj_i(i)=0;
            else
                aa=aaa; % ��ʱ���� % Vg or Sg
                if i>i_m % ������׶Ͳ������ͶӰ����Ŀ�������������Ŀ��ƿ���ʹ��������ĸ��ʱ�ıȲ���������ƴ�һ��
                    aa=aaa*i_m/i; % �����������ƽ��ˮƽ����sg or vg����������
                end
                if aaa>A % ������׶Ͳ������ͶӰ����Ŀ�������������Ŀ��Ʋ������ã���û��ǰ��Ŀ���ʱ������úܴ�
                    aa=A;  % �������Ƕ��µ���Ӱ������������ص�����������������������Ӱ���
                end
                PNN=px(i+1); % ��ʱ��������������i�����ĸ���
                if i~=j
                    for v=1:(i-j) % ѭ���������ֲ�
                        PNN=PNN*(j+v)/v*(1-aa/(A*cos(za)));
                    end
                end
                Ptj_i(i)=PNN*(aa/(A*cos(za)))^j;
            end
            Ptj=Ptj_i(i)+Ptj;
        end
        %���¿�ʼ����Ptj֮���
        if choice==1
            P_s_up=P_s_up+Ptj*(1-pgaps^j);
            P_s_down=P_s_down+Ptj*(1-pgaps)*j;
        end
        if choice==2
            P_v_up=P_v_up+Ptj*(1-pgapv^j);
            P_v_down=P_v_down+Ptj*(1-pgapv)*j;
        end
    end
end

Qs=P_s_up/P_s_down; % ���������б��������ƽ��һ��
Qv=P_v_up/P_v_down; % ���������б��������ƽ��һ��
Qsv=min(Qs,Qv); % �����б��������ƽ��һ��

pigneyman=overlap_pvgorpig(b,n,1,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g); % neymanʱ��pig
pigpoisson=overlap_pvgorpig(b,n,2,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,szap,azimuth_dem,azimuth_sun,theta_g);% poissonʱ��pig
OmegaT = log(pigneyman)/log(pigpoisson); % ��eΪ��
% Ϊʲô��pig������pvg���ۼ�������̫���Ƕȸı��𣿲���۲�Ƕȸı��𣿣�������������������������������������������������������
Lt=OmegaT*pi*r*r*d/b; % �ڴ�ֱ�����ϵ�ͶӰ
Wt=(pi*r*r)^0.5; % �ڴ�ֱ�����ϵ�ͶӰ
E_r=Wt/Lt;
thetae=atan(2*r/E_r); % ������������в�ͬ!!��б���ϵĽǶȣ���������������������������
if thetae<0 % ��֤theta��С��0
    thetae=-thetae;
end

phi_r_sv=abs(azimuth_sun-azimuth_view);
if phi_r_sv>180
    phi_r_sv=360-phi_r_sv;
end
F=1-phi_r_sv*pi/180/thetae;%????????????????????????????????????????
if F<0
    F=0;
end
if F>1
    F=1;
end

Qcb=Qv*Qs+F*(Qsv-Qv*Qs);

if cors==1
    [tib,tic]=cone_ticorb(szap,vzap,alpha,phi_r_sv,r,hb);
    ti=tic+tib;
    pti=Qcb*ti/(Qv*ta_v);
end
if cors==2
    ti=spheroid_ti(vzap,szap,hb,r,phi_r_sv);
    pti=Qcb*ti/(Qv*ta_v);
end