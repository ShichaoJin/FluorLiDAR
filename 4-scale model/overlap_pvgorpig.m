% [pig]=overlap_pvgorpig(10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,szap,azimuth_dem,azimuth_view_or_sun,theta_g)
% [pvg]=overlap_pvgorpig(10000,40,1,1,3500,5,350,0.75,6.5,13,0.5,0,0.8,1.41,3.5,vzap,azimuth_dem,azimuth_view_or_sun,theta_g)
% Subroutine that calculates the gap fraction (Pvg and Pig)

function [pvgorpig pvg_c]=overlap_pvgorpig(b,n,norp,cors,d,m2,k,r,hb,alpha,a,c,omegae,gamae,l,zap,azimuth_dem,azimuth_view_or_sun,theta_g) % ����Ƕ�
% b����Ԫ���,n����Ԫ�е���������,norp,corsֻ�ܵ���1��2,����3��������ֲ�,d����Ԫ�е�����,m2��һ��Ⱥ���ƽ������,
% k����������������е�����(�涨Ϊ350),rͲ�뾶�������ˮƽ�뾶,za̫���춥�ǻ�۲��춥��(�Ƕ�),
% hbͲ�߻�����ֱ��,alpha׶���ǵ�һ��,a��c�Ǽ���ƽ��Ҷ��ǵĲ���(����c���������е�b),
% omegae��֦�ۼ���ָ��,gamae���֦���������,lҶ���ָ��,theta_g�¶�,azimuth_view�۲ⷽλ��

za=acos(cos(theta_g*pi/180)*cos(zap*pi/180)+sin(theta_g*pi/180)*sin(zap*pi/180)*cos(azimuth_view_or_sun*pi/180-azimuth_dem*pi/180));

% ѡ������������ĳ�ֲַ���Ŧ���ֲ�(1)���ǲ��ɷֲ�(2)
if norp==1
    px=neyman(d,n,m2,k);
end
if norp==2
    px=poissonf(d,n,k);
end

% ѡ��������ֲ������״��Ͳ��(1)����������(2)������crown projection on the ground,vg��sg
if cors==1
    [tab,tac]=cone_taborc(r,zap,hb,alpha);
    ta=tab+tac;
end
if cors==2
    ta=spheroid_ta(zap,hb,r);
end

pgap=pgap_ax_b(a,c,omegae,gamae,l,d,b,zap,r,hb,alpha,cors); % ����һ�����Ŀ�϶��
A=b/n*cos(za)/cos(theta_g*pi/180); % б�����
i_m=d/n; % ����ƽ������

% ------------------------------------- 4-scale gap fraction -------------------------------------
% ����pvg_c,��������0������ס���ߵ����
pvg_c=0;
tat=ta; % ��ʱ����
for i=1:k % ��������0������ס����ʱ��1-350���������
    if i>i_m % ������ʱ��Сta % ������׶Ͳ������ͶӰ����Ŀ�������������Ŀ��ƿ���ʹ��������ĸ��ʱ�ıȲ���������ƴ�һ��
        tat=ta*i_m/i;
    end
    if tat>A
        tat=A;
    end
    pvg_c=pvg_c+px(i+1)*(1-tat/A)^i;
end

aaa=ta; % ��������Vg or Sg
P=0;
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
                    PNN=PNN*(j+v)/v*(1-aa/A);
                end
            end
            Ptj_i(i)=PNN*(aa/A)^j;
        end
        Ptj=Ptj_i(i)+Ptj;
    end
    %���¿�ʼ����Ptj֮���
    P=P+Ptj*(pgap^j);
end

pvgorpig=pvg_c+P;
if pvgorpig<0
    pvgorpig=0;
end
if pvgorpig>1
    pvgorpig=1;
end
