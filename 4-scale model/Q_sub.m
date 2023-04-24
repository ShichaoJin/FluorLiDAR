% Ptreev=Q_sub(10000,40,1,1,3500,5,350,0.75,6.5,13,vzap,azimuth_dem,azimuth_view,theta_g)
% Subroutine that calculates the gap fraction (Pvg and Pig)

function [Ptreev]=Q_sub(b,n,norp,cors,d,m2,k,r,hb,alpha,zap,azimuth_dem,azimuth_view_sun,theta_g) % ����Ƕ�
% b����Ԫ���,n����Ԫ�е���������,norp,corsֻ�ܵ���1��2,d����Ԫ�е�����,m2��һ��Ⱥ���ƽ������,
% k����������������е�����(�涨Ϊ350),rͲ�뾶�������ˮƽ�뾶,za̫���춥�ǻ�۲��춥��(�Ƕ�),
% hbͲ�߻�����ֱ��,alpha׶���ǵ�һ��,a��c�Ǽ���ƽ��Ҷ��ǵĲ���(����c���������е�b),
% omegae��֦�ۼ���ָ��,gamae���֦���������,lҶ���ָ��,theta_g�¶�,azimuth_view�۲ⷽλ��

za=acos(cos(theta_g*pi/180)*cos(zap*pi/180)+sin(theta_g*pi/180)*sin(zap*pi/180)*cos(azimuth_view_sun*pi/180-azimuth_dem*pi/180));

if norp==1
    px=neyman(d,n,m2,k);
end
if norp==2
    px=poissonf(d,n,k);
end

if cors==1
    [tab,tac]=cone_taborc(r,zap,hb,alpha);
    ta=tab+tac;
end
if cors==2
    ta=spheroid_ta(zap,hb,r);
end

NN=k;
A=b/n*cos(za)/cos(theta_g*pi/180); % б�����
i_m=d/n;

for j=1:NN
    Pt=0;
    for i=j:NN % ѭ�� Ptj
        Ptj_i(i)=0;
        if i>i_m && Ptj_i(i-1)<0.000000001 && i~=j % ���i-1����������ס���߻����ߵĸ��ʺ�С����������ƽ��ˮƽ������Ptj(i-1)����Ptj(i)���ü���ֱ�Ӹ�ֵ
            Ptj_i(i)=0;
        else
            aa=ta; % ��ʱ���� % Vg or Sg
            if i>i_m % ������׶Ͳ������ͶӰ����Ŀ�������������Ŀ��ƿ���ʹ��������ĸ��ʱ�ıȲ���������ƴ�һ��
                aa=ta*i_m/i; % �����������ƽ��ˮƽ����sg or vg����������
            end
            if ta>A % ������׶Ͳ������ͶӰ����Ŀ�������������Ŀ��Ʋ������ã���û��ǰ��Ŀ���ʱ������úܴ�
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
        Pt=Ptj_i(i)+Pt;
    end
    Ptreev(j)=Pt;
end