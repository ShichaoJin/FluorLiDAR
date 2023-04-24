% spheroid_ti(vzap,30,6.5,0.45,0) % ����Ƕ�
% ���߿����Ĺ����������������(�����ǹڲ��еĿ�϶)
% �����Ƶ�����Ĺ�ʽ��ֱ�ӳ�����ʦ���ĳ��򣿣���������������������������

function ti=spheroid_ti(vzap,szap,hb,r,phi_r_sv)
% vza�۲��춥��,sza̫���춥��,hb����ֱ��,r�����ˮƽ�뾶,phi��Է�λ��

vza=vzap*pi/180; % �Ƕ�ת����
sza=szap*pi/180;
phi=phi_r_sv*pi/180;

vza_prime=atan(hb/(2*r)*tan(vza));
sza_prime=atan(hb/(2*r)*tan(sza));
cs_prime=cos(sza_prime)*cos(vza_prime)+sin(sza_prime)*sin(vza_prime)*cos(phi);

vza=vza*180/pi; % spheroid_ta��Ҫ����Ƕȣ������Ҫ��һ��ת��
ti=spheroid_ta(vza,hb,r)*0.5*(1+cs_prime);