% spheroid_ta(zap,6.5,0.45) % ����Ƕ�
% ���߿������������������(�����ǹڲ��еĿ�϶)
% �����Ƶ�����Ĺ�ʽ��ֱ�ӳ���.c�ĳ��򣿣���������������������������

function tab=spheroid_ta(zap,hb,r)
% vza�۲��춥��,hb����ֱ��,r�����ˮƽ�뾶

zap=zap*pi/180; % �Ƕ�ת����
tab=pi*r*(hb/2*sin(zap)+r*cos(zap));

% .c��������Ĺ�ʽδд