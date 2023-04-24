% pgapval=pgap_ax_b(0.5,0,0.8,1.41,3.5,3500,10000,zap,0.45,6.5,13,1) % ����Ƕ�
% ���������ڵĿ�϶��(����/���ߴ���һ�����Ŀ�϶�������ĸ���)

function pgapval=pgap_ax_b(a,c,omegae,gamae,LAI_canopy,d,b,zap,r,hb,alpha,cors)
% a��c�Ǽ���ƽ��Ҷ��ǵĲ���(����c���������е�b),omegae��֦�ۼ���ָ��,gamae���֦���������,n����Ԫ�е���������
% lҶ���ָ��,d��Ԫ����������,b��Ԫ���(ͶӰ���),za̫���춥�ǻ������춥��,rͲ�뾶��hbͲ��,alpha׶���ǵ�һ��
% shape==1��Ҷ�Σ�shape==2������

if cors==1
    [tab,tac]=cone_taborc(r,zap,hb,alpha);
    ta=tab+tac;
end

if cors==2
    ta=spheroid_ta(zap,hb,r);
end

LAI_tree=LAI_canopy.*b./(d.*pi.*r.^2); % QiuFeng, 20191016
% LAI_tree=LAI_canopy; % QiuFeng, 20191016
% pgapval=exp(-(a+c*zap*pi/180)*omegae*l*b/(gamae*d*ta));
pgapval=exp(-(a+c*zap*pi/180)*omegae*LAI_tree*pi*r^2/(gamae*ta));