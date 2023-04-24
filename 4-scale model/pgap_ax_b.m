% pgapval=pgap_ax_b(0.5,0,0.8,1.41,3.5,3500,10000,zap,0.45,6.5,13,1) % 输入角度
% 计算树冠内的孔隙率(光线/视线穿过一棵树的空隙到达地面的概率)

function pgapval=pgap_ax_b(a,c,omegae,gamae,LAI_canopy,d,b,zap,r,hb,alpha,cors)
% a和c是计算平均叶倾角的参数(其中c等于论文中的b),omegae嫩枝聚集度指数,gamae针对枝的面积比例,n是像元中的样地数量
% l叶面积指数,d像元内树冠总数,b像元面积(投影面积),za太阳天顶角或视线天顶角,r筒半径，hb筒高,alpha锥顶角的一半
% shape==1针叶形，shape==2椭球形

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