function [find_index,canopycover] = bidirection_cylinder(filename,SIFdata, obs_theta, x0, y0, treetop, r, solar_azimuth)
format short g;
shape = size(SIFdata);
SIFdata_raw = ones(shape(1), shape(2)+1) ;% Add an ID Number column; 添加ID编号列
SIFdata_raw(:,1:3) = SIFdata;
SIFdata_raw(:,4) = 1:shape(1);

%workdir = "\\Jinlab309\d\JinShichao\3DSIF\Malixia_14Layer30Points\leafindex";
dirname = "D:\\SIF_bidirection_temp_cylinder"

%v = "raw_DART"
%% Load data 
%Loda SIF data 荧光数据
%SIFdata_raw = load(sprintf("%s\\%s\\xyz_PAR_001_targetID.txt",workdir,v));
% sensor location, 传感器位置
% x0=1.33;y0=3.98;z0=4.39; %z值在小一点，取光纤最底部
% treetop = 3.15;
% r = 5;

zr = treetop + r;
sensor_raw = [x0, y0, zr];
% Produce a set of observation points; 产生一组观测点
dim = 2; % the side length of observing area; the value should not be too large 值小一点，保证不出界
xres = 0.1; yres=0.1;     % sensor resolution in x and y direction; 传感器分辨率
xs = x0-dim:xres:x0+dim;  
ys = y0-dim:yres:y0+dim;
L = size(xs);
sensor_points = zeros(L(2)*L(2), 3);
k =1
for x = xs
    for y = ys
        sensor_points(k,:) = [x, y, zr];
        k =k+1;
    end
end

%%
fi0 = 0; %The azimuth is 90 after rotation;  旋转后方位角为 90 
% theta0=30;fi0=169.58;
% 太阳位置
%solar_azimuth = 45 %228.06; 
%solar_zenith = 15.91;

%% Data normolization to the origin; 数据标准化到原点 (这一步可以不要。主要是为了移到类似一个场景范围内，可视化后对比看下旋转角度可对）
SIFdata_raw_std = zeros(size(SIFdata_raw));
minx = min(SIFdata_raw(:,1)); miny = min(SIFdata_raw(:,2));
SIFdata_raw_std(:,1) = SIFdata_raw(:,1)-minx;
SIFdata_raw_std(:,2) = SIFdata_raw(:,2)-miny;
SIFdata_raw_std(:,3:4) = SIFdata_raw(:,3:4);
if ~exist(sprintf("%s\\SIFdata_raw_std.txt",dirname))
    save(sprintf("%s\\SIFdata_raw_std.txt",dirname), 'SIFdata_raw_std', '-ascii')
end
%sensor_raw_std = [x0-minx, y0-miny, zr]
sensor_points_std = [sensor_points(:,1)-minx, sensor_points(:,2)-miny,sensor_points(:,3)];
if ~exist(sprintf("%s\\sensor_raw_std.txt",dirname))
    save(sprintf("%s\\sensor_raw_std.txt",dirname), 'sensor_points_std', '-ascii')
end
%% The point cloud and sensor rotate to the main plane; 点云和传感器旋转到主平面
angle_parrel = (90-solar_azimuth) ; %
M = makehgtform('zrotate',angle_parrel*pi/180);  
SIFdata_r = SIFdata_raw(:,1:3)*M(1:3,1:3);
SIFdata_r = [SIFdata_r,SIFdata_raw(:,4)];
%sensor_r = sensor_raw*M(1:3,1:3);
sensor_points_r = sensor_points*M(1:3,1:3); 

%% the rotated data is normalized to the origin; 旋转后数据标准化到原点
SIFdata_r_std = zeros(size(SIFdata_r));
minx=min(SIFdata_r(:,1)) ; miny =min(SIFdata_r(:,2));
SIFdata_r_std(:,1) = SIFdata_r(:,1)-minx;
SIFdata_r_std(:,2) = SIFdata_r(:,2)-miny;
SIFdata_r_std(:,3:4) = SIFdata_r(:,3:4);
% sensor_r_std = zeros(size(sensor_r));
% sensor_r_std(1) = sensor_r(1) - minx;
% sensor_r_std(2) = sensor_r(2) - miny;
% sensor_r_std(3) = sensor_r(3);
sensor_points_r_std = [sensor_points_r(:,1)-minx, sensor_points_r(:,2)-miny, sensor_points_r(:,3)];


%dirname = sprintf("%s\\%s\\%s",workdir,v,num2str(solar_azimuth))

% if ~exist(strcat(dirname,'\',num2str(solar_azimuth)), 'dir')
%        mkdir(dirname)
% end
if ~exist(sprintf("%s\\%s_SIFdata_r_std_%s.txt",dirname,filename,num2str(solar_azimuth)))
    save(sprintf("%s\\%s_SIFdata_r_std_%s.txt",dirname,filename,num2str(solar_azimuth)), 'SIFdata_r_std', '-ascii')
end
if ~exist(sprintf("%s\\%s_sensor_r_std_%s.txt",dirname,filename,num2str(solar_azimuth)))
    save(sprintf("%s\\%s_sensor_r_std_%s.txt",dirname,filename,num2str(solar_azimuth)), 'sensor_points_r_std', '-ascii')
end

%% Generate a new observation location for a given sensor; 生成给定传感器，新的观测位置（绕观测顶部的中心点旋转观测的传感器点阵）
%for obs_theta = -30:30:60
    % center point
%     x_obs = r*sind(obs_theta)+ sensor_points_r_std(:,1);
%     y_obs = sensor_points_r_std(:,2); % y keep the same as the rotated value
%     z_obs = r*cosd(obs_theta)+treetop;
%     sensor_obs = [x_obs, y_obs, z_obs];
    

    %Pan to a given point as the center of the sphere; 平移到给定点为球心
      
    a = size(sensor_points_r_std);
    b = (a(1)+1)/2;
    center = sensor_points_r_std(b,:);
%     minx = min(sensor_points_r_std(:,1));
%     miny = min(sensor_points_r_std(:,2));
    minx = center(:,1);
    miny = center(:,2);
    minz = center(:,3)-r; %The previous data rotation is only rotated in the z-axis, and the z-coordinate does not change; 前面的数据旋转只是在z轴旋转，z坐标没变化

    sensor_points_r_std_temp = [sensor_points_r_std(:,1)-minx,sensor_points_r_std(:,2)-miny,sensor_points_r_std(:,3)-minz]; % Set the crown top to a rotating sphere; 设置树冠顶为旋转球心
    %Rotate the given angle with the center of the sphere; 以球心旋转给定角度
    M = makehgtform('yrotate',-1*obs_theta*pi/180); % Take a given point as the center of the sphere, rotating around the y-axis； 以给定点为球心，围绕y轴旋转（第二个参数是旋转角度，当角度为正数往逆时针转，而我们需要正角度顺时针，所以前面乘以-1）
    sensor_points_obs = sensor_points_r_std_temp(:,1:3)*M(1:3,1:3);
    %Restore the horizontal position before panning； 恢复平移前的水平位置 
    sensor_points_obs = [sensor_points_obs(:,1)+minx, sensor_points_obs(:,2)+miny,sensor_points_obs(:,3)+minz ];
    save(sprintf("%s\\%s_sensor_obs_%s_%s.txt",dirname,filename,num2str(solar_azimuth),num2str(obs_theta)), 'sensor_points_obs', '-ascii')


    %% Filter the canopy point cloud used for lookup; 筛选用于查找的冠层点云
    find_index=zeros(a(1),1);
    row=0;
    tic
    for sensor_obs = sensor_points_obs'
        % Remove the point within 0.4m of the sensor hemisphere 去除传感器半球0.4m范围的点
        %arrdif = SIFdata_(:,1:3)-[x0, y0, z0];
        sensor_obs = sensor_obs';
        arrdif = SIFdata_r_std(:,1:3)-sensor_obs;
        dis = sqrt(arrdif(:,1).^2 +arrdif(:,2).^2 +arrdif(:,3).^2 );
        filteridx2 = dis>0.4 ;

        %Select the portion with a z-value greater than 2 (for speed up), which also almost remove the ground and trunk points选择z值大于2的部分（加快速度） ,且大体去除了地面点和树干点
        filteridx1 = SIFdata_r_std(:,3)> (min(SIFdata_r_std(:,3))+2); %SIFdata_r_std(:,3)>2;

        % Remove the point facing backward of the sensor; 去除传感器后向的点
        if obs_theta<0
            filteridx3 = SIFdata_r_std(:,1)>sensor_obs(1);
        elseif obs_theta>0
            filteridx3 = SIFdata_r_std(:,1)<sensor_obs(1);
        else
            filteridx3 = filteridx1; % 相当于没做
        end
        % satisfied at the same time同时满足
        filteridx = filteridx1 & filteridx2 & filteridx3;
        %filteridx(filteridx(:)==0)=[];
        SIFdata_filter= SIFdata_r_std(filteridx,:);

        %%A beam of light is emitted from sensor point to sensor, and the distance to the point in the current differential cell is calculated, each beam retains 1 closest point 逐个传感器点发射一束光，计算当前微分单元中的点距离,每一束保留1个最近点
        shape =  size(SIFdata_filter);
        SIFdata_filter_copy = zeros(shape(1), shape(2)+1);
        SIFdata_filter_copy(:,1:4)=SIFdata_filter(:,1:4);


        x1=sind(obs_theta)*cosd(fi0); %fi0 =0??
        y1=sind(obs_theta)*sind(fi0);
        z1=cosd(obs_theta);
        %x1=-0.22-x0;y1=3.87-y0;z1=0.23-z0;
        %nv=sqrt(x1^2+y1^2+z1^2);
        x_obs =sensor_obs(1) ; y_obs=sensor_obs(2); z_obs=sensor_obs(3);
        x_all=SIFdata_filter_copy(:,1)-x_obs;
        y_all=SIFdata_filter_copy(:,2)-y_obs;
        z_all=SIFdata_filter_copy(:,3)-z_obs;
        
        theta_plus=0;
        theta=obs_theta+theta_plus;
        xd=sind(theta)*cosd(fi0);
        yd=sind(theta)*sind(fi0);
        zd=cosd(theta);
        basev=[x1,y1,z1];
        roatev=[xd,yd,zd];
        %xd=-x1;yd=-y1;zd=-z1;
        %for t=1:1:360
        t_angle=1*pi/180;
        h = makehgtform('axisrotate',basev,t_angle);
        v11 = h(1:3,1:3)*roatev(:);
        x_roate=v11(1);
        y_roate=v11(2);
        z_roate=v11(3);
        


        for i=1:length(SIFdata_filter_copy)
            xt=x_all(i);
            yt=y_all(i);
            zt=z_all(i);
            c=sqrt((y_roate*zt-yt*z_roate)^2+(z_roate*xt-zt*x_roate)^2+(x_roate*yt-xt*y_roate)^2);
            demo=sqrt(xd^2+yd^2+zd^2);
            distance=c/demo;
            SIFdata_filter_copy(i,5)=distance;
        end
        index_short=find(SIFdata_filter_copy(:,5)<0.01); % revised for calculated birdirection sif. the distance threshold will affect the value of canopy cover

        results=SIFdata_filter_copy(index_short,:);
        [m,maxnumber]=max(results(:,3));
        row=row+1;
        
        if isempty(maxnumber)
            find_index(row)=0;
            %sifpoints(row,:)=0;
        else
            index = index_short(maxnumber);
            find_index(row)= SIFdata_filter_copy(index,4); %Column 4 is the ID number; 第4列是编号
            %sifpoints(row,:)=results(maxnumber,:);
        end

        
        
        % sifpoints(:,1)=sifpoints(:,1)+x0;
        % sifpoints(:,2)=sifpoints(:,2)+y0;
        % sifpoints(:,3)=sifpoints(:,3)+z0;
        % sifpoints(sifpoints(:)==0)=[];
        

    end
    toc
    find_index(find_index(:)==0)=[];
    canopycover = length(find_index)/a(1);
    find_index=unique(find_index,'rows');

    %fullname = sprintf("%s\\leafindex_%s.txt",dirname,num2str(obs_theta))
    %save(fullname, 'find_index', '-ascii')

    datatest=SIFdata_r_std(find_index,:);
    leafdatasetName = sprintf("%s\\%s_leafdatatest_%s_%s.txt",dirname,filename,num2str(solar_azimuth),num2str(obs_theta));
    save(leafdatasetName, 'datatest', '-ascii')
    %SIFnadir=sum(sifpoints(:,4))/(length(sifpoints));

end
%% 







