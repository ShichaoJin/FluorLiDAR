clear all; clc;
a = 0;
disp("simulating each point")

%% load SIF data
workdir = "\\Jinlab309\d\JinShichao\3DSIF\Malixia_14Layer30Points";                                     % work directory
%filepath = strcat(workdir,'\SIF_Malixia\','raw','\');
filepath = strcat(workdir,'\MLX1107_3angles\','PAR2SIF','\');                                           % SIF data directory
Files=dir(['\\Jinlab309\d\JinShichao\3DSIF\Malixia_14Layer30Points\MLX1107_3angles\PAR2SIF\*.txt']);    % SIF files

%% directional sif simulation
for file=1:1 
    parfile = Files(file).name;
    filename = sprintf('%s\\%s',filepath,parfile)
    % ����ģ��õ�SIF 
    sifdistribution=load(filename); % (n,10); x,y,z,ֱ�������ҶSIF��ɢ�������ҶSIF�� Xxx
    sifdistribution = sifdistribution(:,1:6); 

    %% ���ҹ۲ⷶΧ�ڵĵ㣬���㷽����SIF
    % ģ����ƽ��ʹ�ֱ��ƽ��2��ƽ�棬 ÿ��ƽ�����19���춥��(-90��10:90)
    %������λ��
    x0=1.33;y0=3.98; % ����ƽ��λ��
    treetop = 3.15; %���������ƽ�� �����߶�
    r = 18;  %������������������CMOS�ľ���

    number = 38;  % ��̫������ƽ��1����ֱ̫������ƽ��2�� *��-90,90,10�� = 2*19=38
    sif_canopy=zeros(number/2,4); %sif, cc; sif2,cc2
    % sif_direct=zeros(number,1);
    % sif_ms=zeros(number,1);

    j=1;
    for azimuth = [0, 90] % ̫����ƽ��ʹ�ֱƽ��
        i=1
        for zenith = -90:5:90 %�۲��춥��
            % ���� bidirection����ȡ�ڲ�SIF
            [leafindex,canopycover] = bidirection_cylinder(parfile(1:end-4), sifdistribution(:,1:3), zenith, x0, y0, treetop, r, azimuth);
            sif = sifdistribution(leafindex,:);
            sif_canopy(i,j) = (mean(sif(:,4))+mean(sif(:,5)));
            sif_canopy(i,j+1) =canopycover;    
            i=i+1;
        end
        j = j+2;
    end

    filename_path=strcat('E:\NJU\3DSIF_Zch\Figure5-8CanopySIF_vs_SimulatedSIF_newData_totalPAR_reviseLeafSIF_RSErevise\','BiSIF_of_001_cylinder',parfile(1:end-4),'.csv');
    csvwrite(filename_path,sif_canopy)

    disp("great!")


end



    



