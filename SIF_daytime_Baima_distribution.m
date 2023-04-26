clear all;clc;
filepath_SIF='E:\NJU\3DSIF_Zch\PAR2SIF\'; %SIF输出点云保存路径
Files=dir(['E:\NJU\3DSIF_Zch\PAR8_18\*.txt']); %辐射输入点云所在目录
filepath='E:\NJU\3DSIF_Zch\PAR8_18\';
RT_ms=load('RT_ms3.txt');  % the RT_ms3.txt file was derived by 4-scale model
RZT_ms=load('RZT_ms2.txt'); 
number=length(Files); 
%sifnadir=zeros(number,1);  % SIF simulation at nadir direction

%获取V
X                           = textread('../inputdata.txt','%s');
N                           = str2double(X);
V                           = assignvarnames; % 
options.Cca_function_of_Cab = 0;

for i = 1:length(V)  % assign values to variables using string matching
    j = find(strcmp(strtok(X(:,1)),V(i).Name));
    if isempty(j) || isnan(N(j+1))                       % N的值既不是nan,J的值也不为空
        if i==2
            fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input spreadsheet...');
            fprintf(1,'%s %s %s\n', 'I will use 0.25*Cab instead');
            options.Cca_function_of_Cab = 1;
        else            
            if ~(options.simulation==1) && (i==30 || i==32)
                fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input spreadsheet...');
                fprintf(1,'%s %s %s\n', 'I will use the MODTRAN spectrum as it is');
            else
                if (options.simulation == 1 || (options.simulation~=1 && (i<46 || i>50)))
                    fprintf(1,'%s %s %s \n','warning: input "', V(i).Name, '" not provided in input spreadsheet');
                    if (options.simulation ==1 && (i==1 ||i==9||i==22||i==23||i==54 || (i>29 && i<37)))
                        fprintf(1,'%s %s %s\n', 'I will look for the values in Dataset Directory "',char(F(5).FileName),'"');
                    else
                        if (i== 24 || i==25)
                            fprintf(1,'%s %s %s\n', 'will estimate it from LAI, CR, CD1, Psicor, and CSSOIL');
                            options.calc_zo = 1;
                        else
                            if (i>38 && i<44)
                                fprintf(1,'%s %s %s\n', 'will use the provided zo and d');
                                options.calc_zo = 0;
                            else
                                if ~(options.simulation ==1 && (i==30 ||i==32))
                                    fprintf(1,'%s \n', 'this input is required: SCOPE ends');
                                    return
                                else
                                    fprintf(1,'%s %s %s\n', '... no problem, I will find it in Dataset Directory "',char(F(5).FileName), '"');
                                end
                            end
                        end
                    end
                end
            end
        end
    end
    
    j2 = []; j1 = j+1;
    while 1
        if isnan(N(j1)), break, end
        j2 = [j2; j1]; %#ok<AGROW>
        j1 = j1+1;
    end
    if isempty(j2)
        V(i).Val            = -999;
    else
        V(i).Val            = N(j2);
    end
end
%获得leafbio
vi = ones(length(V),1);
leafbio.Cab = 40;   %叶绿素浓度
leafbio.Cca = V(2).Val(vi(2));
leafbio.Cdm  = V(3).Val(vi(3));
leafbio.Cw = V(4).Val(vi(4));
leafbio.Cs   = V(5).Val(vi(5));
leafbio.N = V(6).Val(vi(6));
leafbio.Vcmo  = V(9).Val(vi(9));
leafbio.m  = V(10).Val(vi(10));
leafbio.Type = V(11).Val(vi(11));
leafbio.Tparam  = V(14).Val(:); % this is correct (: instead of 14)
leafbio.fqe(2)   = V(15).Val(vi(15));
leafbio.fqe(1)   = leafbio.fqe(2)/5;
leafbio.Rdparam  = V(13).Val(vi(13));

leafbio.rho_thermal   = V(7).Val(vi(7));
leafbio.tau_thermal  = V(8).Val(vi(8));

leafbio.Tyear = V(55).Val(vi(55));
leafbio.beta   = V(56).Val(vi(56));
leafbio.kNPQs   = V(57).Val(vi(57));
leafbio.qLs  = V(58).Val(vi(58));
leafbio.stressfactor  = V(59).Val(vi(59));
%获取optipar
path_input      = '../../data/input/';          % path of all inputs
opticoef    = load([path_input,'fluspect_parameters/Optipar_fluspect_2014.txt']);  % file with leaf spectral parameters
% Optical coefficient data used by fluspect
optipar.nr    = opticoef(:,2);
optipar.Kab   = opticoef(:,3);
optipar.Kca   = opticoef(:,4);
optipar.Ks    = opticoef(:,5);
optipar.Kw    = opticoef(:,6);
optipar.Kdm   = opticoef(:,7);
optipar.phiI  = opticoef(:,9);
optipar.phiII = opticoef(:,10);


%计算吸收率
[spectral] = define_bands;
%wlS  = spectral.wlS;    % SCOPE 1.40 definition
%wlP  = spectral.wlP;    % PROSPECT (fluspect) range
%wlT  = spectral.wlT;    % Thermal range
%wlF  = spectral.wlF;    % Fluorescence range
%nwlP = length(wlP);
%nwlT = length(wlT);
%[rho,tau,rs] = deal(zeros(nwlP + nwlT,1));
%IwlP     = spectral.IwlP;
%IwlT     = spectral.IwlT;
%rho(IwlP)  = leafopt.refl;
%tau(IwlP)  = leafopt.tran;
fversion = @fluspect_bcar;
[leafopt] = fversion(spectral,leafbio,optipar);
rho     = leafopt.refl;
tau     = leafopt.tran;
abs     =1-rho-tau;  %absorption

VISresolution   = 00.001;                       %[um]               Resolution of spectrum in VISible
NIRresolution   = 00.010;                       %[um]               Resolution of spectrum in Near InfraRed
MWIRresolution  = 05.000;                       %[um]               Resolution of spectrum in Mid Wave InfraRed
TIRresolution   = 30.000;                       %[um]               Resolution of spectrum in Thermal InfraRed
VIS         =  0.40:0.001         : 0.85;
NIR             =  0.85:NIRresolution : 2.50;   % [um]              Near Infrared 
MWIR            =  2.50:MWIRresolution:10.00;   % [um]              Spectrum Mid Wave Infra Red
TIR             = 10.00:TIRresolution: 50.00;   % [um]              Spectrum Thermal  Infra Red
wl              = [VIS,NIR, MWIR, TIR]';        % [um]              Spectrum of wavelenghts 
wl              = unique(wl);
nwl             = size(wl       ,1); 
resolution      = [wl(2:end)-wl(1:end-1);wl(nwl)- wl(nwl-1)];

rad         = load('rad.txt');                  %    MODTRAN output 
wavel           = rad(:,1);                     %[nm]               %wavelength in Modtran file
wavel           = wavel/1000;                   %[um]               %wavelength in Modtran file
Esun_M          = rad(:,2);                     %[W m-2 um-1]       %Spectral Solar Radiation for Modtran wavelengths
Esky_M          = rad(:,3);                     %[W m-2 um-1]       %Spectral Sky   Radiation for Modtran wavelengths
[wavel,I]       = sort(wavel);                  %[um]               %reverse (from low to high lambda)
Esky_M          = Esky_M(I);                    %[W m-2 um-1]       %reverse (from low to high lambda)
Esun_M          = Esun_M(I);                    %[W m-2 um-1]       %reverse (from low to high lambda)
Esky_M          = interp1(wavel,Esky_M,wl);     %[W m-2 um-1]       %Spectral Solar Radiation for model wavelengths
Esun_M          = interp1(wavel,Esun_M,wl);     %[W m-2 um-1]       %Spectral Solar Radiation for modelwavelengths
% calculate a normalised spectrum of Esun and Esky

[fEsun ,...
 fEsky ,...
 fEskyt]        = deal(zeros(size(wl)));

%J_o             = wl<0.701;                                         %find optical spectrum 400-700nm
J_o             = wl<0.680;                                         %find optical spectrum 400-680nm
Esunt           = Esun_M(J_o)'*resolution(J_o);                     %Calculate optical sun fluxes (by Integration)
Eskyt           = Esky_M(J_o)'*resolution(J_o);                     %Calculate optical sun fluxes (by Integration)

%fEsun(J_o)      = Esun_M(J_o)/Esunt;                                 %fraction of contribution of Sun fluxes to total light
%fEsky(J_o)      = Esky_M(J_o)/Eskyt;                                 %fraction of contribution of Sky fluxes to total light
Etot= Esunt+Eskyt; %Shichao Added
fEsun(J_o)      = Esun_M(J_o)/Etot;                                 %fraction of contribution of Sun fluxes to total light
fEsky(J_o)      = Esky_M(J_o)/Etot;                                 %fraction of contribution of Sky fluxes to total light
%荧光波段
wlS          = spectral.wlS';       % SCOPE wavelengths, make column vectors
wlF          = spectral.wlF';       % Fluorescence wavelengths
wlE          = spectral.wlE';       % Excitation wavelengths
[dummy,iwlfi]    = intersect(wlS,wlE); %#ok<ASGLU>
[dummy,iwlfo]    = intersect(wlS,wlF);
MbI          = leafopt.MbI;
MbII         = leafopt.MbII;
MfI          = leafopt.MfI;
MfII         = leafopt.MfII;
MpluI                = 0.5*(MbI+MfI);              % [nwlfo,nwlfi]             the factor 1000 is needed for conversion of W m-2 sr-1 nm-1 to W m-2 sr-1 um-1
MminI                = 0.5*(MbI-MfI);              % [nwlfo,nwlfi]
MpluII                = 0.5*(MbII+MfII);              % [nwlfo,nwlfi]             the factor 1000 is needed for conversion of W m-2 sr-1 nm-1 to W m-2 sr-1 um-1
MminII                = 0.5*(MbII-MfII); 

for file=1:number

%filename_path=strcat(filepath_SIF,num2str(daytime));  %保存SIFtxt的路径和命名

filename=Files(file).name;

endnum = length(filename);
time = filename(endnum-5:endnum-4);
hour = str2num(time);
daytime = (hour-8)*60 + 1 ;
filename_path=strcat(filepath_SIF,'SIF_',filename(5:length(filename)-4));
LiDARdata=load(strcat(filepath,filename));  %导入该时刻辐射数据, speed of this step can be optimized. if we organize the simulated PAR at different times in a same liDAR, then we can read a same LiDAR data and select the right columes here

%计算APAR
PAR=LiDARdata(:,6);
PAR_sun = LiDARdata(:,4); % 此处的PAR是最终正确值
PAR_sky = LiDARdata(:,5);

sifdistribution=zeros(length(PAR),10);
sifdistribution(:,1:3)=LiDARdata(:,1:3);  %获取点云坐标
SIF_MS=zeros(length(PAR),1);
%lambda=400:1:700;   %PAR波段范围
lambda=400:1:680; 

for i=1:length(PAR)  

    leafsun=fEsun(J_o) .*PAR(i);
    leafdiff=fEsky(J_o).*PAR(i);
    %leafsun=fEsun(J_o) .*PAR_sun(i); %Shichao, the result will be small
    
    %leafdiff=fEsky(J_o).*PAR_sky(i); %Shichao
    
    Esunf_=leafsun.*abs(J_o);
    Epluf_=leafdiff.*abs(J_o);
    Esunf_=(phot2e(lambda,Esunf_'))'; %umol/m2/nm转w/m2
    Epluf_=(phot2e(lambda,Epluf_'))'; %umol/m2/nm转w/m2

    MpluEsunI  = MpluI * Esunf_;  %MpluI 为光系统1的后向发射矩阵
    MminEsunI  = MminI * Esunf_;
    MpluEsunII  = MpluII * Esunf_; 
    MminEsunII  = MminII * Esunf_;
    
    MpluEpluI  = MpluI * Epluf_;
    MminEpluI  = MminI * Epluf_;	
    MpluEpluII  = MpluII * Epluf_;
    MminEpluII  = MminII * Epluf_;	
	
	sifdistribution(i,4)=(MpluEsunI(98)+MpluEsunII(98)+MpluEpluI(98)+MpluEpluII(98))/pi+ (MminEsunI(98)+MminEsunII(98)+MminEpluI(98)+MminEpluII(98))/pi;%98对应737nm SIF发射  %640-850; MpluEsunI(92)最大，对应640+91=731nm
    SIF_MS(i)=RT_ms(daytime)*(MpluEsunI(98)+MpluEsunII(98)+MminEsunI(98)+MminEsunII(98))/pi+RZT_ms(daytime)*(MpluEpluI(98)+MpluEpluII(98)+MminEpluI(98)+MminEpluII(98))/pi;  %SIF多次散射
    
    sifdistribution(i,4) = sifdistribution(i,4);             % leaf-level SIF induced by direct PAR and diffuse PAR
    sifdistribution(i,5) = SIF_MS(i);                        % multi-scattered SIF 
    sifdistribution(i,6) = sifdistribution(i,4)+ SIF_MS(i);  % Total SIF (used for direct SIF calculation)
    
    sifdistribution(i,7)=(MpluEsunI(98)+MpluEsunII(98))/pi+(MminEsunI(98)+MminEsunII(98))/pi;%760nm SIF sunLeaf SIF   (exclude msSIF)
    sifdistribution(i,8)=(MpluEpluI(98)+MpluEpluII(98))/pi+ (MminEpluI(98)+MminEpluII(98))/pi;%760nm SIF shadeLeaf SIF (exclude msSIF)
      
    sifdistribution(i,9)= sifdistribution(i,7) + RT_ms(daytime)*(MpluEsunI(98)+MpluEsunII(98)+MminEsunI(98)+MminEsunII(98))/pi ;%760nm SIF sunLeaf SIF      (include msSIF)
    sifdistribution(i,10)=sifdistribution(i,8)  + RZT_ms(daytime)*(MpluEpluI(98)+MpluEpluII(98)+MminEpluI(98)+MminEpluII(98))/pi;%760nm SIF shadeLeaf SIF  (include msSIF)
    
end
eval(['save ' filename_path '.txt sifdistribution -ASCII']) %保存SIF点云


end
