% function run2test
close all
clc
%%
plantDataType = 'xlsFormat';
switch 'Original'
    case 'Original'
        plant = Agito2Plant('Linear Stage Far Resonance Set 1.xls');%MeasuredFreqs
        f=plant.f(:);
        pg=20*log10(abs(plant.r(:)));
        pp=angle(plant.r(:))*180/pi;
        n = length(f);
        %%plant02
        plant2 = Agito2Plant('Linear Stage Far Resonance Set 2.xls');
        f2=plant2.f(:);
        pg2=20*log10(abs(plant2.r(:)));
        pp2=angle(plant2.r(:))*180/pi;
        n2 = length(f2);
    case 'xlsFormat'
        plant = xlsread('OKCase1.xls');
        pg= plant(:,1);
        pp=plant(:,2);
        f=plant(:,3);
        n = length(f);
        FigOn=1;
        if FigOn
            figure;
            subplot(211);semilogx(f,pg);
            subplot(212);semilogx(f,pp);
        end
        f2 = f;
        pg2=pg;
        pp2 =pp;
        n2 = n;
end

kpp = 112%42;
kvp = 1709%1848;
kvi = 42;
vPosFilt = [0 0 0 0 0];

VelFilt1 = [2,485,80,0,0];
VelFilt2 = [8,384,20,50,0];
VelFilt3 = [2,30,0.5,0,0];
VelFilt4 = [2,3000,0.5,0,0];
Ts = 1/16384;
AccFFW =0;
VelFilt1_Usage =1; %%sophia chg from value to one variable;
VelFilt2_Usage =1;
VelFilt3_Usage =0;
VelFilt4_Usage =0;
PosFilt_Usage = 0; %%not optimizae the filter
PIV_Usage =1;

LoopType= 1;% with Pos&Vel OpenLoop, disturbance;
NumberOfPlants = 3;
requiredFreq = 20; %Hz;
GreaterLessFlag =1  %%>
f3  = zeros(1,1000);
pg3 = zeros(1,1000);
pp3 = zeros(1,1000);
n3  = 0;
[ GM_Pos, GM_Neg, PM,  PeakValue,nErrCode]= MultiPlantFreqResp( ...
    f, pg,pp,n,f2,pg2,pp2,n2,f3,pg3,pp3,n3,f,pg,pp,n,f,pg,pp,n,f,pg,pp,n,f,pg,pp,n,f,pg,pp,n,0,0,0,0,...
    NumberOfPlants, ...
    vPosFilt,PosFilt_Usage,kpp, AccFFW, kvp, kvi,PIV_Usage, ...
    VelFilt1, VelFilt1_Usage, ...
    VelFilt2, VelFilt2_Usage, ...
    VelFilt3, VelFilt3_Usage, ...
    VelFilt4, VelFilt4_Usage, ...
    Ts, LoopType,requiredFreq,GreaterLessFlag)
if nErrCode <0
    disp('err in MultiPlantFreqResp');
else
    nErrCode
end
Fig2Show=1;
if (Fig2Show)
    close all;
    if NumberOfPlants >=1
    figure(10);
    subplot(211);semilogx(f,pg);hold on;
    subplot(212);semilogx(f,pp);hold on;
    legend('Phase');
    end
    if NumberOfPlants>=2
    figure(10);
    subplot(211);semilogx(f2,pg2,'r');
    subplot(212);semilogx(f2,pp2,'r');legend('Phase');title('Plant02');
    end
    
end

PeakValue
% end
