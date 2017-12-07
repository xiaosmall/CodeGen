clear all
close all
clc
% SetLSpath
%%
%%plant
savePlanton =0;
if savePlanton
Freq = [10.00,15.00,20.00,25.00,30.00,35.00,40.00,50.00,60.00,70.00,80.00,90.00,100.00,120.00,150.00,200.00,250.00,300.00,400.00,500.00]
Gain = [10.48,-6.03,-21.78,-29.38,-25.18,-29.03,-23.31,-24.01,-26.18,-27.91,-29.07,-31.60,-31.76,-34.43,-37.13,-42.46,-46.01,-49.04,-54.11,-57.78]
MinFreqRange = 29.66;
MaxFreqRange = 136.82;
Freq = Freq';
Gain = Gain';
end
if savePlanton ==0
    
mdlNo =1
switch mdlNo 
 case 1
          FileName = 'Linear Stage Far Resonance Set 1.xls';
        plant = Agito2Plant(FileName); %%raw xls data in dB, plant data nondB,plant f inHz
          
        
case 10
   
      Filename = 'MeasuredFreqs.xls'%'MeasuredFreqs_4mat2011.xls';
       plant = Agito2Plant( Filename );
  case 2
  
    load PaperExample02Plantv14.mat;
    Filename ='';
    ff = plant.Frequency;
    frequency = ff/2/pi; %%this is for the data with rad/s unit
    Mag = plant.ResponseData;
    plant = frd(Mag,frequency,'unit','Hz');
 case 3
          FileName = 'Linear Stage Far Resonance Set 2.xls';
        plant = Agito2Plant(FileName); %%raw xls data in dB, plant data nondB,plant f inHz
  case 4
        addpath('D:\07_Ctr\10_plantData');
        FileName = 'Plant_withLF_Res_realcase.xls';
        plant = Agito2Plant(FileName); %%raw xls data in dB, plant data nondB,plant f inHz
end
end
f=plant.f(:);
rdb=20*log10(abs(plant.r(:)));
pp=angle(plant.r(:))*180/pi;
n = length(f);%%pay attention that the n muse be the min and max  points of no

plotflag =1;
if plotflag
    figure(300);
    subplot(211);semilogx(f,rdb);grid on;ylabel('Mag,dB');hold on;xlabel('Freq,Hz')

 end


MinFreqRange =   1.7567      %8.9%5 %2%10;
MaxFreqRange =343.3811      %44%150 %90;
[SlopeResult, PlantGain] = LinearSlope4Manual(f, rdb, n, MinFreqRange, MaxFreqRange)
% %%later uncomment it 
OnlyTest_LinearSlope4Manual =0;
if OnlyTest_LinearSlope4Manual
Slope = -40;
Tolerance = 0.02
[StartFrequency, EndFrequency, SlopeResult, NumberOfPoints,PlantGain] = PGain_from_LinearSlope(f, rdb, n, MinFreqRange, MaxFreqRange, Slope, Tolerance)

s= tf('s');
H=PlantGain/s^2;		

        
      freq=[5:700];%[10:500];%%rad/s
      plantIDEN = frd(H,freq,'unit','Hz');%%AutoTuning required rad/s
      rr = 20*log10(abs(squeeze(plantIDEN.ResponseData)));
   
%      H22=PlantGain22/s^2;		
%      plantIDEN22 = frd(H22,freq,'unit','Hz');%%AutoTuning required rad/s
%      rr22 = 20*log10(abs(squeeze(plantIDEN22.ResponseData)));
%      subplot(211);semilogx(plantIDEN22.Frequency,rr22,'g');grid on;xlabel('Mag,dB');hold on;

if plotflag
    figure(300);hold on;
    subplot(211);semilogx(plantIDEN.Frequency,rr,'r');grid on;xlabel('Mag,dB');hold on;
end

end
% hplt = bodeplot(H,'r.');grid on; setoptions(hplt,'FreqUnits','Hz','PhaseVisible','on');
%         legend('Measured','ParamatricMdlH');
%         xlabel('Freq,Hz');
%         rr = (abs(squeeze(plantIDEN.ResponseData)));
% %%%original piv para cal
%%    P = Agito2Plant(FileName);
%    
% 	rdb = 20*log10(abs(P.r));
% 	f = P.f;
% 	p = angle(r)*180/pi;
% 	n = length(r);
%     subplot(212);semilogx(f,pp);grid on;xlabel('Freq,Hz')
%     bodeplot(plant);grid on;