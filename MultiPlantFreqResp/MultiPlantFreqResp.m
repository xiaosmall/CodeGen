function    [ GM_Pos, GM_Neg, PM,  ValuedBatRequiredFreq,nErrCode]= MultiPlantFreqResp( ...
    vPlantFreq1000_1,vPlantGain1000_1,vPlantPhase1000_1,Lenght_1, ...
    vPlantFreq1000_2,vPlantGain1000_2,vPlantPhase1000_2,Lenght_2,...
    vPlantFreq1000_3,vPlantGain1000_3,vPlantPhase1000_3,Lenght_3,...
    vPlantFreq1000_4,vPlantGain1000_4,vPlantPhase1000_4,Lenght_4,...
    vPlantFreq1000_5,vPlantGain1000_5,vPlantPhase1000_5,Lenght_5,...
    vPlantFreq1000_6,vPlantGain1000_6,vPlantPhase1000_6,Lenght_6,...
    vPlantFreq1000_7,vPlantGain1000_7,vPlantPhase1000_7,Lenght_7,...
    vPlantFreq1000_8,vPlantGain1000_8,vPlantPhase1000_8,Lenght_8,...
    vPlantFreq1000_9,vPlantGain1000_9,vPlantPhase1000_9,Lenght_9, ...
    NumberOfPlants, ...
    PosFilt,PosFilt_Usage,PosGain, AccFFW, VelGain, VelKi,PIV_Usage, ...
    VelFilt1, VelFilt1_Usage, ...
    VelFilt2, VelFilt2_Usage, ...
    VelFilt3, VelFilt3_Usage, ...
    VelFilt4, VelFilt4_Usage, ...
    Ts, RequireLoopType,requiredFreq,GreaterLessSign)
%%Ver02, @20161204 ;Sophia, add nErrCode;
%%Ver03, @20161207 ;if freq has repeatable value,we remove the value,not
%%only feedback nErrCode;

% % input argument% % ====
%       vPlantGain1000 - Vector of gains in dB of the plant frequency response.
%       vPlantPhase1000 -Vector of phase in Deg of the plant frequency response.
%       ...
%       NumberOfPlants:  how many plants to calculate;
%       PosFilt,PosFilt_Usage,PosGain, AccFFW, VelGain, VelKi,PIV_Usage, :-
% %     VelFilt1, VelFilt1_Usage, ...
% %     VelFilt2, VelFilt2_Usage, ...
% %     VelFilt3, VelFilt3_Usage, ...
% %     VelFilt4, VelFilt4_Usage, ...%       Ts,
%       control parameter used;
%          RequireLoopType :   OL   1  %open loop
%                              CL   2  %close loop
%                              DT   3  %disturbance require loop
%       requiredFreq,:- Hz
%       GreaterLessSign:-   >, +1
%                           <, -1

% % output :% % ====
% % % % % % % %  GM_Pos:- dB;
% % % % % % % %  GM_Neg:- dB;
% % % % % % % %   PM:- deg
% % % % % % % % % ValuedBatRequiredFreq -: dB
%%nErrCode = -31;%used plant lenght is not positive value
% #codegen
%%init value of output for codegen
GM_Pos = 0;
GM_Neg = 0;
PM = 0;
ValuedBatRequiredFreq = 0;
nErrCode=0; %init no err;
GM_Pos01=zeros(1,NumberOfPlants);
GM_Neg01=zeros(1,NumberOfPlants);
PM01=zeros(1,NumberOfPlants);
fLg = zeros(1,NumberOfPlants);

para_structure = struct('Freq',zeros(1000,1), 'Gain', zeros(1000,1),...
    'Phase',zeros(1000,1),'Length',0);
para = repmat(para_structure,1,9);

%%%% calculate each plant magnitude & phase vector;
if VelFilt1_Usage ==0
    VelFilt1 = [0,0,0,0,0];
end
if VelFilt2_Usage ==0
    VelFilt2 = [0,0,0,0,0];
end
if VelFilt3_Usage ==0
    VelFilt3 = [0,0,0,0,0];
end
if VelFilt4_Usage ==0
    VelFilt4 = [0,0,0,0,0];
end
if PosFilt_Usage ==0
    PosFilt = [0,0,0,0,0];
end
if PIV_Usage==0  % PIV_Usage
    PosGain=0;
    AccFFW=0; VelGain=0; VelKi=0;
    
end
if RequireLoopType == 3 %%disturbance
    LoopType = 4;
elseif RequireLoopType == 1 %%openloop% with Pos&Vel OpenLoop
    LoopType = 1;
else
    LoopType = 2; %%pos+vel close loop
end
%%if the length of the freq and gain ,phase is the same, then
para(1).Freq = vPlantFreq1000_1;
para(1).Gain = vPlantGain1000_1;
para(1).Phase = vPlantPhase1000_1;
para(1).Length = Lenght_1;

para(2).Freq = vPlantFreq1000_2;
para(2).Gain = vPlantGain1000_2;
para(2).Phase = vPlantPhase1000_2;
para(2).Length = Lenght_2;

para(3).Freq = vPlantFreq1000_3;
para(3).Gain = vPlantGain1000_3;
para(3).Phase = vPlantPhase1000_3;
para(3).Length = Lenght_3;

para(4).Freq = vPlantFreq1000_4;
para(4).Gain = vPlantGain1000_4;
para(4).Phase = vPlantPhase1000_4;
para(4).Length = Lenght_4;

para(5).Freq = vPlantFreq1000_5;
para(5).Gain = vPlantGain1000_5;
para(5).Phase = vPlantPhase1000_5;
para(5).Length = Lenght_5;

para(6).Freq = vPlantFreq1000_6;
para(6).Gain = vPlantGain1000_6;
para(6).Phase = vPlantPhase1000_6;
para(6).Length = Lenght_6;

para(7).Freq = vPlantFreq1000_7;
para(7).Gain = vPlantGain1000_7;
para(7).Phase = vPlantPhase1000_7;
para(7).Length = Lenght_7;

para(8).Freq = vPlantFreq1000_8;
para(8).Gain = vPlantGain1000_8;
para(8).Phase = vPlantPhase1000_8;
para(8).Length = Lenght_8;

para(9).Freq = vPlantFreq1000_9;
para(9).Gain = vPlantGain1000_9;
para(9).Phase = vPlantPhase1000_9;
para(9).Length = Lenght_9;
fvfreq = zeros(1000,1);
for ii =1 :NumberOfPlants
    [vLg01, vLp, GM_Pos01(ii), GM_Neg01(ii), PM01(ii), ~, ~, ~, ...
        PeakValueFreq, DisturbancePeakValue01, DisturbancePeakValueFreq, ...
        ~, ~, ~] ...
        =CalculateAnyFrequencyResponse( para(ii).Freq,para(ii).Gain,para(ii).Phase,para(ii).Length, PosFilt, PosGain, AccFFW, VelGain, VelKi, ...
        VelFilt1, VelFilt2, VelFilt3, VelFilt4, Ts, LoopType);
    nErrCode = inputPositiveCheck(para(ii).Length,-31);
    if nErrCode<0
        return;
    end

    fvfreq = para(ii).Freq(1:para(ii).Length);
    %%after found un-unique,chg to solve f is non-unique problem;    
    [f_afterUnique, index_f_afterUnique] = unique(fvfreq); %sort and remove duplicates
    vLg01_afterUnique = vLg01(index_f_afterUnique);
    nLength = length(index_f_afterUnique);
    fLg(ii) = interp1(f_afterUnique(1:nLength),vLg01_afterUnique(1:nLength),requiredFreq);
end
if  GreaterLessSign >0
    ValuedBatRequiredFreq = min(fLg(1:NumberOfPlants));
else
    ValuedBatRequiredFreq = max(fLg(1:NumberOfPlants));
end
GM_Pos = min(GM_Pos01(1:NumberOfPlants)); %worst of GM_Pos
GM_Neg = min(GM_Neg01(1:NumberOfPlants)); %worst of GM_Neg

PM = min(PM01(1:NumberOfPlants));

end