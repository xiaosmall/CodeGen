%History: add pos filter to pos OL; July08-2015; 
%         chg to two GM output 
%V03  @20161129, if not found CrossOver Freq,CrossOver = NaNRepresent(
%999999) previous is 999999/2/pi£»
function [vOutGain1000, vOutPhase1000, GM_Pos,GM_Neg, PM, CrossOver, dB3Freq, PeakValue, PeakValueFreq, DisturbancePeakValue, DisturbancePeakValueFreq, ...
    NoiseHighestFreqGainValue, NoiseHighestFreqValue, Noise3dBFreqValue] = CalculateAnyFrequencyResponse( vPlantFreq1000, ...
    vPlantGain1000, vPlantPhase1000, ActualLength, vPosFilt1, PosGain, AccFFW, VelGain, VelKi, ...
    vVelFilt1, vVelFilt2, vVelFilt3, vVelFilt4, Ts, Type)

% [vOutGain1000, vOutPhase1000, GM, PM, CrossOver, dB3Freq, PeakValue, PeakValueFreq, DisturbancePeakValue, DisturbancePeakValueFreq, ...
%    NoiseHighestFreqGainValue, NoiseHighestFreqValue, Noise3dBFreqValue] = CalculateAnyFrequencyResponse( vPlantFreq1000, ...
%    vPlantGain1000, vPlantPhase1000, ActualLength, vPosFilt1, PosGain, AccFFW, VelGain, VelKi, ...
%    vVelFilt1, vVelFilt2, vVelFilt3, vVelFilt4, Ts, Type)
%
%   Calculate any frequency response of AGITO control loop.
%
%   Inputs:
%
%       vPlantFreq1000 - Vector of the frequencies that the Plant frequency
%           response was measured. Fixed size 1000 (code is generated for this
%           fix size). Use trailing zeros.
%       vPlantGain1000 - Vector of gains in dB of the plant frequency response.
%           As above.
%       vPlantPhase1000 -Vector of phase in Deg of the plant frequency response.
%           As Above.
%       ActualLength - Actual length of valid values in above input
%           vectors. output vectors also include valid values up to this
%           index.
%       vPosFilt1 - Vector that define Agito position Prefilter
%           All filer vectors are of fixed length 5.
%           (1) is a number that defines the filter, see below.
%           (2)-(5) are the filter parameters.
%       PosGain - Position feedback gain.
%       AccFFW - Acceleration feed forward gain.
%       VelGain - Velocity feedback proportional gain.
%       VelKi - Velocity feedback integration gain.
%       vVelFilt1 - Velocity backword feedback filter type and parameters.
%           See above for length and content.
%       vVelFilt2 - Velocity forward feedback filter.
%           See above for length and content.
%       vVelFilt3 - Velocity forward feedback filter.
%           See above for length and content.
%       vVelFilt4 - Velocity forward feedback filter.
%           See above for length and content.
%       Ts - Sampling time, in seconds
%       Type - Type of required loop output frequency response:
%%     this type is original type which is different from the  type inside of the code
%           Type 0: Plant.
%               Shows the plant as it was measured. CurrRef to Pos.
%
%           Type 1: Velocity Control Only: Open Loop.
%               Shows the open loop of the velocity control, assuming no position control at all.
%               Open loop is calculated by "opening the loop" at the
%               current command (CurrRef), which is the input to the Plant.
%
%           Type 2: Position Over Velocity: Closed Loop.
%               Shows the closed loop of the position control. From PosRef to
%               Pos. Of course, including the built in velocity FFW.
%               AccFFW parameter is not included in this calculation.
%
%           Type 3: Position Over Velocity: Closed Loop with Anti Vibrations Pre-Filter.
%               Shows the closed loop of the position control. Including the
%               pre-filer on the position reference. Of course, including the built in velocity FFW.
%               AccFFW parameter is not included in this calculation.
%
%           Type 4: Position Over Velocity: Input Disturbances.
%               Shows the response to input disturbances. Input disturbances
%               are at the input to the Plant. Calculated from disturbances
%               to Position.
%
%           Type 5: Position Over Velocity: Measurement Noises.
%               Shows the response to measurement noises. Measurement noises
%               are on the position reading. Calculated from noises to Position.
%
%           Type 6: Position Over Velocity: Open Loop.
%               Shows the true open loop of the overall system.
%               Open loop is calculated by "opening the loop" at the current
%               command (CurrRef), which is the input to the Plant.
%
%           Type 7: Velocity Control Only: Closed Loop.
%               Shows the closed loop of the velocity control, assuming no
%               position control at all. VelRef to Vel[2] (Vel[2] is the output of the derivative).
%
%           Type 8: Velocity Control Only: Closed Loop, to Filtered Velocity.
%               Shows the closed loop of the velocity control, assuming no position control at all. VelRef to Vel[1] (Vel[1] is the
%               output the bi-quad filter on the velocity reading).
%
%           Type 9: Velocity Control Only: Input Disturbances.
%               Shows the response to input disturbances. Input
%               disturbances are at the input to the Plant. Calculated from disturbances
%               to Vel[2] (derivative of position).
%
%           Type 10: Velocity Control Only: Measurement Noises.
%               Shows the response to measurement noises. Measurement noises
%               are on the position reading.. Calculated from noises to Vel[2]
%               (derivative of position).
%
%   Supported types of filters:
%
%           0           1           2               3           4           5                   6                        7             8                 9
%   ['NONE', 'LPF1', 'LPF2', 'LPF3', 'LDLG2', 'LDLG1', 'LDLG1FP', 'LDLG2FP', 'NOTCH', 'CLDLG'];
%
%   For the parameters for each type of filter, see Help for Nesection2coeff fucntion.
%
%   Outputs:
%
%       vOutGain1000 - Vector of gains in dB that represents the gain of the output frequency response
%           Take only first Actuallength elements. All others are zeros.
%       vOutPhase1000 - Vector of phases in Deg that represents the phase of
%           As above.
%       the output frequency response
%       A set of values with results of relevant calculations, as follows:
%
%       (Calcualtions refer to the Velocity control only of the requested
%       transfer function is for the velocity control only)
%
%           GM
%           PM
%           CrossOver
%           dB3Freq
%           PeakValue
%           PeakValueFreq
%           DisturbancePeakValue
%           DisturbancePeakValueFreq
%           NoiseHighestFreqGainValue
%           NoiseHighestFreqValue
%           Noise3dBFreqValue
%
%   Note:
%
%   From the matlab point of vew, input frequency/gain/[hase vector can be
%   of any length. However, when generating the C code, we use definition
%   fo vector of 1000 points, so that thw code will be generated for vector
%   of size 1000 for the frequency/gain/phase and 5 for each filter
%   definition.
%
%   Example:
%
%   plant = ...
%   Agito2Plant('LinearStageFarResonanceSingleSet.xls');f=plant.f(:);pg=20*log10(abs(plant.r(:)));pp=angle(plant.r(:))*180/pi;len=length(f);
%   [g, p, GM, PM, CrossOver, dB3Freq, PeakValue, PeakValueFreq, DisturbancePeakValue, ...
%       DisturbancePeakValueFreq, NoiseHighestFreqGainValue, NoiseHighestFreqValue, Noise3dBFreqValue] = ...
%       CalculateAnyFrequencyResponse(f,pg, ...
%       pp,len,[0,0,0,0,0],100,0,1000,20,[2,300,0.5,0,0], ...
%       [8,350,20,10,0],[0,0,0,0,0],[0,0,0,0,0],1/16384,1);
%
%#codegen
%
%   Take only valid vector parts
NaNRepresent =  999999;
vPlantFreq = vPlantFreq1000(1:ActualLength);
vPlantGain = vPlantGain1000(1:ActualLength);
vPlantPhase = vPlantPhase1000(1:ActualLength);
%

%
%   The function was written for transfer function types that are defined
%   diffrently than those defined now. The following is used to change the
%   type definition, so it will much the fucntion code.
%
%   The function was written for this type definitions (in contrast to what
%   defined above in the function header
%
%   0 ?Plant
%   1 ?Open Loop, Vel. control only
%   2 ?Closed Loop, Vel. control only
%   3 ?Closed Loop to Vel[1], Vel. control only
%   4 ?Input Disturbances, Vel. control only
%   5 ?Measurement Noises, Vel. control only
%   6 ?Open Loop
%   7 ?Closed Loop
%   8 ?Closed Loop with Pos. Filter 1
%   9 ?Input Disturbances Response
%   10 ?Measurement Noises Response
%
OriginalTypes = [0 6 7 8 9 10 1 2 3 4 5];
Type = OriginalTypes(Type+1);
%
% Calculate the plant frequency response and put it in FRD variable P
P = 10.^(vPlantGain/20).*exp(1i*vPlantPhase*pi/180) + (AccFFW - AccFFW); % +/-AccFFW just for Matlab to make the code fully valid: need to use all input arguments
s = vPlantFreq*2*pi*1i;
Lp = 1;
Lv = 1;
D = 1;

if Type==0
    % Case 0 => Plant
    G = P;
   
else
    [AgitoFilt01,N1] = Newsection2coeff(vVelFilt1(1), vVelFilt1(2), vVelFilt1(3), vVelFilt1(4), vVelFilt1(5), 1/Ts);
    [AgitoFilt02,N2] = Newsection2coeff(vVelFilt2(1), vVelFilt2(2), vVelFilt2(3), vVelFilt2(4), vVelFilt2(5), 1/Ts);
    [AgitoFilt03,N3] = Newsection2coeff(vVelFilt3(1), vVelFilt3(2), vVelFilt3(3), vVelFilt3(4), vVelFilt3(5), 1/Ts);
    [AgitoFilt04,N4] = Newsection2coeff(vVelFilt4(1), vVelFilt4(2), vVelFilt4(3), vVelFilt4(4), vVelFilt4(5), 1/Ts);
    VF01 = AgitoFiltFreqs(AgitoFilt01,vPlantFreq*2*pi,1/Ts,N1);
    VF02 = AgitoFiltFreqs(AgitoFilt02,vPlantFreq*2*pi,1/Ts,N2);
    VF03 = AgitoFiltFreqs(AgitoFilt03,vPlantFreq*2*pi,1/Ts,N3);
    VF04 = AgitoFiltFreqs(AgitoFilt04,vPlantFreq*2*pi,1/Ts,N4);

    Cv = VelGain./2^16*(VelKi./s+1);
    D = s./(s*Ts/2+1);
    Lv = Cv.*VF01.*VF02.*VF03.*VF04.*P.*D;
    if Type == 1,
        % Case 1 => Open velocity loop  Cv*VF2*VF3*VF4*P*D*VF1
        
        G = Lv;
    else
        if Type == 2,
            % Case 2 => veloty command to unfilered velocity  Lv/(1+Lv)
            
            G=Cv.*VF02.*VF03.*VF04.*P.*D./(1+Lv);
        elseif Type == 3,
            % Case 3 => Veloty command to velocity feedabck Lv/(1+Lv)
            G = Lv./(1+Lv);
        elseif Type == 4,
            % Case 4 => Input disturbance to veloty PD/(1+Lv)
            G = P.*D./(1+Lv);
        elseif Type == 5,
            % Case 5 => Noise to veloty D/(1+Lv)
            G = D./(1+Lv);
        else
            Lp = (D+PosGain).*Cv.*VF01.*VF02.*VF03.*VF04.*P;
            if Type == 6,
                % Case 6 => Posiotn open loop (DF1+Cp)*Cv*F2*P
                G = Lp;
            else
                Gp = Lp./(1+Lp);
                if Type == 7
                    % Case 7 => Closed position loop Lp/(1+Lp)
                    G = Gp;
                elseif Type == 8,
                    % Case 8 => Closed position loop with prefilter PF*Lp/(1+Lp)
                    [PreFilt01,Np] = Newsection2coeff(vPosFilt1(1), vPosFilt1(2), vPosFilt1(3), vPosFilt1(4), vPosFilt1(5), 1/Ts);
                     PF = AgitoFiltFreqs(PreFilt01,vPlantFreq*2*pi,1/Ts,Np);

                    G = PF.*Gp;
                elseif Type == 9,
                    % Case 9 => Input disturbance to position P/(1+Lp)
                    G = P./(1+Lp);
                else
                    % Case 10 => Noise to position -Lp/(1+Lp)
                    G = -Lp./(1+Lp);
                    
                end 

            end
        end
    end
end
% vOutGain = 20*log10(abs(G.r(:)));
% vOutPhase = angle(G.r(:))*180/pi;

vOutGain = 20*log10(abs(G));
vOutPhase = unwrap(angle(G))*180/pi;

%
% Output variables
%
vOutGain1000 = 0 * vPlantFreq1000;      % generate zeros vector.
vOutPhase1000 = 0 * vPlantFreq1000;      % generate zeros vector.
vOutGain1000(1:ActualLength) = vOutGain(1:ActualLength);
vOutPhase1000(1:ActualLength) = vOutPhase(1:ActualLength);
%
GM_Pos = 0;
GM_Neg = 0;
PM = 0;
CrossOver = 0;
dB3Freq = 0;
PeakValue = 0;
PeakValueFreq = 0;
DisturbancePeakValue = 0;
DisturbancePeakValueFreq = 0;
NoiseHighestFreqGainValue = 0;
NoiseHighestFreqValue = 0;
Noise3dBFreqValue = 0;

if Type >= 1 && Type<=5,
    L = Lv;
    H = D./(1+Lv);
    P2 = P;
elseif Type >= 6 && Type<=10,
    L = Lp;
    H = Lp./(1+Lp);
    P2 = P./Lp;
else
    return
end
[GM_Pos,GM_Neg,Pm,Wpm] = Complex2Margin(L,vPlantFreq);%%GM_Pos in dB, GM_Neg in dB
debugon=0; %%sophia,will set to zero;
if debugon
  GM_Pos
  GM_Neg
   Pm
   Wpm
   pause
end
Q = L./(1+L);

[PeakVal,Fpeak,F3dB] = PeakAndFreq(Q,vPlantFreq);
[PeakValD,FpeakD,~] = PeakAndFreq(P2.*H,vPlantFreq);
[~,~,F3dBN] = PeakAndFreq(H,vPlantFreq);
debugon=0; %%sophia,will set to zero;
if debugon
    ComplexBodePlot(Q,vPlantFreq)
end
PM = Pm;
if Wpm ==  NaNRepresent
    CrossOver=  NaNRepresent;
else
    CrossOver = Wpm/2/pi;
end
dB3Freq = F3dB;
PeakValue = PeakVal;
PeakValueFreq = Fpeak;

DisturbancePeakValue = PeakValD;
DisturbancePeakValueFreq = FpeakD;

NoiseHighestFreqGainValue = 20*log10(abs(H(end)));
NoiseHighestFreqValue = vPlantFreq(end);
Noise3dBFreqValue = F3dBN;


end


function [PeakVal,Fpeak,F3dB] = PeakAndFreq(Q,f)

g = 20*log10(abs(Q));
gp3 = g+3;
sgp3 = sign(gp3);
ind2cross = find(diff(sgp3)<0);
if ~isempty(ind2cross)
    Index2Db3 = ind2cross(1);
    F3dB = interp1(g(Index2Db3:Index2Db3+1),f(Index2Db3:Index2Db3+1),-3);
else
    F3dB = 0;
end
[PeakVal,PeakIndex]=max(g);
Fpeak = f(PeakIndex);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

