function [SlopeResult, PlantGain,ErrCode] = LinearSlope4Manual(f, r, n, MinFreqRange, MaxFreqRange)
%   Input parameters:
%
%   f - Frequencies vector. Assumed to have fixed length of 1000 (needs fixed length for codegen) 
%   r - Gains vector, in [dB]. As above.
%   n - Actual length of the f and r vectors.
%   MinFreqRange - Start frequency range.
%   MaxFreqRange - End frequency range.
%   Output parameters:
%Output---------
%   SlopeResult - the calculated slope for this frequency range.
%   PlantGain - Calculated gain of plant.
%   ErrCode --- if negative value, Daniel will display the ErrCode;
%           --- [-3]--- input value out of range;
%   All output parameters are 0 if the function fails to find a solution.
%   plantGain calculation is based on each freq point(w) and G/(jw)^2,then make average of the whole interpolated freq points 
%#codegen
%init output
SlopeResult = 0;
PlantGain = 0;
ErrCode  = 0;
% var = f;
% if isempty(var) || isnan(var) 
%     ErrCode  = -3;
% end
ErrCode = inputVarCheck(f,-3);
ErrCode = inputVarCheck(r,-3);
ErrCode = inputVarCheck(n,-3);
ErrCode = inputVarCheck(MinFreqRange,-3);
ErrCode = inputVarCheck(MaxFreqRange,-3);

if ErrCode <0
    return;
end

df = 0.1;
if MinFreqRange <= 0 | MaxFreqRange <=0
    ErrCode = -4;
end
if ErrCode <0
    return;
end

fnew = f(1):df:f(n); %% 1*XX 
%%
% nLength = length(fnew);
% % fnew = fnew(:);%% do we need to XX*1 ? temp, will uncomment later;
% % % // Den = reshape(Den,1,numel(Den));%%1*XX matrix transopse eg
rdb = interp1(f(1:n),r(1:n),fnew);
idxG = find(  fnew > MinFreqRange &  fnew < MaxFreqRange);
if isempty(idxG)
	ErrCode = -1;
    return;
end
ActiveLength = length(idxG);
if ActiveLength > length(fnew)
	ErrCode = -2;
    return;
end
xx= log10(fnew(idxG(1):idxG(ActiveLength)));

[pCoefficient,~] = polyfit(xx,rdb(idxG(1):idxG(end)),1);
SlopeResult = pCoefficient(1);
%PlantGain22 = (2*pi)^2*10^(pCoefficient(2)/20)

%%%for  plant gain calculation;
gainLineRange = zeros(1000,1);
rnondb = 10.^(rdb/20); %%non-dB

for ii=1:ActiveLength
      gainLineRange(ii) = abs(rnondb(idxG(ii)))*fnew(idxG(ii))*2*pi*fnew(idxG(ii))*2*pi;
end
vPlantGain  = sum(gainLineRange)/length(idxG); %%can not use mean of matlab for the zero value 
PlantGain = vPlantGain(1);