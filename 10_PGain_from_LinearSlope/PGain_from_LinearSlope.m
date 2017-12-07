function [StartFrequency, EndFrequency, SlopeResult, NumberOfPoints,PlantGain, nErrCode] = PGain_from_LinearSlope(inf, r, n, MinFreqRange, MaxFreqRange, Slope, Tolerance,vPh)
% [StartFrequency, EndFrequency, SlopeResult, NumberOfPoints,PlantGain] = PGain_from_LinearSlope(f, r, n, MinFreqRange, MaxFreqRange, Slope, Tolerance)
%Ver02 - add nErrCode, other bug is fixed at the LinearSlopeRange(Ver06)
%Ver03 -  df chg to 0.05 from 0.01,otherwise size too large;
%         MaxMemorySize is chged from 1000 to 1e7;
%#codegen
%%nErrCode
% -1 : input parameter wrong: MinFreqRange or MaxFreqRange outof range;
% -2 : input parameter (MinFreqRange> MaxFreqRange) conflict;
% -3: not find any period within the tolorance;
% -5; calculation fail;
% -10: calculate failed;
% -11: freq range is out of middle variable size;
%to use phase info to limit the slope searching range, any freq whose phase
%smaller than -240deg, will be not the linear slope range ,will not be the
%search range:
% idxPhLarger180 = find( phdeg <=( -180-60 ));
% EndIndex = idxPhLarger180(1);
%%%
% 8 input:
% inf, r, n, MinFreqRange, MaxFreqRange, Slope, Tolerance,vPh
% 6 output
% % StartFrequency, EndFrequency, SlopeResult, NumberOfPoints,PlantGain, nErrCode
%%%%%%%%~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
StartFrequency = 0;         % to be ready if we return due to some error
EndFrequency = 0;
SlopeResult = 0;
NumberOfPoints = 0;
PlantGain = 0;
nErrCode = 0;
MaxMemorySize = 1e7;
gainLineRange = zeros(MaxMemorySize,1);%%1000 is for the c codegen,since c code will assign other n number of point to other non-zero value, this will affect to the following sum calculation

%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%input value protection;
%-------------------------------------------------------------------------%
f = inf(1:n);
if MinFreqRange> MaxFreqRange
    nErrCode = -2;
end
if nErrCode<0
    return;
end
if MinFreqRange< f(1) || MaxFreqRange > f(end)
    nErrCode = -1;
end
if nErrCode<0
    return;
end
%-------------------------------------------------------------------------%
% Generate fix resolution response vector (0.1Hz)
%
df = 0.05;
fnew = f(1):df:f(n);fnew = fnew(:);
% rdb = interp1(f(1:n),r(1:n),fnew);
%to solve f is non-unique problem;
[f_afterUnique, index_f_afterUnique] = unique(f); %sort and remove duplicates
r_afterUnique = r(index_f_afterUnique);
n = length(index_f_afterUnique);
rdb = interp1(f_afterUnique(1:n),r_afterUnique(1:n),fnew);
%EndIndex = length(fnew);
%improve @20160728
phdeg = interp1(f_afterUnique(1:n),vPh(1:n),fnew);%%improve program,add vPh info
%%%not use the ending freq;but use the vPh=-180deg to set the EndIndex;
idxPhLarger180 = find( phdeg <=( -180-60 ));
if isempty(idxPhLarger180)
    idxPhLarger180 = find( phdeg == min(phdeg));
end
EndIndex = idxPhLarger180(1);
%
%
tol = abs(Tolerance*Slope);
log10f = log10(fnew);
MinIndexRange = fix(MinFreqRange/df);
MaxIndexRange = fix(MaxFreqRange/df);
%
% Must be at least the minimal number of elements.
%
if MinFreqRange > (f(n)-f(1)),
    nErrCode = -5; %%???
    return;
end
%
y = zeros(EndIndex,1);
%
% Claculate the slopes of all the gain responses in MinIndexRange window
% and put the results in y vector.
%
for k=1:EndIndex-MinIndexRange
    H = [log10f(k:k+MinIndexRange-1) ones(MinIndexRange,1)];
    x = H\rdb(k:k+MinIndexRange-1);
    y(k) = x(1);
end
%
% Calculate the longest series of successive
% indexes where the slope meets the required
% accuracy.
% i2mean is index into the mean value of these
% series
%
% kk - contains the indexes of all the slopes that meet the required
% accuracy
%
kk = find(abs(y-Slope)<tol);
kk = kk(diff(kk)==1);
if isempty(kk),
    nErrCode = -3; %%tolerance too tighten,pls increase tolerance
    return;             % None found, error
end
%
% jj - is set of indexes of kk. kk(jj(1)), kk(jj(2), ..., kk(jj(N)
% are the indexes that point to the starting and the ending of sequence of kk.
jj=find(diff(kk)~=1);
%%modify @20160728
if isempty(jj), %means only one stage inside of the tolorance
    StartFrequency = fnew(kk(1));
    EndFrequency = fnew(kk(end));
    SlopeResult = mean(y(kk(1:end)));
    NumberOfPoints = length(kk);
    nErrCode =0; %means success;
    %%if return here,calculate the plantGain;!!!!!
    for ii=1:length(kk(1:end))
        gainLineRange(ii) = abs(10.^(rdb(kk(ii))/20))*fnew(kk(ii))*2*pi*fnew(kk(ii))*2*pi;
    end
    PlantGain  = sum(gainLineRange)/length(kk); %%can not use mean of matlab for the zero value
    PlantGain = PlantGain(1);
    %%%finish cal plant gain
    return;             % No need continue,following to find the best one during many stages;
end
jj = [1;jj;jj+1];
jj = sort(jj);
jj(end) = [];
%
% ii - is an index of jj. kk(jj(ii)) is the index to starting point of the
% longest sequence.
%
% Note that we are looking for the longest sequence from the point of view
% of a user looking at the Bode graph. This means, when the X axis is
% logarithmic.
% As a result, we use the "/" (division) of end and start frequencies of
% each range (and not "-").
%
[~,ii]=max(kk(jj(2:2:end))./kk(jj(1:2:end-1)));%%bug:if the tolorance too small,the ii maybe empty
i2mean = round(sqrt(kk(jj(ii+1))*kk(jj(ii))));   % Geometrical mean
%
% Calculate the permmited range that the min range can be extended to
%
l = zeros(3,1);
l(1) = 2*(i2mean-1)+MinIndexRange;                                            % limitation by the first index of the vectors.
l(2) = 2*(EndIndex - ((i2mean-1)+MinIndexRange))+MinIndexRange;      % limitation by the last index of the vectors.
l(3) = MaxIndexRange;                                                               % limitation as set by the user.
min_l = min(l);
%
% Claculate the slopes of all gain responses around the optimal point
% and put the results in z vector.
%
z = zeros((min_l-MinIndexRange)/2,1);
for k=1:(min_l-MinIndexRange)/2
    IndexRange=i2mean-k:i2mean+MinIndexRange+k-1;
    H = [log10f(IndexRange) ones(length(IndexRange),1)];
    x = H\rdb(IndexRange);
    z(k) = x(1);
end
%
zms = z-Slope;
%
% Takes the last one, meaning the wider window that still sarisfies the slope demand
%
k = find(abs(zms)<tol,1,'last');
flag = (k~=0);
%
% Take the values only if a solution was found
%
if ((~isempty(flag)) && all(flag)),                 % all() is required by codegen, checking isempty is required as all([]) = 1, surprisingly
    k1 = i2mean-k;
    k2 = i2mean+MinIndexRange+k-1;
    StartFrequency = fnew(k1);
    EndFrequency = fnew(k2);
    SlopeResult = z(k);
    NumberOfPoints = k2 - k1 + 1;
else
    %
    % We didn't suceed to extend the initial solution, so we return the initial solution
    %
    StartFrequency = fnew(i2mean);
    EndFrequency = fnew(i2mean+MinIndexRange-1);
    SlopeResult = y(i2mean);
    NumberOfPoints = MinIndexRange;
end
%
% Required for codegen to create scalar output parameters
%
StartFrequency = StartFrequency(1);
EndFrequency = EndFrequency(1);
SlopeResult = SlopeResult(1);
NumberOfPoints = NumberOfPoints(1);

%%%for  plant gain calculation;
idxG = find(  fnew > StartFrequency &  fnew < EndFrequency);
if isempty(idxG)
    nErrCode = -10;%program inside wrong;calculate failed;
    return;
end
if length(idxG) > MaxMemorySize
	nErrCode = -11;%freq range is out of middle variable size;
	return;
end
for ii=1:length(idxG)
    gainLineRange(ii) = abs(10.^(rdb(idxG(ii))/20))*fnew(idxG(ii))*2*pi*fnew(idxG(ii))*2*pi;
end
PlantGain  = sum(gainLineRange)/length(idxG); %%can not use mean of matlab for the zero value
PlantGain = PlantGain(1);

