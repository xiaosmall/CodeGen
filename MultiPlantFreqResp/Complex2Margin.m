%%this one is sub function of ALS; 
%%histoy: Jul-08-2016   add one more ouput(GmNeg) in this function,chg unit to
%%       Sep-21-2016   Output Inf if GM is infinity; this version different from ALS folder one ;
%%        Nov-30-2016    set GmPos=NaN if not find;   
%%GmPos in dB , %GmNeg in dB for codegen ;
%%[input]
%%          L-- complex number
%%          f-- Frequency in Hz
%[output:] 
%          GmPos -- Gain margin positive in dB
%          GmNeg -- Gain margin negative in dB
%          Pm -- phase margin  in deg
%          Wpm -- cross over freq in rad/s;
% 3 type of result for Gm,Pm:(only two types of results,no infinity)
%(1) not found(Nan):
%                   not cross 0 dB, Pm is Nan; 
%                   not cross -180 deg, Gm(pos& neg) is Nan;
% 
%(2) infinity;-- if GmPositive is normal number,GmNegative is
%(3) some number value;

function [GmPos,GmNeg,Pm,Wpm] = Complex2Margin(L,f)
InfValueHack = 123454321;
GmPos = InfValueHack;
GmNeg = InfValueHack;
Pm = InfValueHack;
Wpm = InfValueHack;
%%
NaNRepresent =  999999;
aL = abs(L);                    % aL is the open loop Gain vectoe
pL = unwrap(angle(L))*180/pi;   % pL is the open loop phase in deg.
while pL(1) > 90,
    pL = pL - 360;
end
% Find indexes of all cases that to open loop phase cross the -180 deg
pLm180 = pL+180;
spLm180 = sign(pLm180);
ind2m180cross = find(diff(spLm180)~=0);

% If there is no -180 degree cross the gim amrign is infinity!
if isempty(ind2m180cross),
    GmPos = NaNRepresent;   %%chg by sophia for c
    GmNeg = NaNRepresent;
%     GmPos = Inf;   %%when this case need to remove GMPos and GmNeg in ctrl-cost calculation;
%     GmNeg = Inf;
    % Estimate the gain at -180 degree by gain vector interpolation
else
    Ga = zeros(length(ind2m180cross),1);
    for i=1:length(ind2m180cross),
        Ga(i) = aL(ind2m180cross(i)) + (-180-pL(ind2m180cross(i))) * ...
            diff(aL(ind2m180cross(i):ind2m180cross(i)+1))...
            /diff(pL(ind2m180cross(i):ind2m180cross(i)+1));
    end
    % Search for gains that are the closest to 0 dB
    Ind2MinGm = [];
    Gdb = 20*log10(Ga);
    ind2pos = find(Gdb>=0);
    if ~isempty(ind2pos),
        [~,TempIndex] = min(Gdb(ind2pos));
        Ind2MinGm = ind2pos(TempIndex);
        GmNeg =  1/Ga(Ind2MinGm);%%chg by sophia to fix c code
        GmNeg = 20*log10(GmNeg);
        GmNeg = GmNeg(1);
        %%no need else,  GmNeg is infinit(123454321);the is conditional stable ,
    end
    % For the conditional stable case look for negative gain margin as well
    ind2neg = find(Gdb<0);
    if ~isempty(ind2neg),
        [~,TempIndex] = max(Gdb(ind2neg));
        %Ind2MinGm = [Ind2MinGm ind2neg(TempIndex)];
        GmPos = 1/Ga(ind2neg(TempIndex));%%chg by sophia to fix c code
        GmPos = 20*log10(GmPos);
        GmPos =GmPos(1);
     else %2016-11-30; need to release;
         GmPos = NaNRepresent; %2016-11-30;
    end
    % The gain margin is one over the gain!
   % Gm = 1./Ga(Ind2MinGm);
end

% Find indexes of all cases that to open loop gain cross the 0 dB and the 
% gain curve has negative slope.
saLm1 = sign(aL-1);
ind2cross = find(diff(saLm1)<0);
% If there is no 0 dB crossing the phase margin can't be measured!
if isempty(ind2cross),
    Pm = NaNRepresent;%2016-0927
    Wpm =NaNRepresent;
    return;
    % In case there is only one crossing: calculate the phase at the 
    % crossing point by an interpoltation
elseif length(ind2cross) == 1,    
    Ph = pL(ind2cross(end)) + (1-aL(ind2cross(end))) * ...
        diff(pL(ind2cross(end):ind2cross(end)+1))...
        /diff(aL(ind2cross(end):ind2cross(end)+1));
    Wpm = interp1(aL(ind2cross(end):ind2cross(end)+1),f(ind2cross(end):ind2cross(end)+1),1)*2*pi;

    % In case there are multiply 0 dB crossing, we need to check the range 
    % of the phases at the crossing points in order to determine the stability.
else
    % Find the indexes of all the cases that the phase is higher than -180 degree
    % at the 0 dB crossing points.
    ind2pos = find(pL(ind2cross)>=-180);
    % In case that there is no crossing at phase higher than -180 degree, 
    % the phase margin is calculated by the 1st 0 dB crossing point
    if isempty(ind2pos),
        Ph = pL(ind2cross(1)) + (1-aL(ind2cross(1))) * ...
            diff(pL(ind2cross(1):ind2cross(1)+1))...
            /diff(aL(ind2cross(1):ind2cross(1)+1));
        Wpm = interp1(aL(ind2cross(1):ind2cross(1)+1),f(ind2cross(1):ind2cross(1)+1),1)*2*pi;
        % If there is 0 dB crossing with phase higher than -180 degree.
        % We need to look also if there are 0 dB crossing with phase lower than
        % -180 degree
    else
        ind2neg = find(pL(ind2cross)<-180);
        % No crossing at phase lower than -180 degree, the phase margin is
        % calculated by the last 0 dB crossing point
        if isempty(ind2neg)
            Ph = pL(ind2cross(end)) + (1-aL(ind2cross(end))) * ...
                diff(pL(ind2cross(end):ind2cross(end)+1))...
                /diff(aL(ind2cross(end):ind2cross(end)+1));
            Wpm = interp1(aL(ind2cross(end):ind2cross(end)+1),f(ind2cross(end):ind2cross(end)+1),1)*2*pi;
            % There are 0 dB crossing at phase higher & lower than -180 degree,
            % the stability is determined by the phase of the inetrmidate
            % crossing point.
            % At this point the gain slop is positive at the 0 dB crossing
        else
            i = ones(3,1);
            Phv = ones(3,1);
            % Closet point to 180 degree from the higher range 
            [~,ind2min] = min(pL(ind2cross(ind2pos)));
             % Closet point to 180 degree from the lower range
            [~,ind2max] = max(pL(ind2cross(ind2neg)));
            % Indexes calculation
            i(1) = ind2cross(ind2pos(ind2min));
            i(2) = ind2cross(ind2neg(ind2max));
            % Intermidate crossing point with positive slope
            ind2cross_p = find(diff(saLm1(i(1):i(2) ))>0);
            i(3) = i(1)+ind2cross_p-1;
            % In case the phase is lower than -180 degree at the positive
            % slope ccrossing -> No stable
            if pL(i(3)) > -180,
                Stability = -1;
            % In case the phase is higher than -180 degree at the positive
            % slope crossing -> Stable
            else
                Stability = 1;
            end
            % Calculate the phase at the three critical points
            Phv(1) = pL(i(1))+(1-aL(i(1)))*diff(pL(i(1):i(1)+1))/diff(aL(i(1):i(1)+1));
            
            Phv(2) = pL(i(2))+(1-aL(i(2)))*diff(pL(i(2):i(2)+1))/diff(aL(i(2):i(2)+1));
            
            Phv(3) = pL(i(3))+(1-aL(i(3)))*diff(pL(i(3):i(3)+1))/diff(aL(i(3):i(3)+1));
            % The phase margin is determine accorringo the one is closet to 
            % the instability point (-180 degree)
            [TempPm,ii]=min(abs(Phv+180));
            % Technical correction to enable correct calculation of the 
            % Phase margin at line 133
            Ph = TempPm*Stability-180;
            Wpm = interp1(aL(i(ii):i(ii)+1),f(i(ii):i(ii)+1),1)*2*pi;

        end
    end
end
% Fix the phase value, so it will be between -360 and 0

while Ph>0
    Ph = Ph-360;
end
while Ph<-360
    Ph = Ph+360;
end
% The phase margin is the difference between the phase at the 0 dB
% crossing and -180 degree.
Pm = Ph+180;

end
