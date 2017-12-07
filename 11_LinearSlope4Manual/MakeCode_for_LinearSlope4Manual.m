Z = 1.0;                        % just a dummy float double value;
vZ = zeros(1000,1);     % just a dummy float double vector of size 1000, for inputfrequency/gain/phase;
ToRun = 'NotRun'%'ToRun'
if strcmp(ToRun,'ToRun')
    f =  vZ;
    r = vZ;
    n= length(f);
    MinFreqRange = f(1);
    MaxFreqRange = f(end);
    [SlopeResult, PlantGain,ErrCode] = LinearSlope4Manual(f, r, n, MinFreqRange, MaxFreqRange);
end
codegen -args {vZ, vZ, Z, Z, Z} -report -config:lib     LinearSlope4Manual.m

