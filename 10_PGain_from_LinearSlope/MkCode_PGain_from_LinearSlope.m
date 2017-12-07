Z = 1.0;                        % just a dummy float double value;
vZ = zeros(1000,1);     % just a dummy float double vector of size 1000, for inputfrequency/gain/phase;
codegen -args {vZ, vZ, Z, Z, Z, Z, Z, vZ} -report -config:lib     PGain_from_LinearSlope.m

