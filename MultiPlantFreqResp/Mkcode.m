clear;clc;close all;
Z = 1.0;                        % just a dummy float double value;
vZ = zeros(1000,1);     % just a dummy float double vector of size 1000, for inputfrequency/gain/phase;
f = vZ;
pg= vZ;pp= vZ;n = length(f);
kpp = 100;
kvp = 1000;
kvi = 20;
PosFilt = [0 0 0 0 0];
VelFilt1 = [2,3000,0.5,0,0];
VelFilt2 = [8,350,20,10,0];
VelFilt3 = [2,30,0.5,0,0];
VelFilt4 = [2,3000,0.5,0,0];
LoopType=1;requiredFreq=80;GreaterLessFlag=1;
Ts = 1/16384;
codegen -args {f, pg,pp,n,f,pg,pp,n,f,pg,pp,n,f,pg,pp,n,f,pg,pp,n,f,pg,pp,n,f,pg,pp,n,f,pg,pp,n,f,pg,pp,n,Z,PosFilt,0,kpp,0,kvp,kvi, 1,VelFilt1, 1,VelFilt2, 1, VelFilt3, 1,VelFilt4, 1,Ts,LoopType,requiredFreq,GreaterLessFlag} -report -config:lib MultiPlantFreqResp.m
%%
