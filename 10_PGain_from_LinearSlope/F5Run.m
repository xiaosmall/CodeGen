% FileName = 'Linear Stage Far Resonance Set 1.xls'%'MeasuredFreqs.xls'%'MeasuredFreqs_4mat2011.xls';
clear;clc;close all;
% FileName = 'Plant_withLF_Res_realcase.xls'%
plantCase =1;
switch plantCase
    case 1,
        FileName = 'Plant_withLF_Res_realcase.xls'%
        P = Agito2Plant(FileName);
        figure;bodeplot(P);grid on;
        fmaxSlope = 200;%420;
        fminSlope = 10;
        Slope = -40;
        r = P.r;r=r(:);
        rdb = 20*log10(abs(r));
        f = P.f;       
        p = angle(r)*180/pi;
        n = length(r);
        SlopeTolerancePercentage = 0.9;
    case 2
        filenamecase = 1;
        switch filenamecase
            case 'FailedCase'
                addpath('D:\09_servo\09_rotoryMotor');
                %                 FileName = 'a_f1.xls';
            case 4
                FileName = 'B_Axis_f1.xls';
            case 5
                addpath('D:\09_servo\10_AVA40');
                FileName = 'PWM.xls';
            case 6
                addpath('D:\09_servo\02_MGV41-RUN');
                FileName = 'MGV41_NOLoad.xls';
            case 2
                                         FileName = 'Linear Stage Far Resonance Set 2.xls';
                %FileName = 'Linear Stage Far Resonance Set 3.xls';
                addpath('D:\07_Ctr\42_Complex2Margin');
            case 1
                addpath('D:\09_servo\12_AKBdemoSet');
                FileName = 'A_Axis1_design-f1.xls';
                
        end
        P = Agito2Plant(FileName);
        Slope = -40;
        r = P.r;r=r(:);
        rdb = 20*log10(abs(r));
        f = P.f;f=f(:);
        p = angle(r)*180/pi;
        n = length(r);
        d = 1000-n;
        rdb = [rdb;zeros(d,1)];
        f = [f;zeros(d,1)];
        p = [p;zeros(d,1)];
        debugFlag = 1;
        if(debugFlag==1)
            figure(100);semilogx(f,rdb,'r');grid on;ylabel('Mag,dB');hold on;
            figure(101);semilogx(f,p,'r');grid on;ylabel('Ph,deg')
            % %             figure(100);plot(f,rdb,'or');grid on;ylabel('Mag,dB')
            % %             figure(101);plot(f,p,'or');grid on;ylabel('Ph,deg')
            %     plot(fnew);grid on;hold on;
        end
        switch filenamecase
            case 6
                fminSlope =f(1);
                fmaxSlope =750;%Hz
                SlopeTolerancePercentage = 0.3
            case 5
                fminSlope =f(1);
                fmaxSlope =300;%Hz
                SlopeTolerancePercentage = 0.2
            case 4
                fminSlope =10%%%f(1);
                fmaxSlope =200;%Hz
                SlopeTolerancePercentage = 0.4
            case 3
                fminSlope =f(1);
                fmaxSlope =200;%Hz
                SlopeTolerancePercentage = 0.9
            case 1
                fminSlope =f(1);
                fmaxSlope =295;%Hz
                SlopeTolerancePercentage = 0.1
                
            case 2
                fminSlope =f(1);
                fmaxSlope =550;%Hz
                SlopeTolerancePercentage = 0.1
                
            otherwise
                fminSlope =f(1);
                fmaxSlope =f(end);
                SlopeTolerancePercentage =0.1;
                
        end
        Slope = -40;
    case 3
        
end
[SlopeStartFrequency, SlopeEndFrequency, SlopeResult, NumberOfPoints,PlantGain, nErrCode] = PGain_from_LinearSlope(f,rdb,n,fminSlope,fmaxSlope,Slope,SlopeTolerancePercentage,p)
% % % [SlopeStartFrequency, SlopeEndFrequency, SlopeResult, NumberOfPoints,nErrCode] = LinearSlopeRange(f,rdb,n,fminSlope,fmaxSlope,Slope,SlopeTolerancePercentage,p)
% %             figure(100);semilogx([SlopeStartFrequency:SlopeEndFrequency],rdb,'r');grid on;ylabel('Mag,dB');hold on;
% Verify_PGain_from_LinearSlope