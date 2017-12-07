function [AgitoFilt,N] = Newsection2coeff(NumType, par1, par2, par3, par4, SampleFreq)
%
% [AgitoFilt, Filter] = Newsection2coeff(NumType, par1, par2, par3, par4, SampleFreq)
%
% The function gets as input the filters type and the filters parameters 
%  in term of user language and return discrete filter coefficients vector 
%  that are used in the following filter implementation:
%
% Y = (X * AgitoFilt[1] + XPrev * AgitoFilt[2] + XPrevPrev * AgitoFilt[3] - Y * AgitoFilt[4] - YPrev * AgitoFilt[5]) / 65536.0;
% AgitoFilt[1-5] are integer variables. They are the filter coefficients, multiplied by a constant of 65536.
% X  The filetr input.
% Y  The filetr output.
%
% It also returns the filter as transfer function, to be used by CalculateAnyFrequencyResponse()
%
% Inputs:
%
% NumType - value that defines the filter type.
% par# - are the filter parameters in user language.
% SamplFreq = control sampling frequency in Hz
% 
% Outputs:
% AgitoFilt - Vector of the filetr coefficinents
% Filter - transfer function of the filter
%   
%
% The filters types and their corresponding parameters are detailed below:
%
%       0 - NONE: no filter. Filter is set to be 1.
%
%       1 - LPF1 - 1st order Low Pass Filter the par1 is the pole frequency in
%              Hz, par2,3, and 4 aren't applicable.
%
%       2 - LPF2 - 2nd order Low Pass Filter, par1 is the pole frequency in
%              Hz, par2 is the pole damping ratio. par3, and 4 aren't applicable.
%
%       3 - LPF3 - 2nd order Low Pass Filter with zero,
%              par1 is the pole frequency in
%              Hz, par2 is the pole damping ratio and par3 is the zero frequency
%              in Hz. par4 isn't applicable.
%
%       4 - LDLG2 -2nd order Lead/Lag filter
%              par1 is the lead first zero , par2
%              is the lead second zero, par3 is the leg first pole and
%              par4 is the leg second pole all the frequency values are in Hz
%
%       5 - LDLG1 - 1st order Lead/Lag filter
%              par1 is the lead zero and par2 is the leg pole and
%              all the frequency values are in Hz
%       
%       6- LDLG1FP - 1st order Lead/Lag filter defined by the frequency at the 
%               phase peak/min and the phase level.
%               par1 is the frequency at the phase's peak/min in Hz, 
%               par2 in the required phase peak in deg.
%
%       7 - LDLG2FP - 2nd order Lead/Lag filter defined by the frequencies at
%               the phases peak/min and the phases level (Approximately)
%               par1 is the 1st frequency at the phase peak/min in Hz, 
%               par2 in the required phase peak in deg at the 1st frequency.
%               par3 is the 2nd frequency at the phase peak/min in Hz, 
%               par4 in the required phase in deg at the 2nd frequency.
%
%
%       8 - NOTCH - Notch filter
%               par1 is the notch frequency in Hz, par2 is
%               the notch depth in dB and par 3 is the Notch width in Hz par4 isn't applicable.
%
%       9 - CLDLG - Complex Lead/Leg filter.
%               par1 - frequency of the complex zero in Hz, par2 damping
%               ratio of the complex zero , par3 is the frequency of the
%               pole in Hz and par4 is the damping ratio of pole.
%
%       10 - BQUD - General bequad filter
%              par1 contains the continue time domain denominator coefficients vector and
%              par2 is the continue time domain nominator coefficients vector.
%
% Output:
%
% AgitoFilt - vector of discrete filter coefficients that defined by the
% 				  above equation.
%
% Examples:
%
%       AgitoFilt = section2coeff('LPF1',30,[],[],[],16384);    % 1st order lowpass with cut-off at 30 Hz       
%       AgitoFilt = section2coeff('LPF2',300,0.5,[],[],16384);  % 2nd order lowpass with cut-off at 300 Hz          
%       AgitoFilt = section2coeff('LPF3',300,0.5,300,[],16384);  % 2nd order with zero lowpass with cut-off at 300 Hz       
%       AgitoFilt = section2coeff('LDLG2',80,120,90,100,16384);         
%       AgitoFilt = section2coeff('LDLG1',80,120,[],[],16384);         
%       AgitoFilt = section2coeff('NOTCH',100,20,5,[],16384);          
%       AgitoFilt = section2coeff('CLDLG',100,0.1,150,0.5,16384);      
%       AgitoFilt = section2coeff('LDLG1FP',100,20,[],[],16384);      
%
% See also: Filt2FrequencyResponse

%       Author(s): Y. Zimmerman 28-12-2012,
%       Copyright 2012 The Agito, Inc.
%       $Revision: 1.0.0 $ $Date: 2012/12/25 16:52:48 $

%
% Define the coefficients in continue domain (Laplace transfer function)
%
N = 2;
nd = zeros(1,3);
dd = zeros(1,3);
nc = zeros(1,3);
dc = zeros(1,3);

switch NumType
    case 1 %'LPF1'
        p = par1*2*pi;
        nc = [0 0 p];
        dc = [0 1 p];
        N = 2;
        [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N);
%
    case 2 %'LPF2'
        p = par1*2*pi;
        xi = par2;
        nc = [0 0 p^2];
        dc = [1 2*xi*p p^2];
        N = 3;
        [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N);
%
    case 3 %'LPF3'
        p = par1*2*pi;
        xi = par2;
        z = par3*2*pi;
        nc = [0 p^2/z p^2];
        dc = [1 2*xi*p p^2];
        N = 3;
        [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N);
%
    case 4 %'LDLG2'
        z1 = par1*2*pi;
        z2 = par2*2*pi;
        p1 = par3*2*pi;
        p2 = par4*2*pi;
        nc = [1/(z1*z2) 1/z1+1/z2 1];
        dc = [1/(p1*p2) 1/p1+1/p2 1];
        N = 3;
        [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N);
%
    case 5 %'LDLG1'
        z1 = par1*2*pi;
        p1 = par2*2*pi;
        nc = [0 1/(z1) 1];
        dc = [0 1/(p1) 1];
        N = 2;
        [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N);
%
    case 6 %'LDLG1FP'
        wp = par1*2*pi;
        pp = par2*pi/180;
        w1 = wp*(sqrt(tan(pp)^2+1)-tan(pp));
        w2 = wp*(sqrt(tan(pp)^2+1)+tan(pp));
        nc = [0 1/(w1) 1];
        dc = [0 1/(w2) 1];
        N = 2;
        [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N,wp);
%
    case 7 %'LDLG2FP'
        wp1 = par1*2*pi;
        pp1 = par2*pi/180;
        wp2 = par3*2*pi;
        pp2 = par4*pi/180;
        w1 = wp1*(sqrt(tan(pp1)^2+1)-tan(pp1));
        w2 = wp1*(sqrt(tan(pp1)^2+1)+tan(pp1));
        w3 = wp2*(sqrt(tan(pp2)^2+1)-tan(pp2));
        w4 = wp2*(sqrt(tan(pp2)^2+1)+tan(pp2));
        nc = [1/(w1*w3) 1/w1+1/w3 1];
        dc = [1/(w2*w4) 1/w2+1/w4 1];
        N = 3;
        [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N);
%
    case 8 %'NOTCH'
        w = par1*2*pi;
        d = 10^(par2/20);
        dw = par3*2*pi;
        xi = dw/(2*w);
        nc = [1 2*xi*w w^2];
        dc = [1 2*xi*d*w w^2];
        N = 3;
        [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N,w);
%
    case 9 %'CLDLG'
        wz = par1*2*pi;
        xiz = par2;
        wp = par3*2*pi;
        xip = par4;
        nc = [1/wz^2 2*xiz/wz 1];
        dc = [1/wp^2 2*xip/wp 1];
        N = 3;
        [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N);
%
    case 10 %'BQUD'
        nc = par1;
        dc = par2;
        N = 3;
       [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N);
       
    case 0 %'NONE'
        nc = [1 0 0];
        dc = [1 0 0];
        N = 3;
       [nd,dd] = ExplicitC2D(nc,dc,SampleFreq,N);

end
%
AgitoFilt = [0 0 0 0 0];
AgitoFilt(1:N) = nd;
AgitoFilt(4:2+N) = dd(2:N);
%
AgitoFilt = round(AgitoFilt * 65536);

