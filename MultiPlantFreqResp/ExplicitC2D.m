function [a,b] = ExplicitC2D(n,d,SampleFreq,N,w)
%
% [a,b] = ExplicitC2D(n,d,SampleFreq,N)
%
% Inputs:
%
% n - N X 1 array that conatins the continue transfer function numerator coeffcients.
% d = N X 1 array that conatinscontinue transfer function denominator coefficients.
%
%                               n(1)*s^2+n(2)*s+n(1)/d(1)*s^2+d(2)*s+d(1)
%
% SampleFreq - Sample Frequency in Hz.
% N - filter order+1 can have the value 3 or 2;
%
% Outputs:
%
% a - N X 1 array that conatins the discrete transfer function numerator coeffcients.
% b - N X 1 array that conatins the discrete transfer function numerator coeffcients.
%
% TODO: consider using prewrap for relatively high frequencies.
%
k = 2 * SampleFreq;              % Default bilinear method
%
if nargin==5,
    k = w/tan(w / (2 * SampleFreq));% Blinear with pre warping at frequency w
end
%
k2 = k^2;
%
% Applying bilinear transformation on the continue coefficient to have the
% discrete filter coefficients
%
a = zeros(N,1);
b = zeros(N,1);
if N == 3,
    a(1) = k2*n(1) + k*n(2) + n(3);
    a(2) = -2*k2*n(1) + 2*n(3);
    a(3) = k2*n(1) - k*n(2) + n(3);
    %
    b(1) = k2*d(1) + k*d(2) + d(3);
    b(2) = -2*k2*d(1) + 2*d(3);
    b(3) = k2*d(1) - k*d(2) + d(3);
else
    a(1) = k*n(2) + n(3);
    a(2) = -k*n(2) + n(3);
    %
    b(1) =  k*d(2) + d(3);
    b(2) = -k*d(2) + d(3);
end
%
% Normalize the coefficients so b(1) = 1
%
for i=N:-1:1,
    a(i) = a(i)/b(1);
    b(i) = b(i)/b(1);
end
debugon = 0;
if debugon
    figure(320); hold off;
    atransv = reshape(a,1,3);
    btransv = reshape(b,1,3);
    bodeplot(tf(atransv,btransv,1/SampleFreq));grid on;
end
