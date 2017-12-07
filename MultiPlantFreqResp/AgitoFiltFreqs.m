function Freqs = AgitoFiltFreqs(AgitoFilt,w,Fs,N)

Hd=[(exp(w*1i/Fs)).^2 exp(w*1i/Fs) ones(size(w))];
nd = zeros(1,3);
dd = zeros(1,3);
if N == 3
    nd = AgitoFilt(1:3)/65536;
    dd = [1 AgitoFilt(4:5)/65536];
elseif N==2,
    nd = [0 AgitoFilt(1:2)/65536];
    dd = [0 1 AgitoFilt(4)/65536];
end

nd = nd(:);
dd = dd(:);

Freqs = (Hd*nd)./(Hd*dd);