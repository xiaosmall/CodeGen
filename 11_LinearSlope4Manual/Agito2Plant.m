function plant = Agito2Plant(ExcelFileName,MatlabFileName)
%
%
% plant = Agito2Plant(ExcelFileName,MatlabFileName)
% Convert Agito format frequency response data file into Auto Loop Shaping
% format Matlab file.
%
% Inputs:
% ExcelFileName - String with Agito frequency response data file name.
% MatlabFileName - String with Matlab file name to save the Auto Loop Shaping 
% plant frequency response data file
% 
%
% Output:
% plant - An FRD Matlab model that contains the frequency response.
%
% Example:
% plant = Agito2Plant('MeasuredFreqs.xls','PlantExample');

data = xlsread(ExcelFileName);
[~,n] = size(data);
m = n/3;

fm = data(:,3:3:m*3);
uf = fm(:);
fv = sort(uf);
i0 = diff(fv)==0;
fv(i0) = [];
plant = frd(zeros(length(fv),1),fv,'unit','hz');
for k=1:m,
    ak = 10.^(data(:,3*(k-1)+1)/20);
    pk = data(:,3*(k-1)+2)*pi/180;
    fk = data(:,3*(k-1)+3);
    i0 = diff(fk)==0;
    fk(i0) = [];
    ak(i0) = [];
    pk(i0) = [];
    
    an = interp1(fk,ak,fv,'linear','extrap' );
    pn = interp1(fk,pk,fv,'linear','extrap' );
    
    z = an.*exp(1i*pn);
    plant(k) = frd(z,fv,'unit','hz');
end
if nargin==2
    save(MatlabFileName,'plant');
end