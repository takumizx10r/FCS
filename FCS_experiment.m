clear
close all
dt=4.0*10^(-6);
MaxIter=1*10^5;
Resolution=0.170;
[imfile, path, indx]=uigetfile(strcat(pwd,'\*.tif'));
imdata=imread(strcat(path,'/',imfile));

if MaxIter==inf
    TotalIteration=size(imdata,1)*size(imdata,2)-1;
elseif MaxIter>size(imdata,1)*size(imdata,2)
    display('MaxIter exceeds image size.');
    return;
else
    TotalIteration=MaxIter;
end

for i=0:TotalIteration
    tdata(i+1)=i*dt;
    N_photon(i+1)=imdata(fix(i/size(imdata,2))+1,rem(i,size(imdata,2))+1);
end
figure
plot(tdata,N_photon)

% % % Fitting Autocorrelation function
display('Fitting Autocorrelation function')
[FitPara]=func_FCS(dt,N_photon,Resolution,0.0);
N_detect=FitPara(1);
D_eff=(Resolution/2.0)*(Resolution/2.0) / (4.0*FitPara(2));
% % % 
cd (path);