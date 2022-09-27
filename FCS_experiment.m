clear
close all
dt=2.0*10^(-6); %%sec
MaxIter=1*10^5;
Resolution=0.170; %%um
[imfile, path, indx]=uigetfile(strcat(pwd,'\*.tif'));
imdata=imread(strcat(path,'/',imfile));

if MaxIter==inf
    TotalIteration=size(imdata,1)*size(imdata,2)-1;
elseif MaxIter>size(imdata,1)*size(imdata,2)
    disp('MaxIter exceeds image size.');
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

% % % Fitting Autocorrelation function to determine DZ
[FitPara]=func_FCS_fit_resolutions(dt,N_photon,Resolution);
N_detect=FitPara(1);
Dz=FitPara(2);

% % % Fitting Autocorrelation function
% [FitPara]=func_FCS(dt,N_photon,Resolution,0.60);
% N_detect=FitPara(1);
% D_eff=FitPara(2);
% % % 
cd (path);