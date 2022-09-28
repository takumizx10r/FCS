clear
close all
dt=4.0*10^(-6); %%sec
MaxIter=inf;
Resolution=0.170; %%um
[imfile, path, indx]=uigetfile(strcat(pwd,'\*.tif'));
imdata=imread(strcat(path,'/',imfile));
[~, name, ext]=fileparts(imfile);

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
ax=gca;
xlabel('Time \itt \rm(s)','FontSize',20)
ylabel('Photon','FontSize',20)
ax.FontSize=18;
axtoolbar('Visible','off');
exportgraphics(gcf, ...
    strcat(path,'/',name,'-photon.png'), ...
    'Resolution',600)
% % % Fitting Autocorrelation function to determine DZ
[FitPara,G]=func_FCS_fit_resolutions(dt,N_photon,Resolution);
N_detect=FitPara(1);
Dz=FitPara(2);
exportgraphics(gcf, ...
    strcat(path,'/',name,'-autocorrelation.png'), ...
    'Resolution',600)
% % % Fitting Autocorrelation function
% [FitPara,G]=func_FCS(dt,N_photon,Resolution,Dz);
% N_detect=FitPara(1);
% D_eff=FitPara(2);
% exportgraphics(gcf, ...
%     strcat(path,'/',name,'-autocorrelation.png'), ...
%     'Resolution',600)
% % % 
cd (path);
save(name)