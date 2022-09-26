function [x] = func_FCS(interval,photon)
N_max=size(photon,2);

prompt = {'Initial N:','Initial tau_D:', 'FitRange(frame):'};
dlgtitle = 'Initial values';
dims = [1 35];
definput = {'10','0.00001',num2str(N_max)};
answer = inputdlg(prompt,dlgtitle,dims,definput);
Ni=str2num(answer{1});
Di=str2num(answer{2});
FitRange=str2num(answer{3});

ave_photon=mean(photon);
photon=photon-ave_photon;
tau=zeros(1,N_max);
tic
for j=0:N_max-1
    keep(j+1:N_max)=photon(1:N_max-j).*photon(j+1:N_max);
    g(j+1)=mean(keep);
    tau(j+1)=interval*j;
    clear keep
end
toc
g=1.0 + g./( ave_photon*ave_photon );
% g=g./( mean(photon.*photon) );
% % % Fitting 
display('Fit Autocorrelation Function');
n=1;
while n>0
    G = @(x,time)  1.0 + 1.0/x(1) ./ (1+time/x(2));
    x0=[Ni Di];
    [x,resnorm,~,exitflag,output] = lsqcurvefit(G,x0,tau(1:FitRange),g(1:FitRange));
    semilogx(tau,g)
    hold on
    semilogx(tau(1:FitRange),G(x,tau(1:FitRange)))
    hold off
%     ylim([-0.1 1.1]);

    answer = questdlg('Do you like to fit again with different I.V.?', ...
    	'Yes or No', ...
    	'Yes','No','Cancel');
    % Handle response
    switch answer
        case 'Yes'
            prompt = {'Initial N:','Initial tau_D:', 'FitRange(frame):'};
            dlgtitle = 'Initial values';
            dims = [1 35];
            definput = {num2str(Ni),num2str(Di),num2str(FitRange)};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            Ni=str2num(answer{1});
            Di=str2num(answer{2});
            FitRange=str2num(answer{3});
        case 'No'
            n=0;
    end

end

end