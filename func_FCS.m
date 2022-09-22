function [x] = func_FCS(interval,photon)

prompt = {'Initial N:','Initial tau_D:'};
dlgtitle = 'Initial values';
dims = [1 35];
definput = {'10','3'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
Ni=str2num(answer{1});
Di=str2num(answer{2});

ave_photon=mean(photon);
photon=photon-ave_photon;
N_max=size(photon,2);
for j=0:N_max-1
    
    for i=1:N_max-j
        keep(i)=photon(i)*photon(i+j);
    end

    g(j+1)=mean(keep);
    tau(j+1)=interval*j;
    clear keep
end
% g=g./( ave_photon*ave_photon );
g=g./( mean(photon.*photon) );
% % % Fitting 
display('Fit Autocorrelation Function');
n=1;
while n>0
    G = @(x,time)   1.0/x(1) ./ (1+time/x(2));
    x0=[Ni Di];
    [x,resnorm,~,exitflag,output] = lsqcurvefit(G,x0,tau,g);
    semilogx(tau,g)
    hold on
    semilogx(tau,G(x,tau))
    hold off
    ylim([-0.1 1.1]);

    answer = questdlg('Do you like to fit again with different I.V.?', ...
    	'Yes or No', ...
    	'Yes','No','Cancel');
    % Handle response
    switch answer
        case 'Yes'
            prompt = {'Initial N:','Initial tau_D:'};
            dlgtitle = 'Initial values';
            dims = [1 35];
            definput = {'10','3'};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            Ni=str2num(answer{1});
            Di=str2num(answer{2});
        case 'No'
            n=0;
    end

end

end