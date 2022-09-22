function [x] = func_FCS(interval,photon)
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
G = @(x,time)   1.0/x(1) ./ (1+time/x(2));
x0=[10.0 0.00010];
[x,resnorm,~,exitflag,output] = lsqcurvefit(G,x0,tau,g);
semilogx(tau,g)
hold on
semilogx(tau,G(x,tau))
hold off
ylim([-0.1 1.1]);
end