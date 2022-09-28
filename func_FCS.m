function [FitPara,g] = func_FCS(interval,photon,R,Dz)
N_max=size(photon,2);
w=0.5*R;
prompt = {'Initial N:','Initial D:', 'FitRange(frame):', 'Dimension:','Move sum:'};
dlgtitle = 'Initial values';
dims = [1 35];
definput = {'10','100',num2str(N_max), '3','3'};
answer = inputdlg(prompt,dlgtitle,dims,definput);
Ni=str2num(answer{1});
Di=str2double(answer{2});
FitRange=str2num(answer{3});
dimension=str2num(answer{4});
movesum=str2double(answer{5});
wz=Dz;
if movesum~=0
    photon=movsum(photon,movesum);
end

ave_photon=mean(photon);
photon=photon-ave_photon;
tau=zeros(1,N_max);
g=zeros(1,N_max);
tic
for j=0:N_max-1
    keep=zeros(1,N_max);
    keep(j+1:N_max)=photon(1:N_max-j).*photon(j+1:N_max);
    tau(j+1)=interval*j;
    %     g(j+1)=mean(keep);
    g(j+1)=mean(keep(j+1:N_max));
    %     clear keep
end
toc
g = 1.0 + g./( ave_photon*ave_photon );

% g=g./( mean(photon.*photon) );
% % % Fitting
disp('Fit Autocorrelation Function...');
n=1;
while n>0
    if dimension==3
        G = @(x,time)  1.0 + 1.0/x(1) ./ (1+4.0*x(2).*time./w^2).*sqrt(1.0./(1.0+4.0*x(2).*time./Dz^2));
    elseif dimension==2
        G = @(x,time)  1.0 + 1.0/x(1) ./ (1+4.0*x(2).*time./w^2);
    end

    x0=[Ni Di];
%     [x,resnorm,~,exitflag,output] = lsqcurvefit(G,x0,tau(1:FitRange),g(1:FitRange));
%     FitPara=x;
    A=[];
    b=[];
    Aeq=[];
    beq=[];
    nonlcon=[];
    lb=[0.0   ,   0.0];
    ub=[1000,   1000.0];
    options = optimset('fmincon');
    options.Algorithm=('interior-point');
    FitPara = fmincon(@LeastSquare,x0,A,b,Aeq,beq,lb,ub,nonlcon,options);

% % % % % % % % % % % %     Make log plot
    disp(FitPara);
    semilogx(tau(1:FitRange),g(1:FitRange),'k*');
    hold on
    semilogx(tau(1:FitRange),G(FitPara,tau(1:FitRange)),'b-');
    hold off
    ax=gca;
    axtoolbar('Visible','off');
    xlabel('Lag \it\tau \rm(s)','FontSize',20)
    ylabel('Autocorrelation \it','FontSize',20)
    ax.FontSize=18;
    legend('Exp.','\itG(\tau)','fontsize',18)
    
% % % % % % % % % % % % % % % % % 
    answer = questdlg('Do you like to fit again with different I.V.?', ...
        'Yes or No', ...
        'Yes','No','Cancel');
    % Handle response
    switch answer
        case 'Yes'
            prompt = {'Initial N:','Initial D:', 'FitRange(frame):', 'Dimension:','Move sum:'};
            dlgtitle = 'Initial values';
            dims = [1 35];
            definput = {num2str(Ni),num2str(Di),num2str(FitRange),num2str(dimension),num2str(movesum)};
            answer = inputdlg(prompt,dlgtitle,dims,definput);
            Ni=str2num(answer{1});
            Di=str2double(answer{2});
            FitRange=str2num(answer{3});
            dimension=str2num(answer{4});
            movesum=str2double(answer{5});
        case 'No'
            n=0;
    end

end
    function Q=LeastSquare(x)
        function G_=func_autocorrelation_function(x,t)
            if dimension==3
                G_=1.0 + 1.0/x(1) ./ (1+4.0*x(2)*t/w^2)./sqrt((1.0+4.0*x(2)*t/wz/wz));
            elseif dimension==2
                G_=1.0 + 1.0/x(1) ./ (1+4.0*x(2)*t/w^2);
            end
        end
        sum=0.0;
        for i=1:FitRange
            sum=sum+( g(i) - func_autocorrelation_function( x, tau(i)) )^2;
        end
        Q=sum;
    end

end