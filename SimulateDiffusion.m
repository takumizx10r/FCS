clear
close all
dt=2*10^-6;
N_particle=2E3;
MaxFrame=1*10^5;
Resolution=0.2;
noise_level=1.0;
D=400;
mu=0;
sigma=sqrt(2.0*D*dt);
% % % Initial condition
x_g=-1 + (1+1)*rand(2,N_particle);
x_l=x_g;
x0_g=x_g;
rng(0,"twister");
% % % Initial output t=0
iter=0;
MSD(iter+1)=0.0;
Time(iter+1)=(iter-1)*dt;
[N_photon(iter+1),BF]=count_particle(x_l,Resolution);
% % % Make output folders
outputfigurefile=strcat(pwd,'/im_particle/');
mkdir(outputfigurefile);
% % %
for iter=1:MaxFrame
    % % %    Update
    dx=sigma.*randn(2,N_particle)+mu;
    x_g=x_g+dx;
    x_l=x_l+dx;
    % % %    Calculate MSD
    R2=( (x_g(1,:)-x0_g(1,:)).*(x_g(1,:)-x0_g(1,:)) ...
        +(x_g(2,:)-x0_g(2,:)).*(x_g(2,:)-x0_g(2,:)));
    MSD(iter+1)=mean(R2);
    Time(iter+1)=iter*dt;
    % % %
    % % %     Boundary condition - continuous -
    keep=find(x_l(1,:)>1.0);
    x_l(1,keep)=x_l(1,keep)-2.0;
    keep=find(x_l(1,:)<-1.0);
    x_l(1,keep)=x_l(1,keep)+2.0;
    keep=find(x_l(2,:)>1.0);
    x_l(2,keep)=x_l(2,keep)-2.0;
    keep=find(x_l(2,:)<-1.0);
    x_l(2,keep)=x_l(2,keep)+2.0;
    % % %
    % % % Count num particle
    [N_photon(iter+1),BF]=count_particle(x_l,Resolution);
    % % %
    % % % output
%     if rem(iter,10^3)==0
%         scatter(x_l(1,:),x_l(2,:),'k.');
%         hold on
%         scatter(x_l(1,BF),x_l(2,BF),'green.');
%         pbaspect([1 1 1])
%         cir=viscircles([0. 0.], Resolution/2.,'color','b','linestyle','--','linewidth',1);
%         xlim([-1 1]);
%         ylim([-1 1]);
%         ax=gca;
%         xlabel('\itx','FontSize',20)
%         ylabel('\ity','FontSize',20)
%         ax.FontSize=18;
%         axtoolbar('Visible','off');
%         f(fix(iter/10^3))=getframe(gcf);
%         hold off
%         exportgraphics(gcf, ...
%             strcat(outputfigurefile,sprintf('%05d.png',fix(iter/10^3))), ...
%             'Resolution',600)
%     end
    % % %
end

% % % Fitting MSD
F = @(x,xdata)4.0*x(1)*xdata;
x0=D;
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,Time,MSD);
figure
hold on
plot(Time,MSD,'k*','MarkerSize',1)
plot(Time,F(x,Time),'b-')
ax=gca;
xlabel('Time','FontSize',20)
ylabel('MSD','FontSize',20)
ax.FontSize=18;
axtoolbar('Visible','off');
legend('\itMSD_N','\it4Dt','fontsize',18,'location','northwest')
xlim([0 inf]);
ylim([0 inf]);
hold off
exportgraphics(gcf, ...
    strcat(pwd,'/MSD.png'), ...
    'Resolution',600)
% % % % % % % % % % % % % %

% % % Fitting Autocorrelation function
noise=-noise_level + (noise_level+noise_level)*rand(1,numel(N_photon));
N_photon=N_photon+noise;
[FitPara]=func_FCS(dt,N_photon,Resolution,inf);
N_detect=FitPara(1);
D_eff=FitPara(2);
exportgraphics(gcf, ...
    strcat(pwd,'/AutoCoFunc.png'), ...
    'Resolution',600)
% % % % % % % % % % % % % %












% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

function [num_particle,count]=count_particle(x_position,r)
count=find( x_position(1,:).*x_position(1,:) + x_position(2,:).*x_position(2,:) ...
    <= (r*0.50)*(r*0.50) );
num_particle=size(count,2);
end