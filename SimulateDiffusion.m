clear
close all
dt=2*10^-6;
N_particle=10^3;
MaxFrame=1*10^4;
Resolution=0.2;
D=200;
mu=0;
sigma=sqrt(2.0*D*dt);
% % % Initial condition
x_g=-1 + (1+1)*rand(2,N_particle);
x_l=x_g;
x0_g=x_g;
% % % Initial output t=0
iter=0;
MSD(iter+1)=0.0;
Time(iter+1)=(iter-1)*dt;
N_photon(iter+1)=count_particle(x_l,Resolution);
% % %

rng(0,"twister");
for iter=1:MaxFrame
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
    N_photon(iter+1)=count_particle(x_l,Resolution);
    % % %
    % % % output
    if rem(iter,1000)==0
        scatter(x_l(1,:),x_l(2,:),'k.');
        pbaspect([1 1 1])
        hold on
        cir=viscircles([0. 0.], Resolution/2.,'color','b','linestyle','--');
        xlim([-1 1]);
        ylim([-1 1]);
        f(fix(iter/100))=getframe(gcf);
        hold off
    end
    % % %
end

% % % Fitting MSD
F = @(x,xdata)4.0*x(1)*xdata;
x0=D;
[x,resnorm,~,exitflag,output] = lsqcurvefit(F,x0,Time,MSD);
figure
hold on
plot(Time,MSD)
plot(Time,F(x,Time))
xlabel  'Time';
ylabel  'MSD';
hold off
% % % 

% % % Fitting Autocorrelation function
[FitPara]=func_FCS(dt,N_photon);
N_detect=FitPara(1);
D_eff=(Resolution/2.0)*(Resolution/2.0) / (4.0*FitPara(2));
% % % 


function [num_particle]=count_particle(x_position,r)
    count=find( x_position(1,:).*x_position(1,:) + x_position(2,:).*x_position(2,:) ...
        <= (r/2.0)*(r/2.0) );
    num_particle=size(count,2);
end