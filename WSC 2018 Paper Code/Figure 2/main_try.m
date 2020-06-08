clear all, close all, clc
scrX =[1020:1119]';
Kt = length(scrX);
mu = oracle_queue(scrX)';

[~,istar] = min(mu);
X =scrX(1:5:end);
Y = mmk(X);
plot(X,Y,'.');
hold on
plot(scrX,mu,'k-')

Nv = [100 200 500 1000];
num_sim = 100;
SC = zeros(length(Nv),length(scrX),num_sim);
SL = zeros(length(Nv),length(scrX),num_sim);
SSS = zeros(length(Nv),length(scrX),num_sim);
tvals = zeros(size(Nv));
for s = 1:length(Nv)
    N = Nv(s);
    if s > 1
        fprintf('\n Elaspsed Time (min) = %0.1f',sum(tvals(1:(s-1))))
        fprintf('\n Estimated Total Time (min) = %0.1f',mean((tvals(1:(s-1))./Nv(1:(s-1))))*sum(Nv(1:end)))
        fprintf('\n Estimated Time Left (min) = %0.1f \n',mean((tvals(1:(s-1))./Nv(1:(s-1))))*sum(Nv(s:end)))
    end
    fprintf('\n \n Case N =  %d \n',N);
    tic
    lcv = ones(1,num_sim);
    parfor_progress(num_sim); 
    parfor k = 1:num_sim
        Xs =repmat(X,N/length(X),1);
        Y = mmk(Xs);
        SC(s,:,k) = SC(s,:,k)+RS_C(scrX,Xs,Y,0.05)';
        SL(s,:,k) = SL(s,:,k)+RS_L(scrX,Xs,Y,0.05,0.03)';
        
        X2 =repmat(scrX,N/length(scrX),1);
        Y2 = mmk(X2);
        SSS(s,:,k) = SSS(s,:,k)+RS_SS(scrX,X2,Y2,0.05)';
        parfor_progress; 
    end
    parfor_progress(0);
    tvals(s) = toc/60;
end

save('Data_03_13_2018','SC','SSS','SL','s','num_sim','scrX','X','mu','istar','Nv')

clear all, close all

load('Data_03_13_2018')
s = 3
figure(1)
set(gca,'Position',[0,0,600,650])
for slcv = 1:s
subplot(3,s, slcv)
bar(scrX,sum(SSS(slcv,:,:),3)/num_sim,1,'linestyle','none','facecolor',[0.6 0.6 0.6])
hold on
plot(scrX,0.95*ones(size(scrX)),'k-')
plot(scrX(istar)*ones(2,1),[0,sum(SSS(slcv,istar,:),3)/num_sim],'r:','linewidth',2)
xlim([1020,1120])
set(gca,'xtick',[1025:25:1125],'ytick',[0:0.2:1])
if slcv ==1
ylabel({'Subset Selection','frequency'},'interpreter','latex')
else
    set(gca,'ytick',[])
end

title(['$N = ' num2str(Nv(slcv)) '$'],'interpreter','latex')


subplot(3, s, slcv+s)
bar(scrX,sum(SL(slcv,:,:),3)/num_sim,1,'linestyle','none','facecolor',[0 0.6 0.6])
hold on
plot(scrX,0.95*ones(size(scrX)),'k-')
plot(scrX(istar)*ones(2,1),[0,sum(SL(slcv,istar,:),3)/num_sim],'r:','linewidth',2)
xlim([1020,1120])
set(gca,'xtick',[1025:25:1125],'ytick',[0:0.2:1])
if slcv ==1
ylabel({'PO, Lipschitz','frequency'},'interpreter','latex')
else
    set(gca,'ytick',[])
end


subplot(3,s,2*s+slcv)
bar(scrX,sum(SC(slcv,:,:),3)/num_sim,1,'linestyle','none','facecolor',[0.6 0 0.6])
hold on
plot(scrX,0.95*ones(size(scrX)),'k-')
plot(scrX(istar)*ones(2,1),[0,sum(SC(slcv,istar,:),3)/num_sim],'r:','linewidth',2)
xlim([1020,1120])
set(gca,'xtick',[1025:25:1125],'ytick',[0:0.2:1])
if slcv ==1
ylabel({'PO, convex','frequency'},'interpreter','latex')
else
    set(gca,'ytick',[])
end
end

print_correctly('mmc_simulation_results')
