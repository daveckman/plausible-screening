close all, clear all, clc

x1 = -1;
x2 = 1.45;

muhat1 = -1;
muhat2 = 1;

x0 = 0;
[p L] = quadprog([1 0 0 0 0;
    0 1 0 0 0;
    0 0 10^(-9) 0 0;
    0 0 0 10^(-9) 0;
    0 0 0 0 10^(-12)],-[muhat1;
    muhat2;
    0;
    0;
    0],[-1 1 0 -(x2-x1) 0 ;
    1 -1 -(x1-x2) 0 0;
    1 0 -(x1-x0) 0 -1;
    0 1 0 -(x2-x0) -1;
    0 -1 0 0 1;
    -1 0 0 0 1],[0;
    0;
    0;
    0;
    0;
    0]);

m =@(x) max(p(1)+p(3)*(x-x1),p(2)+p(4)*(x-x2));
figure(1)
set(gca,'position',[0 0 800 200])
xtry = -1.5:0.01:3.5
subplot(1,2,1)
plot(xtry,muhat2*ones(size(xtry)),'k:','linewidth',0.5,'color',[0.7 0.7 0.7]);
hold on
plot(xtry,muhat1*ones(size(xtry)),'k:','linewidth',0.5,'color',[0.7 0.7 0.7]);
plot(x0*[1 1],[-10,10],'k:','linewidth',0.5,'color',[0.7 0.7 0.7]);
plot(x1*[1 1],[-10,10],'k:','linewidth',0.5,'color',[0.7 0.7 0.7]);
plot(x2*[1 1],[-10,10],'k:','linewidth',0.5,'color',[0.7 0.7 0.7]);
plot(x1,muhat1,'ko','markerfacecolor',[0.6 0.6 0.6])
plot(x2,muhat2,'ko','markerfacecolor',[0.6 0.6 0.6])
% 
%plot(xtry,m(xtry),'k--','linewidth',1.5);
R1 = min(m(xtry));
R2 = max(m(xtry));
ylim([R1-0.05*(R2-R1), R2-0.3*(R2-R1)])
diff_y = diff([R1-0.05*(R2-R1), R2-0.3*(R2-R1)]);
diff_x = range(xtry);
 
%text(0.4+0.1,m(0.4)-0.1,'convex','rotation',atand(diff_x/diff_y*(m(0.4)-m(0.39))/0.01),...
%    'horizontalalignment','center')

hold on

[p L] = quadprog([1 0 0;
    0 1 0 ;
    0 0 10^(-4)],-[muhat1;
    muhat2;
    0;],[-1 1 0  ;
    1 -1 0;
    1 0 -1;
    0 1 -1;
    0 -1 1;
    -1 0 1],[2*abs(x2-x1);
    2*abs(x2-x1);
    2*abs(x0-x1);
    2*abs(x0-x2);
    0;
    0]);

mp =@(x) max(max(p(1)-2*abs(x-x1),p(2)-2*abs(x-x2)),p(3)-2*abs(x-x0));
mm =@(x) min(min(p(1)+2*abs(x-x1),p(2)+2*abs(x-x2)),p(3)+2*abs(x-x0));
hold on
plot(xtry,0.5*mp(xtry)+0.5*mm(xtry),'k-','linewidth',1.5);
xlim([min(xtry) max(xtry)])
set(gca,'TickLabelInterpreter', 'latex','FontName','Times','fontsize',12)
set(gca,'xtick',[x1 x0 x2],'xticklabel',{'$x_1$','$x_0$','$x_2$'})
set(gca,'ytick',[muhat1 muhat2],'yticklabel',{'$\hat{\mu}_1$','$\hat{\mu}_2$'})
text(0.75-0.2,0.5*mp(0.75)+0.5*mm(0.75)+0.1,'Lipschitz','rotation',atand(diff_x/diff_y*(0.5*mp(0.75)+0.5*mm(0.75)-0.5*mp(0.74)-0.5*mm(0.74))/0.01),...
    'horizontalalignment','center')
axis square
ylabel('simulation response')
xlabel('decision variable')


x0 = 2;
[p L] = quadprog([1 0 0 0 0;
    0 1 0 0 0;
    0 0 10^(-9) 0 0;
    0 0 0 10^(-9) 0;
    0 0 0 0 10^(-12)],-[muhat1;
    muhat2;
    0;
    0;
    0],[-1 1 0 -(x2-x1) 0 ;
    1 -1 -(x1-x2) 0 0;
    1 0 -(x1-x0) 0 -1;
    0 1 0 -(x2-x0) -1;
    0 -1 0 0 1;
    -1 0 0 0 1],[0;
    0;
    0;
    0;
    0;
    0]);

m =@(x) max(p(1)+p(3)*(x-x1),p(2)+p(4)*(x-x2));

subplot(1,2,2)
plot(xtry,muhat2*ones(size(xtry)),'k:','linewidth',0.5,'color',[0.7 0.7 0.7]);
hold on
plot(xtry,muhat1*ones(size(xtry)),'k:','linewidth',0.5,'color',[0.7 0.7 0.7]);
plot(x0*[1 1],[-10,10],'k:','linewidth',0.5,'color',[0.7 0.7 0.7]);
plot(x1*[1 1],[-10,10],'k:','linewidth',0.5,'color',[0.7 0.7 0.7]);
plot(x2*[1 1],[-10,10],'k:','linewidth',0.5,'color',[0.7 0.7 0.7]);
plot(x1,muhat1,'ko','markerfacecolor',[0.6 0.6 0.6])
plot(x2,muhat2,'ko','markerfacecolor',[0.6 0.6 0.6])
%text(2.5,m(2)+0.25,'convex','rotation',0,'horizontalalignment','center')


%drawbrace([x1-0.1 muhat1], [x1-0.1 p(1)])
%drawbrace([x2-0.1 p(2)],[x2-0.1 muhat2])

%plot(xtry,m(xtry),'k--','linewidth',1.5);
ylim([R1-0.05*(R2-R1), R2-0.3*(R2-R1)])
hold on

[p L] = quadprog([1 0 0;
    0 1 0 ;
    0 0 10^(-12)],-[muhat1;
    muhat2;
    0;],[-1 1 0  ;
    1 -1 0;
    1 0 -1;
    0 1 -1;
    0 -1 1;
    -1 0 1],[2*abs(x2-x1);
    2*abs(x2-x1);
    2*abs(x0-x1);
    2*abs(x0-x2);
    0;
    0]);

mp =@(x) max(max(p(1)-2*abs(x-x1),p(2)-2*abs(x-x2)),p(3)-2*abs(x-x0));
mm =@(x) min(min(p(1)+2*abs(x-x1),p(2)+2*abs(x-x2)),p(3)+2*abs(x-x0));
hold on
plot(xtry,0.5*mp(xtry)+0.5*mm(xtry),'k-','linewidth',1.5);
xlim([min(xtry) max(xtry)])
set(gca,'TickLabelInterpreter', 'latex','FontName','Times','fontsize',12)
set(gca,'xtick',[x1  x2 x0],'xticklabel',{'$x_1$','$x_2$','$x_0$'})
set(gca,'ytick',[muhat1 muhat2],'yticklabel',{'$\hat{\mu}_1$','$\hat{\mu}_2$'})
text(2.7,0.5*mp(3)+0.5*mm(3)-0.2,'Lipschitz','rotation',0,...
    'horizontalalignment','center')
% 
% drawbrace([x2+0.05 muhat2],[x2+0.05 p(2)])
% drawbrace([x1+0.1 p(1)],[x1+0.1 muhat1])
% text(x1+0.25,(p(1)+muhat1)/2+0.025,'$|m(x_1) - \widehat{\mu}_1|$','fontsize',7.5,'color','red','Interpreter','latex')
% text(x2+0.21,(p(2)+muhat2)/2+0.025,'$|m(x_2) - \widehat{\mu}_2|$','fontsize',7.5,'color','red','Interpreter','latex')

axis square
xlabel('decision variable')

p(3)
(muhat1+muhat2-2*abs(x0-x2))/2

2*L+muhat1.^2+muhat2.^2
1/2*(muhat1-muhat2+2*abs(x0-x2)).^2

print_correctly('minimizing_functions')