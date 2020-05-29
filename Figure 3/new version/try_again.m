clear all, close all, clc
[A,B] =meshgrid(10:80,10:100);

scrX = [A(:),B(:)];
scrX = scrX(scrX(:,1)+14.5 < scrX(:,2),:);
scrXn = [scrX;...
    repmat(scrX(scrX(:,1)-scrX(:,2)==-15,:),[5,1]);
    repmat(scrX(scrX(:,1)==10,:),[5,1]);
    repmat(scrX(scrX(:,1)==80,:),[5,1]);
    repmat(scrX(scrX(:,2)==10,:),[5,1]);
    repmat(scrX(scrX(:,2)==100,:),[5,1])];

num_sim = 10;
K = 25;

[IDX, C] = kmeans(scrXn,K);
X = round(C);
plot(X(:,1),X(:,2),'ro')
hold on
plot(scrX(:,1),scrX(:,2),'r.')

Xr = repmat(X,[num_sim ,1]);

mu = simu_oracle(scrX);
cact = 0;
for k = 1:size(scrX,1)
        cprop = abs(mu(k)-mu)./(sqrt(sum((repmat(scrX(k,:),size(scrX,1),1)-scrX).^2,2))');
        if max(cprop) >= cact 
            cact = max(cprop);
        end
end
[~,istar] = min(mu);
cact

M = 100;

L_save = zeros(M,length(scrX));

parfor_progress(100); 
parfor c = 1:M
L_save(c,:)  = RS_L(scrX,Xr,sSsimu(Xr),0.05,3);
        parfor_progress; 
end
parfor_progress(0); 

save('Data_02_27_2018','L_save','num_sim','scrX','X','K','istar','cact','mu','M')
        
clear all, close all, clc
load('Data_02_27_2018')

C = sum(L_save,1)'/M;

for k = 1:size(scrX,1)
    
rectangle('Position',[scrX(k,1)-C(k)/2  scrX(k,2)-C(k)/2 C(k) C(k)],'facecolor',[0.2 0.8 0.8],'linestyle','none','Curvature',1)

end
axis square
hold on

plot(X(:,1),X(:,2),'kd','markerfacecolor','k')
plot(scrX(istar,1),scrX(istar,2),'r*')
plot(scrX(istar,1),scrX(istar,2),'ro')

vch = convhull(scrX(:,1),scrX(:,2));
plot(scrX(vch,1),scrX(vch,2),'k-')
xlim([min(scrX(:,1))-5,max(scrX(:,1))+5])
ylim([min(scrX(:,2))-5,max(scrX(:,2))+5])
xlabel('Reorder quantity','interpreter','latex')
ylabel('Order-to quantity','interpreter','latex')

print_correctly('sS_numerical')

