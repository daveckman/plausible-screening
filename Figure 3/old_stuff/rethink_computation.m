
[A,B] =meshgrid(10:80,10:60);
scrX = [A(:),B(:)];
scrX = scrX(scrX(:,1)+14.5 < scrX(:,2),:);
scrXn = [scrX;...
    repmat(scrX(scrX(:,1)-scrX(:,2)==-15,:),[5,1]);
    repmat(scrX(scrX(:,1)==10,:),[5,1]);
    repmat(scrX(scrX(:,1)==80,:),[5,1]);
    repmat(scrX(scrX(:,2)==10,:),[5,1]);
    repmat(scrX(scrX(:,2)==100,:),[5,1])];

num_sim = 200;
K = 50;


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


warning off
[IDX, C] = kmeans(scrXn,K);
X = round(C);

Xr = repmat(X,[num_sim ,1]);

EWW = sSsimu(Xr);
 profile on
%  

close all, clc
tic
[R t] = RS_L_RT2(scrX,Xr,EWW,0.05);
toc
profile viewer
figure(1)
subplot(1,2,1)
for k = 1:size(scrX,1)
    rectangle('Position',[scrX(k,1)-R(k)/2  scrX(k,2)-R(k)/2 R(k) R(k)],'facecolor',[0.2 0.8 0.8],'linestyle','none','Curvature',1)
end
% tic
% [R2 t2] = RS_L_old(scrX,Xr,EWW,0.05);
% toc
% figure(1)
% subplot(1,2,2)
% for k = 1:size(scrX,1)
%     rectangle('Position',[scrX(k,1)-R2(k)/2  scrX(k,2)-R2(k)/2 R2(k) R2(k)],'facecolor',[0.8 0.2 0.8],'linestyle','none','Curvature',1)
% end

