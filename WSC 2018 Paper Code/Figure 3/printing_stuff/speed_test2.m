close all, clear all, clc
[A,B] =meshgrid(10:80,10:100);
scrX = [A(:),B(:)];
scrX = scrX(scrX(:,1)+14.5 < scrX(:,2),:);
scrXn = [scrX;...
    repmat(scrX(scrX(:,1)-scrX(:,2)==-15,:),[5,1]);
    repmat(scrX(scrX(:,1)==10,:),[5,1]);
    repmat(scrX(scrX(:,1)==80,:),[5,1]);
    repmat(scrX(scrX(:,2)==10,:),[5,1]);
    repmat(scrX(scrX(:,2)==100,:),[5,1])];

num_sim = 50;
K = 35;


warning off
[IDX, C] = kmeans(scrXn,K);
X = round(C);
EWW = sSsimuCRN(X,num_sim);


EWW = sSsimuCRN(X,num_sim);
Xr = repmat(X,[num_sim ,1]);
EWW2 = reshape(EWW,1,[]);
EWW = sSsimu(Xr);
[R3 t2] = RS_L_old(scrX,Xr,EWW,0.05);
figure(1)
subplot(1,2,1)
for k = 1:size(scrX,1)
    rectangle('Position',[scrX(k,1)-R3(k)/2  scrX(k,2)-R3(k)/2 R3(k) R3(k)],'facecolor',[0.2 0.8 0.8],'linestyle','none','Curvature',1)
end

EWW = sSsimuCRN(X,num_sim);
[R2 t2] = RS_L_CRN(scrX,X,EWW,0.05);
subplot(1,2,2)
for k = 1:size(scrX,1)
    rectangle('Position',[scrX(k,1)-R2(k)/2  scrX(k,2)-R2(k)/2 R2(k) R2(k)],'facecolor',[0.8 0.2 0.8],'linestyle','none','Curvature',1)
end
