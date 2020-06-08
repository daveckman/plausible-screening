% close all, clear all, clc
% [A,B] =meshgrid(10:80,10:100);
% scrX = [A(:),B(:)];
% scrX = scrX(scrX(:,1)+14.5 < scrX(:,2),:);
% scrXn = [scrX;...
%     repmat(scrX(scrX(:,1)-scrX(:,2)==-15,:),[5,1]);
%     repmat(scrX(scrX(:,1)==10,:),[5,1]);
%     repmat(scrX(scrX(:,1)==80,:),[5,1]);
%     repmat(scrX(scrX(:,2)==10,:),[5,1]);
%     repmat(scrX(scrX(:,2)==100,:),[5,1])];
% 
% K = 60;
% num_sim = 5*ones(K,1);
% 
% 
% warning off
% [IDX, C] = kmeans(scrXn,K);
% X = round(C);
% 
% 
% num_sim(1:K) = 100;
% for k = 1:size(X,1)
%     for l = 1:num_sim(k)
%         Y(k,l) = sSsimu( X(k,:));
%     end
%     for l = (num_sim(k)+1):(max(num_sim))
%         Y(k,l) = nan;
%     end
% end
% 
% save('done_simu_26')
% disp('done simu')
close all, clear all, clc
load('done_simu_26')

figure(2)
% profile on
tic
R = RS_C_A(scrX,X,Y,0.05);
toc
%profile viewer
subplot(1,2,2)
for k = 1:size(scrX,1)
    rectangle('Position',[scrX(k,1)-R(k)/2  scrX(k,2)-R(k)/2 R(k) R(k)],'facecolor',[0.8 0.2 0.8],'linestyle','none','Curvature',1)
end
%profile viewer


figure(2)
% profile on
tic
R = RS_CL_A(scrX,X,Y,0.05);
toc
%profile viewer
subplot(1,2,1)
for k = 1:size(scrX,1)
    rectangle('Position',[scrX(k,1)-R(k)/2  scrX(k,2)-R(k)/2 R(k) R(k)],'facecolor',[0.8 0.2 0.8],'linestyle','none','Curvature',1)
end
%profile viewer
% 
%  %profile on
% tic
% R2 = RS_L_old(scrX,X,Y,0.05);
% toc
% %profile viewer
% subplot(1,2,1)
% for k = 1:size(scrX,1)
%     rectangle('Position',[scrX(k,1)-R2(k)/2  scrX(k,2)-R2(k)/2 R2(k) R2(k)],'facecolor',[0.2 0.8 0.8],'linestyle','none','Curvature',1)
%     rectangle('Position',[scrX(k,1)-R(k)/2  scrX(k,2)-R(k)/2 R(k) R(k)],'edgecolor',[0.8 0.2 0.8],'Curvature',1)
% end
