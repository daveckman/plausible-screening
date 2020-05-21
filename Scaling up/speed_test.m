% close all, clear all, clc
% scrX = rand(10000,13);
% scrX = 0.5 + 4.5*scrX;
% 
% K = 200;
% 
% warning off
% [IDX, C] = kmeans(scrX,K);
% 
% for c = 1:K
%     dist = sum((repmat(C(c,:),[size(scrX,1),1]) -scrX).^2,2);
%     [~,jfsa] = min(dist);
%     X(c,:) = scrX(jfsa,:);    
% end
% 
% 
% 
% num_sim(1:5) = 1000;
% num_sim(6:K) = 500;
% 
% for k = 1:size(X,1)
%     for l = 1:num_sim(k)
%         Y(k,l) = SAN(X(k,:),10,round(rand*10^8));
%     end
%     
%     for l = (num_sim(k)+1):(max(num_sim))
%         Y(k,l) = nan;
%     end
%     
% end
% 

[R t] = RS_C_A(scrX,X,Y,0.05);

sum(R)
 