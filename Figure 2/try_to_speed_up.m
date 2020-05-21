% 
% scrX =[1020:1119]';
% Kt = length(scrX);
% X =scrX(end:-5:1);
% N = 500;
% Xs =repmat(X,N/length(X),1);
% Y = mmk(Xs);
% 
% save('blah_07_12_2018')
 clear all, close all, clc
 
 load ('blah_07_12_2018')
 profile on
 
tic
[SC V] = RS_C(scrX,Xs,Y,0.01);
toc
tic
[SC V2] = RS_C_old(scrX,Xs,Y,0.01); 
toc
[V' V2]
% 
% [V V2 V<(V2+0.1)]