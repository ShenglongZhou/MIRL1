clc; clear all;close all;
% This file aims at generating recovery success rate.
% Gaussian and  Partial DCT type measurement matrices, that will be tested.
% proname{1} will generate the plot of success rate for Gaussian matrix
% proname{2} will generate the plot of success rate for Partial DCT matrix

% Initialization
addpath('MIRL1'); 
m        = 64; 
n        = 256; 
k0       = 10:2:40; 
Sample   = 100; 
SuccRate = [];
proname  = {'GaussianMat','PartialDCTMat'}; 
problem  = proname{1}; % change here to test another type matrix

% Test examples
for j=1:length(k0) 
    rate  = 0;
    for p = 1:Sample        
    [A,b,x_opt ] = CSMatrix(problem, m,n,k0(j) ); 
    opts.IterOn  = 0;
    x            = MIRL1(A,b,opts);    
    if norm(x-x_opt)/norm(x)<1e-2; rate=rate+1;end 
    end       
    SuccRate=[SuccRate rate/Sample]; clc; SuccRate
end

% Graph design
plot(k0,SuccRate(1,:),'r*-');
set(gca,'FontName','Times','FontSize',10)
ylabel('Success rate') 
xlabel('Sparsity')
title([num2str(problem)])
