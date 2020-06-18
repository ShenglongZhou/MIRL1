clc; clear all;close all;
% This file aims at generating resutls of recovery error and CPU time. 
% Gaussian and  Partial DCT type measurement matrices, that will be tested.
% pbnm{1} generates the plots for Gaussian matrix
% pbnm{2} generates the plots for Partial DCT matrix
% k=0.01*n, 0.05*n or 0.08*n generate different cases under same example.

% Initialization
addpath('MIRL1'); 
n0      = 1000:1000:5000;  
Smpl    = 10;             
test    = 1;     % test=1, 2
rate    = 0.01;  % rate=0.01, 0.05, 0.08
pbnm    = {'GaussianMat','PartialDCTMat'}; 
result  = [];

% Test exmaples
for j = 1:length(n0)                                                                                                                
    n     = n0(j); 
    m     = ceil(n/4); 
    k     = ceil(rate*n);    
    dxx   = 0;
    dAx   = 0;
    time  = 0;    
    for i = 1:Smpl % Smpl=10 examples for each dimension
    [A,b,xopt] = CSMatrix(pbnm{test},m,n,k ); 
    opts.IterOn= 0;
    [x,out]    = MIRL1(A,b,opts);           
    dxx        = dxx  + norm(xopt-x,'fro') ;
    dAx        = dAx  + norm(A*x-b,'fro');
    time       = time + out.time;
    end   
    result     = [result [dxx;dAx;time]/Smpl]
end

% Graph design
figure
ylab = {'||x-x_{opt}||','||\Phix-b||','Time'};
for i= 1:3
    subplot(1,3,i)
    plot(n0,result(i,:),'ro-'),   hold on    
    ylabel(ylab{i}), xlabel('N'), grid on
    axis([min(n0) max(n0) 0 max(result(i,:))])
end
