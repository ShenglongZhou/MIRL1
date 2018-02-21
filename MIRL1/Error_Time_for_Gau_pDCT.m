clc; clear all;close all;
% This file aims at generating resutls of recovery error and CPU time. 
% Gaussian and  Partial DCT type  measurement matrices, that will be tested.
% proname{1} will generate the plots for Gaussian matrix
% proname{2} will generate the plots Partial DCT matrix
% k=floor(0.01*n), floor(0.05*n) or floor(0.08*n) will generate different
% cases under same example.

% Initialization
addpath('MIRL1'); 
n0      = 1000:1000:5000;  
Smaple  = 10;               
proname = {'GaussianMat','PartialDCTMat'}; 
problem = proname{2};         % change here to test another type matrix
result  = [];

% Test exmaples
for j=1:length(n0)                                                                                                                
    n     = n0(j); 
    m     = floor(n/4); 
    k     = floor(0.01*n);    % k=floor(0.05*n) or k=floor(0.08*n)
    gapxx = 0;
    gapAx = 0;
    time  = 0;    
    for i=1:Smaple            % test Sample=10 examples for each dimension
    [A,b,x_opt ] = CSMatrix(problem, m,n,k );            
    t0    = cputime; 
    x     = MIRL1(A,b,[]); 
    t     = cputime-t0;             
    gapxx = gapxx + norm(x_opt-x,'fro') ;
    gapAx = gapAx + norm(A*x-b,'fro');
    time  = time    + t;
    end   
    result=[result   [ gapxx; gapAx ; time ]/Smaple]
end

% Graph design
if k< floor(0.08*n) 
for i=1:3
    subplot(1,3,i)
    plot(n0,result(i,:),'r*-');hold on
    set(gca,'FontName','Times','FontSize',8)
    if     i==1; ylabel('||x-x_{orig}||','FontSize',10 )
    elseif i==2; ylabel('||\Phix-b||','FontSize',10) 
    else         ylabel('Time','FontSize',10) 
    end
    xlabel('N','FontName','Times','FontSize',8)
    axis([min(n0) max(n0) 0 max(result(i,:))])
end
end

