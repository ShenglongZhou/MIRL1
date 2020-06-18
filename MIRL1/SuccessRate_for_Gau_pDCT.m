clc; clear all;close all;
% This file aims at generating recovery success rate.
% Gaussian and  Partial DCT type measurement matrices, that will be tested.
% pbnm{1} will generate the plot of success rate for Gaussian matrix
% pbnm{2} will generate the plot of success rate for Partial DCT matrix

% Initialization
addpath('MIRL1'); 
m      = 64; 
n      = 256; 
k0     = 10:2:40; 
Smpl   = 100; 
test   = 1; % test =1, 2 
pbnm   = {'GaussianMat','PartialDCTMat'}; 
ScRt   = [];

% Test examples
for j     = 1:length(k0) 
    rate  = 0;
    for p = 1: Smpl        
        [A,b,x_opt] = CSMatrix(pbnm{test}, m,n,k0(j) ); 
        opts.IterOn = 0;
        x           = MIRL1(A,b,opts);    
        rate        = rate+( norm(x-x_opt)/norm(x)<1e-2 ); 
    end       
    ScRt        = [ScRt rate/Smpl]; clc; ScRt
end

% Graph design
figure, plot(k0,ScRt,'r*-');
ylabel('Success rate'), xlabel('Sparsity')
title(pbnm{test}), grid on
