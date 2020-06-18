clc; clear all;close all;
% This file aims at generating recovery success rate.
% Toeplitz Correlation and Over Sampled Partial DCT type measurement matrices will be tested.
% pbnm{1} and k0=10:2:40 generate recovery success rate for Toeplitz Correlation matrix.
% pbnm{2} and k0=10:2:32 generate recovery success rate for Over Sampled Partial DCT matrix.

% Initialization 
addpath('MIRL1'); 
m    = 100; 
n    = 1000; 
k0   = 10:2:36; 
Smpl = 100; 
test = 1; 
pbnm = {'ToeplitzCorMat','OverSamDCTMat'}; 
A    = zeros(m,n); 
ScRt = [];  

% Test examples
for j = 1:length(k0) 
    rate  = 0; 
    k     = k0(j); 
    for p = 1:Smpl
        [A,b,xopt ]  = CSMatrix(pbnm{test},m,n,k ); 
        opts.IterOn  = 0;
        opts.k       = k;
        x            = MIRL1(A,b,opts);    
        rate         = rate + (norm(x-xopt)/norm(x)<1e-2);
    end
    ScRt=[ScRt rate/Smpl]; clc; ScRt
end

% Graph design
figure, plot(k0,ScRt(1,:),'r*-');
ylabel('Success rate'), xlabel('Sparsity')
title(pbnm{test}), grid on
