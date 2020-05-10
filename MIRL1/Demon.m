clc; clear all;close all;
% This file aims at demonstrating MIRL1 for recovering the compressed
% sensing problems, and be the supportive material for the article titled
% "A Null-Space-Based Weighted l1 Minimization Approach to Compressed Sensing"
% which can be downloaded at:
% http://imaiai.oxfordjournals.org/content/early/2016/02/11/imaiai.iaw002

% Gaussian, Partial DCT, Toeplitz Correlation and Over Sampled Partial DCT 
% measurement matrices will be tested. 
% Different data dimensions can be tested under same example.

% Initialization 
warning off, addpath('MIRL1');

test = 1;
m    = 1000; n = 4*m;  k = floor(0.01*n);
if test == 4  % This is particularly for 'OverSamDCTMat' 
m    = 100;  n = 2000; k = 10; 
end

proname  = {'GaussianMat','PartialDCTMat','ToeplitzCorMat','OverSamDCTMat'}; 
problem  = proname{test};

% generate the data (A, b, x_opt)
fprintf('\n Start to generate data... ');
fprintf('\n ---------------------------------- ');
[A,b,xopt] = CSMatrix(problem,m,n,k); 

% call MIRL1 solver to solve the problem
fprintf('\n Problem name: %s\n',problem);  
fprintf(' Sample size:  n=%d, m=%d, k=%d\n \n Start to run...',n,m,k);
fprintf('\n ---------------------------------- \n');
opts.IterOn = 1;
%opts.k      = k;
[x, Out]    = MIRL1(A,b,opts);
relerr      = norm(x-xopt)/norm(x);
clear A b xopt

% result output
fprintf('\n ---------------------------------- ');
fprintf('\n Relative error  =  %1.3e ', relerr);
fprintf('\n CPU time        =  %1.3f(sec)', Out.time);
if relerr<1e-2; fprintf('\n Recovery is successful!\n\n');
else;           fprintf('\n Recovery is unsuccessful!\n\n');
end
