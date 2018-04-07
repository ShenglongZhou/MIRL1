clc; clear all;close all;
% This file aims at demonstrating MIRL1 for recovering the compressed
% sensing problems, and be the supportive material for the article titled
% "A Null-Space-Based Weighted l1 Minimization Approach to Compressed Sensing"
% which can be downloaded at:
% http://imaiai.oxfordjournals.org/content/early/2016/02/11/imaiai.iaw002

% Gaussian, Partial DCT, Toeplitz Correlation and Over Sampled Partial DCT 
% measurement matrices will be tested. 
% Different data dimensions can be tested under same example.

% Initialization 
addpath('MIRL1'); 
%m =    64; n  = 4*m;    k =20;
%m =   100; n  = 2000;   k = 10; % This is particularly for 'OverSamDCTMat';
m = 1000; n  = 4*m;     k = floor(0.01*n);

proname  = {'GaussianMat','PartialDCTMat','ToeplitzCorMat','OverSamDCTMat'}; 
problem  = proname{1};
% generate the data (A, b, x_opt)
fprintf('\n Start to generate data... ');
fprintf('\n ----------------------------------- ');
[A,b,x_opt ] =  CSMatrix(problem, m,n,k );
if isequal(b,[]); return; end

% call MIRL1 solver to solve the problem
fprintf('\n Problem name: '); disp(problem);
fprintf(' Sample size:  n=%d, m=%d, k=%d.\n \n Start to run...',n,m,k);
fprintf('\n ----------------------------------- \n');
opts.IterOn = 1;
[x, Out]  = MIRL1(A,b,opts);

% result output
fprintf('\n ----------------------------------- ');
fprintf('\n Relative error  =  %1.3e ', norm(x-x_opt)/norm(x));
fprintf('\n CPU time        =  %1.3f(sec)\n\n', Out.time);
