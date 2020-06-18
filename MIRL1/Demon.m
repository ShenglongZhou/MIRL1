clc; clear all;close all;
% This file aims at demonstrating MIRL1 for recovering the compressed
% sensing problems, and be the supportive material for the article titled
% "A Null-Space-Based Weighted l1 Minimization Approach to Compressed Sensing"
% which is available at:
% http://imaiai.oxfordjournals.org/content/early/2016/02/11/imaiai.iaw002

% Gaussian, Partial DCT, Toeplitz Correlation and Over Sampled Partial DCT 
% measurement matrices will be tested. 
% Different data dimensions can be tested under same example.

% Initialization 
warning off, addpath('MIRL1');

test  = 1;
switch test
  case {1,2,3}; m=1e3; n=4*m;  k=ceil(0.01*n);
  case 4;       m=100; n=2000; k=10; %For 'OverSamDCTMat'     
end
proname  = {'GaussianMat',   'PartialDCTMat',...
            'ToeplitzCorMat','OverSamDCTMat'}; 
problem  = proname{test};

% generate (A, b, x_opt)
fprintf(' Generate data of %s...',problem); 
[A,b,xopt] = CSMatrix(problem,m,n,k); 

% call MIRL1 solver 
opts.IterOn = 1;
[x, Out]    = MIRL1(A,b,opts);
relerr      = norm(x-xopt)/norm(x);
 
fprintf(' CPU time:         %.3f(sec)\n', Out.time);
fprintf(' Objective:        %5.3e\n', Out.obj);
fprintf(' Relative error:   %5.3e\n', relerr);
fprintf(' Sample size:      m=%d,n=%d,k=%d\n',m,n,k);
if relerr<1e-2; fprintf(' Recovery is successful!\n\n');
else;           fprintf(' Recovery is unsuccessful!\n\n');
end
