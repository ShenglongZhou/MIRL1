clc; clear all;close all;
% This file aims at generating recovery success rate.
% Toeplitz  Correlation and Over Sampled Partial DCT type measurement matrices 
% will be tested.
% proname{1} and k0=10:2:40  will generate the plot of recovery success rate
% for Toeplitz Correlation matrix.
% proname{2} and k0=10:2:32  will generate the plot of recovery success rate 
% for Over Sampled Partial DCT matrix.

% Initialization 
addpath('MIRL1'); 
m        = 100; 
n        = 1000; 
k0       = 10:2:36; 
Sample   = 100; 
SuccRate = []; 
proname  = {'ToeplitzCorMat','OverSamDCTMat'}; 
problem  = proname{1};              % change here to test another type matrix
A        = zeros(m,n); 
 
if isequal(problem,'ToeplitzCorMat') % generate the Toeplitz Correlation matrix
   mark  = 0;
   Sig   = zeros(n,n);
   for i = 1:n
   for j = 1:n; Sig(i,j) = (.5)^(abs(i-j));end
   end
   Sig=real(Sig^(1/2));   
else
    mark=1;
end

% Test examples
for j=1:length(k0) 
    rate = 0; 
    k    = k0(j); 
    for p= 1:Sample
        if mark==0;            % generate data for Toeplitz Correlation 
        A     = randn(m,n)*Sig;
        I0    = randperm(n); I=I0(1:k);
        x_opt = zeros(n,1);  
        while nnz(x_opt)~=k; x_opt(I) = randn(k,1); end 
        x_opt = x_opt+0.01*sign(x_opt ); 
        b     = A(:,I)*x_opt(I); 
        else                   % generate data for Over Sampled Partial DCT 
        [A,b,x_opt ] =  CSMatrix(problem, m,n,k ); 
        end
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
