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

if isequal(pbnm{test},'ToeplitzCorMat')  % generate the Toeplitz Correlation matrix
   Sig   = zeros(n,n);
   for i = 1:n; Sig(i,:)=(0.5).^(abs(i-(1:n))); end
   Sig   = real(Sig^(1/2));   
end

% Test examples
for j = 1:length(k0) 
    rate  = 0; 
    k     = k0(j); 
    for p = 1:Smpl
        switch pbnm{test}
            case 'ToeplitzCorMat' % generate data for Toeplitz Correlation                  
                A    = randn(m,n)*Sig;
                I0   = randperm(n); 
                I    = I0(1:k);
                xopt = zeros(n,1);  
                while nnz(xopt)~=k; xopt(I) = randn(k,1); end 
                xopt = xopt + 0.01*sign(xopt); 
                b    = A(:,I)*xopt(I); 
            case 'OverSamDCTMat'  % generate data for Over Sampled Partial DCT 
                [A,b,xopt ] =  CSMatrix(pbnm{test}, m,n,k ); 
        end
        opts.IterOn  = 0;
        x            = MIRL1(A,b,opts);    
        rate         = rate+ (norm(x-xopt)/norm(x)<1e-2);
    end
    ScRt=[ScRt rate/Smpl]; clc; ScRt
end

% Graph design
figure, plot(k0,ScRt(1,:),'r*-');
ylabel('Success rate'), xlabel('Sparsity')
title(pbnm{test}), grid on
