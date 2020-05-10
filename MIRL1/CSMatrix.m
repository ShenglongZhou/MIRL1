function [A,b,x_opt ] = CSMatrix(problemname,m,n,k )
% This file aims at generating data of four examples, Gaussian, Partial DCT, 
% Toeplitz Correlation and Over Sampled Partial DCT type measurement matrices,
% including measurement matrices A, the gound truth sparse solution x_opt and
% the observation vector b, namely, A*x_opt=b.

A = zeros(m,n);

switch problemname     
    case 'GaussianMat'        
        A     = sqrt(1/m)*randn(m,n);
        I0    = randperm(n); 
        I     = I0(1:k);       
    case 'PartialDCTMat'
        r     = rand(m,1);  
        for i = 1:m 
        A(i,:)= sqrt(1/m)*cos(2*pi*r(i)*(0:(n-1)));
        end
        I0    = randperm(n); 
        I     = I0(1:k);        
    case 'ToeplitzCorMat'
        Sig   = zeros(n);
        for i = 1:n; Sig(i,:) = (.5).^(abs(i-(1:n))); end    
        A     = mvnrnd(zeros(n,1),Sig,m); 
        I0    = randperm(n); I=I0(1:k);           
    case 'OverSamDCTMat'
        F     = 10;
        if n-2*F*(k-1)<=0
           fprintf('\n Please re-input k which must be <= n/20 \n'); 
           return;
        end
        r     = rand(m,1); 
        for i = 1:m 
        A(i,:)= sqrt(1/m)*cos(2*pi*r(i)*(0:(n-1))/F);
        end
        % random sampling k integers from 1--n with spacing at least 2F
        I     = randsample(n-2*F*(k-1),k);
        I     = sort(I);
        I     = I + (0:k-1)'*2*F;         
    otherwise
        disp('Please input a problem name !!!');        
end

x_opt    = zeros(n,1);  
while nnz(x_opt)~=k 
x_opt(I) = randn(k,1); 
end
x_opt    = x_opt+0.01*sign(x_opt ); 
b        = A(:,I)*x_opt(I); 
end
