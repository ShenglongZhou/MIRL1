function [x, Out] = MIRL1(A,b,opts)

% A solver for reweighted L1-minimization model:
% min 0.5||Ax-b||^2 + mu||w.*x||_1
% Written by 05/05/2015, Shenglong Zhou
% yall1 solver is taken from http://yall1.blogs.rice.edu/

% --- Input:
%     A    --- an m x n matrix with m<n;
%     b    --- an m-vector;
%     opts --- a structure with fields:
%              opts.tol   -- tolerance for yall1 solver, defaul one is 1e-4
%              opts.rate  -- for updating the sparsity, defaul one is 1/(log(n/m));
%              opts.k     -- for the given sparsity;

[m,n] = size(A);                  At    = A';
x     = zeros(n,1);               w     = ones(n,1); 
mu    = .01*max(abs(At*b));       a     = 0.2; 
i0    = floor(m/(4*log(n/m)));    eps1  = 1e-10; 
tol   = 1e-4;                     theta = mu*m/n/10;

if log(n/m)<=1; rate=.7;    else rate=1/(log(n/m));  end 
if  isfield(opts,'rate');   rate = opts.rate;        end
if  isfield(opts,'tol');    tol  = opts.tol;         end
if n<2000;  itermax=1000;   else   itermax=100;      end
opts_yall1.tol=tol;
tic;
for iter=1:itermax
        
    x0 = x; 
    w0 = w; 
    opts_yall1.pho     = mu; 
    opts_yall1.weights = w; 
    opts_yall1.x0      = x0;
    x  = yall1(A, b, opts_yall1);               % call yall1 solver to solve the weighted L1-minimization
    
    xx0=x-x0; 
    ErrorTol=sqrt(sum(xx0.*xx0))/max(sqrt(sum(x0.*x0)),1); 
    if  isfield(opts,'iter') && opts.iter==1;
    fprintf(' Iter:  %2d   ErrorTol: %1.2e\n',iter,ErrorTol );
    end

    if ErrorTol<1e-2 | iter==itermax;                          % refinement
        T   = find(abs(x)>=1e-3);  B   = A(:,T);
        x = zeros(n,1);            x(T)= linsolve(B,b);             
        break;
    end                                         
     
     sx = sort(abs(x),'descend');  eps2  = max(1e-3,sx(i0));  
     k  = sparsity(sx);            theta = 1.005*theta;       % update the sparsity 
       
     w  = ModWeight(x,abs(xx0),theta,k,eps2);                      % update the weight
     
     beta=sum(w0.*abs(x))/sum(w.*abs(x));
     if beta>1; mu=a*mu; else mu=beta*mu; end                 % update the penalty parameter

end
Out.iter=iter; Out.time=toc;
%------------------Modifeid weights----------------------%
function w = ModWeight(x,h,theta,k,eps2)
    
    [~,Index]=sort(h,'descend'); w=ones(n,1);      
    if k==0;
        w=1./(abs(x)+eps2);
    else
        w(Index(1:k))  = eps1+theta*sum(h(Index(2:k+1)))/sum(h(Index(1:k)));
        w(Index(k+1:n))= eps1+theta+1./(abs(x(Index(k+1:n)))+eps2);
   end  
end

%------------------Update the Sparsity------------------%
function sp = sparsity(x)
    sumx=sum(x); y=0; sp=0;   
    while y<(rate*sumx)
        sp=sp+1; y=y+x(sp,1); 
    end    
end
end
