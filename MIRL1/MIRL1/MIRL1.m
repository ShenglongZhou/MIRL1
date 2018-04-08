
function [x, Out] = MIRL1(A,b,opts)

% A solver for reweighted L1-minimization model:
% min 0.5||Ax-b||^2 + mu||w.*x||_1
% Written by 05/05/2015, Shenglong Zhou
% yall1 solver is taken from http://yall1.blogs.rice.edu/

% --- Inputs:
%     A    --- an m x n order matrix with m<n;
%     b    --- an m x 1  order vector;
%     opts --- a structure with fields:
%              opts.tol    -- tolerance for yall1 solver, default one: 1e-4;
%              opts.rate   -- for updating the sparsity,  default one: 1/(log(n/m));
%              opts.k      -- for the given sparsity  if it is known in advance;
%              opts.IterOn -- display results in each iteration, default one: 1.
% --- Outputs:
%     x    --- recovered solution, an n x 1  order vector;
%     Out ---  a structure with fields:
%              Out.iter    -- number of total iterations;
%              Out.time    -- total computational time.


[m,n] = size(A);                  At    = A';
x     = zeros(n,1);               w     = ones(n,1); 
mu    = 0.01*max(abs(At*b));      a     = 0.2; 
i0    = floor(m/(4*log(n/m)));    IterOn=1;
tol   = 1e-4;                     theta = mu*m/n/10;

if log(n/m)<=1; rate=.7;    else; rate=1/(log(n/m)); end 
if  isfield(opts,'rate');   rate   = opts.rate;      end
if  isfield(opts,'tol');    tol    = opts.tol;       end

if  isfield(opts,'IterOn'); IterOn = opts.IterOn;    end

if n<2000;  itermax=1000;   else;   itermax=100;     end

opts_yall1.tol = tol;

tic;
for iter=1:itermax
        
    x0 = x; 
    w0 = w; 
    opts_yall1.pho     = mu; 
    opts_yall1.weights = w; 
    opts_yall1.x0      = x0; 
    % call yall1 solver to solve the weighted L1-minimization
    x  = yall1(A, b, opts_yall1);              
    
    xx0=x-x0; 
    ErrorTol=sqrt(sum(xx0.*xx0))/max(sqrt(sum(x0.*x0)),1); 
    if IterOn; fprintf(' Iter:%2d   ErrorTol: %1.2e\n',iter,ErrorTol ); end

    if ErrorTol<1e-2 | iter==itermax                                % refinement
        T   = find(abs(x)>=1e-3);  B   = A(:,T);
        x = zeros(n,1);            x(T)= linsolve(B,b);             
        break;
    end                                         
     
     sx    = sort(abs(x),'descend');  
     eps2  = max(1e-3,sx(i0));  
     if  isfield(opts,'k') 
     k     = opts.k; 
     else
     k     = sparsity(sx,rate);                            % update the sparsity   
     end
     theta = 1.005*theta;       
           
     w     = ModWeight(x,abs(xx0),theta,k,eps2);             % update the weight
     
     beta  = sum(w0.*abs(x))/sum(w.*abs(x));
     if beta>1; mu=a*mu; else; mu=beta*mu; end        % update penalty parameter

end
Out.iter=iter; Out.time=toc;

end

%------------------Modifeid weights----------------------------------------
function w = ModWeight(x,h,theta,k,eps2)
    n     = length(x);
    eps1  = 1e-10; 
    [~,Index]=sort(h,'descend'); w=ones(n,1);      
    if k==0
        w=1./(abs(x)+eps2);
    else
        w(Index(1:k))  = eps1+theta*sum(h(Index(2:k+1)))/sum(h(Index(1:k)));
        w(Index(k+1:n))= eps1+theta+1./(abs(x(Index(k+1:n)))+eps2);
   end  
end

%------------------Update the Sparsity-------------------------------------
function sp = sparsity(x,rate)
    sumx=sum(x); y=0; sp=0;   
    while y<(rate*sumx)
        sp=sp+1; y=y+x(sp,1); 
    end    
end

