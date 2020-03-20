
function [x, Out] = MIRL1(A,b,opts)

% A solver for reweighted L1-minimization model:
%
%    min 0.5||Ax-b||^2 + mu||w.*x||_1
%
% Written by 05/05/2015, Shenglong Zhou
%
% Note: yall1 solver is taken from http://yall1.blogs.rice.edu/

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
if nargin<2
   error('Inputs are not enough')
elseif nargin==2
   opts=[];
end

[m,n]       = size(A);                 
[Itmax,rate,tol,IterOn,mu,i0,theta]...
            = Get_parameters(m,n,A,b,opts);     
x           = zeros(n,1);               
w           = ones(n,1);
opts_ya.tol = tol;
t_start     = tic;

for iter=1:Itmax
        
    x0 = x; 
    w0 = w; 
    opts_ya.pho     = mu; 
    opts_ya.weights = w; 
    opts_ya.x0      = x0; 
    
    % call yall1 solver to solve the weighted L1-minimization
    x     = yall1(A, b, opts_ya);                  
    dx    = x-x0; 
    Error = sqrt(sum(dx.*dx))/max(sqrt(sum(x0.*x0)),1); 
    if IterOn; fprintf(' Iter:%3d     Error: %1.2e\n',iter,Error ); end

    if Error<1e-2 | iter==Itmax                                % refinement
        T        = find(abs(x)>=1e-3);  
        B        = A(:,T);
        x        = zeros(n,1);            
        x(T)     = linsolve(B,b);
        Out.iter = iter; 
        Out.time = toc(t_start);
        return;
    end    
    
    sx    = sort(abs(x),'descend');  
    eps2  = max(1e-3,sx(i0));  
    if isfield(opts,'k'); k = opts.k;                % update the sparsity
    else;                 k = sparsity(sx,rate);                              
    end
     
    theta = 1.005*theta;                  
    w     = ModWeight(x,abs(dx),theta,k,eps2);        % update the weight    
    beta  = sum(w0.*abs(x))/sum(w.*abs(x));
     
    if beta>1; mu = 0.2*mu;                     % update penalty parameter
    else ;     mu = beta*mu; 
    end        

end
end
%------------------Set Parameters----------------------------------------
function [itmax,rate,tol,IterOn,mu,i0,theta] = Get_parameters(m,n,A,b,opts)

if n<2000;      itmax=1000; else;  itmax = 100;          end
if log(n/m)<=1; rate=.7;    else;  rate  = 1/(log(n/m)); end 
if isfield(opts,'rate');    rate   = opts.rate;          end
if isfield(opts,'tol');     tol    = opts.tol;    else;  tol    = 1e-4; end
if isfield(opts,'IterOn');  IterOn = opts.IterOn; else;  IterOn = 1;    end 
mu    = 0.01*max(abs(A'*b));     
i0    = ceil(m/(4*log(n/m))); 
theta = mu*m/n/10;

end


%------------------Modifeid weights----------------------------------------
function w = ModWeight(x,h,theta,k,eps2)
    n         = length(x);
    w         = ones(n,1);
    eps1      = 1e-10; 
    [~,Ind] = sort(h,'descend');       
    if k==0
        w=1./(abs(x)+eps2);
    else
        w(Ind(1:k))  = eps1+theta*sum(h(Ind(2:k+1)))/sum(h(Ind(1:k)));
        w(Ind(k+1:n))= eps1+theta+1./(abs(x(Ind(k+1:n)))+eps2);
   end  
end

%------------------Update the Sparsity-------------------------------------
function sp = sparsity(x,rate)
    rs=rate*sum(x); y=0; sp=0;   
    while y < rs
    sp=sp+1; y=y+x(sp,1); 
    end    
end

