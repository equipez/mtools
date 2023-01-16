function [f, g, H] = chrosen(x)
%Chained Rosenbrock function.
%
%   See
%   [1] Powell (2006) The NEWUOA software for unconstrained optimization without derivatives
%   [2] Toint (1978) Some Numerical Results Using a Sparse Matrix Updating Formula in Unconstrained Optimization
%
%   L. N. Vicente, S. Gratton, Z. Zhang, 2019

n=length(x);

%alpha = 10;
alpha = 4;

f=0;
g=zeros(n,1);
H=zeros(n,n);

%for i=1:n-1 % Essentially the same function with permuted x
%  f = f + (x(i)-1)^2+alpha*(x(i)^2-x(i+1))^2;
%%
%  g(i)   = g(i) + 2*(x(i)-1)+alpha*2*(x(i)^2-x(i+1))*2*x(i);
%  g(i+1) = g(i+1) - alpha*2*(x(i)^2-x(i+1));
%%
%  H(i,i)    =  H(i,i)+2+alpha*2*2*(3*x(i)^2-x(i+1));
%  H(i,i+1)  =  H(i,i+1)-alpha*2*2*x(i);
%  H(i+1,i)  =  H(i+1,i) -alpha*2*2*x(i);
%  H(i+1,i+1)=  H(i+1,i+1)+alpha*2;
%end
for i=1:n-1
  f = f + (x(i+1)-1)^2+alpha*(x(i+1)^2-x(i))^2;
%
  g(i+1)   = g(i+1) + 2*(x(i+1)-1)+alpha*2*(x(i+1)^2-x(i))*2*x(i+1);
  g(i) = g(i) - alpha*2*(x(i+1)^2-x(i));
%
  H(i+1,i+1)    =  H(i+1,i+1)+2+alpha*2*2*(3*x(i+1)^2-x(i));
  H(i+1,i)  =  H(i+1,i)-alpha*2*2*x(i+1);
  H(i,i+1)  =  H(i,i+1) -alpha*2*2*x(i+1);
  H(i,i)=  H(i,i)+alpha*2;
end

return
