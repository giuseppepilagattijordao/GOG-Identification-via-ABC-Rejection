function y = tridiag(b,a,c, f )
%catch error
if isnan(b(1))||isnan(a(1))||isnan(c(1))||isnan(f(1))
            display('stop tridiag')
            %return
        end
%  Solve the  n x n  tridiagonal system for y:
%
%  [ a(1)  c(1)                                  ] [  y(1)  ]   [  f(1)  ]
%  [ b(2)  a(2)  c(2)                            ] [  y(2)  ]   [  f(2)  ]
%  [       b(3)  a(3)  c(3)                      ] [        ]   [        ]
%  [            ...   ...   ...                  ] [  ...   ] = [  ...   ]
%  [                    ...    ...    ...        ] [        ]   [        ]
%  [                        b(n-1) a(n-1) c(n-1) ] [ y(n-1) ]   [ f(n-1) ]
%  [                                 b(n)  a(n)  ] [  y(n)  ]   [  f(n)  ]
%
%  f must be a vector (row or column) of length n
%  a, b, c must be vectors of length n (note that b(1) and c(n) are not used)
% some additional information is at the end of the file
% a
% b
% c

% b=[0;b];
% c=[c;0];
n = length(f);
v = zeros(n,1);   
y = v;
w = a(1);
y(1) = f(1)/w;
for i=2:n
    v(i-1) = c(i-1)/w;
    w = a(i) - b(i)*v(i-1);
    y(i) = ( f(i) - b(i)*y(i-1) )/w;
%     if isnan(y(i))||isnan(w)||isnan(v(i-1))||y(i)==Inf||w==Inf||v(i-1)==Inf||y(i)==-Inf||w==-Inf||v(i-1)==-Inf
%        %display('stop')
%        y(i)=1e6;
%    end
end
for j=n-1:-1:1
   y(j) = y(j) - v(j)*y(j+1);
%    if isnan(y(j))
%        %display('stop')
%        y(j)=1e6;
%    end
end



%  This is an implementation of the Thomas algorithm.  It does not overwrite a, b, c, f but 
%  it does introduce a working n-vector (v).
%%%%%  Example
% n = 5; a = 4*ones(n,1); b = ones(n,1); c = 3*ones(n,1);
% f = rand(n,1);
% y = tridiag(a,b,c,f);
%%%%%  check solution
% A = diag(a,0) + diag(ones(n-1,1),-1) + diag(3*ones(n-1,1),1)
% A*y - f
%%%%% Conditions that will guarantee the matrix equation can be solved using this algorithm: 
%%%%%  1. matrix strictly diagonally dominant
%%%%%  2. matrix diagonally dominant, c_i not zero for all i, and abs(b_n) < abs(a_n)
%  It has been tested on MATLAB, version R2010b and version R2012a
%  version: 1.0
%  March 9, 2013
end