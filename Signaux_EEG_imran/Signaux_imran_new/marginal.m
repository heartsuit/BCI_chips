% function y = marginal (x,deb)
% marginales of x, x x being the vector of a dyadic decomposition

function y = marginal (x,deb)

N=length(x);
NL=log2(N);

% square value in place of absolute value
% x=x/1000;
% x = x.^2;

x=x(end:-1:1);
% x=x./sum(x); % normalisation
x=x+1e-10; % for the log of the Kullback distance

for j=1:NL-deb
    c(j) = sum(x(1:N/2^j));
    x(1:N/2^j) = [];
end

%c(j+1) = sum(x); %the approximation
y=c;
