function [x,p,xx,w] = scalewave(h,kk);
%  [x,p,xx,w] = scalewave(h,kk)   calculates samples of the scaling function
%  phi(t) = p by  kk  successive approximations from the
%  scaling coefficients  h.  Initial iteration is a constant.
%  phi_k(t) is plotted at each iteration.     csb 5/19/93
%
%  et les coefficients de la fonction ondelette w

if nargin==1, kk=11; end;        % Default number of iterations
h2= h*2/sum(h);                  % normalize  h(n)
K = length(h2)-1; S = 128;       % Sets sample density
p = [ones(1,3*S*K),0]/(3*K);     % Sets initial iteration
P = p(1:K*S);                    % Store for later plotting
%axis([0 K*S+2 -.5 1.4]);
hu = upsample(h2,S);                % upsample h(n) by S
for iter = 0:kk                  % Successive approx.
   p = dnsample(conv(hu,p));     % convolve and down-sample
%   plot(p); pause;               % plot each iteration
%   P = [P;p(1:K*S)];             % store each iter. for plotting
end
p = p(1:K*S);                    % only the supported part
L = length(p);
x = ([1:L])/(S);
%axis([0 3 -.5 1.4]);
%plot(x,p);                       % Final plot
%title('Fonctions échelle et ondelette');

%ondelette
h2 = h*2/sum(h); NN = length(h2);
LL = length(p);
KK = round((LL)/(NN-1));
h1u = upsam(h2(NN:-1:1).*cos(pi*[0:NN-1]),KK);
w  = dnsample(conv(h1u,p));
w  = w(1:LL);
xx = [0:LL-1]*(NN-1)/(LL-1);
%plot(xx,w,'r');
