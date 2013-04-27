function [c,L] = wavedec_simplif(x,mode,n,Lo_D,Hi_D,shifta,shiftd)
% [c,L] = wavedec_simplif(x,mode,n,Lo_D,Hi_D,shifta,shiftd)
% A.Maitrot - 05/03/2004
% DWT (sur plusieurs niveaux) de la toolbox wavelet simplifiee 
% mode = 'sym' : symmetric-padding
%      = 'sp0' : zero padding
%      = 'sp1' : constante padding (1er et dernier échantillon)
%      = 'per' : periodized extension (sans rallongement en taille)
% pour gagner en temps de calcul
% 
% ATTENTION : compte-tenu de la gestion des effets de bord (ici, mode
% 'sym' : on rajoute 2*length(h)-2 au signal à chaque niveau), le niveau maximum de
% décomposition n'est pas log2(N) !!! 
% (N = length(x)) 
% Le  niveau max  est : fix(log2(lx/(lh-1))) (cf. wmaxlev.m)
%
% wavedec_simplif est une simplification de wavedec (toolbox wavelet de matlab) 
% dont le "help" est le suivant :
%
%   WAVEDEC Multi-level 1-D wavelet decomposition.
%   WAVEDEC performs a multilevel 1-D wavelet analysis
%   using either a specific wavelet 'wname' or a specific set 
%   of wavelet decomposition filters (see WFILTERS).
%   [C,L] = WAVEDEC(X,N,'wname') returns the wavelet
%   decomposition of the signal X at level N, using 'wname'.
%
%   N must be a strictly positive integer (see WMAXLEV).
%   The output decomposition structure contains the wavelet
%   decomposition vector C and the bookkeeping vector L.
%
%   For [C,L] = WAVEDEC(X,N,Lo_D,Hi_D),
%   Lo_D is the decomposition low-pass filter and
%   Hi_D is the decomposition high-pass filter.
%
%   The structure is organized as:
%   C      = [app. coef.(N)|det. coef.(N)|... |det. coef.(1)]
%   L(1)   = length of app. coef.(N)
%   L(i)   = length of det. coef.(N-i+2) for i = 2,...,N+1
%   L(N+2) = length(X).
%
%   See also DWT, WAVEINFO, WAVEREC, WFILTERS, WMAXLEV.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 04-Dec-2001.
%   Copyright 1995-2002 The MathWorks, Inc.
%   $Revision: 1.15 $ $Date: 2002/03/28 17:27:30 $

% Check arguments.



% Initialization.
s = size(x); x = x(:)'; % row vector
c = [];      L = [length(x)];

lev = fix(log2(max(s)/(length(Lo_D)-1)));
% if n > lev, disp('Niveau max de décomposition déraisonnable...')
%     disp(['il vaudrait mieux : ',num2str(lev)])
% end


for k = 1:n
    [x,d] = dwt_simplif(x,mode,Lo_D,Hi_D,shifta,shiftd); % decomposition
    c     = [d c];            % store detail
    L     = [length(d) L];    % store length
end

% Last approximation.
c = [x c];
L = [length(x) L];

if s(1)>1, c = c'; L = L'; end