function [a,d] = dwt_simplif(x,mode,Lo_D,Hi_D,shifta,shiftd)
%DWT Single-level discrete 1-D wavelet transform.
%
% A.Maitrot - 10/10/2003
% DWT (sur 1 seul niveau) de la toolbox wavelet simplifiee 
% mode = 'sym' : symmetric-padding
%      = 'sp0' : zero padding
%      = 'sp1' : constante padding (1er et dernier échantillon)
%      = 'per' : periodized extension (sans rallongement en taille)
% pour gagner en temps de calcul
% [a,d] = dwt_simplif(x,mode,Lo_D,Hi_D,shifta,shiftd)
%
%   DWT performs a single-level 1-D wavelet decomposition
%   with respect to particular wavelet filters
%   (Lo_D and Hi_D) that you specify.
%
%   [CA,CD] = DWT_SIMPLIF(X,Lo_D,Hi_D) computes the wavelet decomposition
%   as above given these filters as input:
%   Lo_D is the decomposition low-pass filter.
%   Hi_D is the decomposition high-pass filter.
%   Lo_D and Hi_D must be the same length.
%

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 12-Mar-96.
%   Last Revision: 02-Aug-2000.
%   Copyright 1995-2002 The MathWorks, Inc.
%   $Revision: 1.15 $

% Check arguments.

% Default: Shift and Extension : mode Symmetrization impose ('sym')


% Compute sizes.
lf = length(Lo_D);
lx = length(x);

% Extend, Decompose &  Extract coefficients.
% disp('lenEXT lenKEPT')

if mode == 'per',
    lenEXT = lf/2; lenKEPT = 2*ceil(lx/2);
else
    lenEXT = lf-1 ; lenKEPT = lx+lf-1 ;
end

if mode == 'sym',
    % y = wextend('1D',dwtEXTM,x,lenEXT); % extension des 2 cotes
    I = [lenEXT:-1:1 , 1:lx , lx:-1:lx-lenEXT+1];
    if lx<lf
        K = (I<1);
        I(K) = 1-I(K);
        J = (I>lx);
        while any(J)
            I(J) = 2*lx+1-I(J);
            K = (I<1);
            I(K) = 1-I(K);
            J = (I>lx);
        end
    end
    % disp('wextend :')
    y  = x(I) ;
    % fin y = wextend
elseif mode == 'sp0',
    ext_V = ones(1,lenEXT);
    ext_L = kron(ext_V,x(1)); ext_R = kron(ext_V,x(length(x)));
    y = [ext_L x ext_R];
elseif mode == 'sp1',
    d_L = x(1)-x(2);
    d_R = x(length(x))-x(length(x)-1);
    ext_V0 = ones(1,lenEXT); 
    ext_VL = [lenEXT:-1:1];  
    ext_VR = [1:lenEXT];
    ext_L = kron(ext_V0,x(1)) + kron(ext_VL,d_L);
    ext_R = kron(ext_V0,x(length(x))) + kron(ext_VR,d_R);
    y = [ext_L x ext_R];
elseif mode == 'per',
    if rem(lx,2) , x(lx+1) = x(lx); lx = lx+1; end
    I = [lx-lenEXT+1:lx , 1:lx , 1:lenEXT];
    if lx<lenEXT
        I = mod(I,lx);
        I(I==0) = lx;
    end
    y = x(I);
   
end    

% a = convdown(y,Lo_D,lenKEPT,shift);
a = conv2(y(:)',Lo_D(:)') ; if size(x,1)>1, a = a' ; end 
% disp('conv :')
% a
% NB : conv2 est + rapide que conv

% a = wkeep(a,lenKEPT); % on garde la partie centrale
sx = length(a);
ok = (lenKEPT>=0) & (lenKEPT<sx) & (lenKEPT == fix(lenKEPT));
if ok==0 , first = 1; last = lenKEPT; return; end
dd = (sx-lenKEPT)/2;
first = 1+floor(dd); last = sx-ceil(dd);
if ok , a = a(first(1):last(1));end
% disp('wkeep : avant sous echantillonnage')
% a
% fin a = wkeep 

% disp('dyaddown')
a = a(2-rem(shifta,2):2:end) ;

% d = convdown(y,Hi_D,lenKEPT,shift);
d = conv2(y(:)',Hi_D(:)') ; if size(x,1)>1, y = y' ; end 

% d = wkeep(d,lenKEPT); % on garde la partie centrale
sx = length(d);
ok = (lenKEPT>=0) & (lenKEPT<sx) & (lenKEPT == fix(lenKEPT));
if ok==0 , first = 1; last = lenKEPT; return; end
dd = (sx-lenKEPT)/2;
first = 1+floor(dd); last = sx-ceil(dd);
if ok , d = d(first(1):last(1));end
% fin d = wkeep

d = d(2-rem(shiftd,2):2:end); 


% %-----------------------------------------------------%
% % Internal Function(s)
% %-----------------------------------------------------%
% function y = convdown(x,f,lenKEPT,shift)
% 
% % y = wconv('1D',x,f); % remplace par :
% y = conv2(x(:)',f(:)') ; if size(x,1)>1, y = y' ; end
%  
% % fin wconv
% 
% y = wkeep(y,lenKEPT); % on garde la partie centrale
% 
% % y = dyaddown(y,shift); % remplace par :
% y = y(2-rem(shift,2):2:end);
%  
% % fin dyaddown
% %-----------------------------------------------------%
