function [fr,z,nf]=analyseFFT(y,Fe)
%analyse temporelle

n=length(y);
T = n/Fe; 
% t : temps discr�tis� pour la dur�e du signal (en n points)
t = (1:n)/n*T;
% trac� du signal temporel
% figure(1);
%subplot(2,2,1),plot(t,y);
%xlabel('Temps'); ylabel('Pression'); 

%analyse fr�quentielle

% Fmax =1024;
z = fft(y,size(y,1)*5);
T = (size(z,1)-1)/Fe ;
% fr = 0 :1/T :Fmax ;
fr = 0 :1/T :Fe ;
nf = length(fr) ;
%subplot(2,2,2),
% plot(fr,abs(z(1 :nf))) ;
% xlabel('fr�quence'); ylabel('Amplitude');

%analyse sonogramme

%subplot(2,2,3),specgram(y,1024,Fe) ;
%soundsc(y,Fe);
end


