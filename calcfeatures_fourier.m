% calc_features
% param=0%uniforme
% param=1%dyalique
function features = calcfeatures_fourier(signal_learn,Fe, Nfft,tabf_uniforme,tabf_dyadique,param)
%calcul de la puissance par bande
signal_learn = signal_learn - mean(signal_learn);
S= (abs(fft(signal_learn,Nfft)/Nfft)).^2;
f=(0:Nfft-1)*Fe/Nfft;
if param ==0
    for i=1:size(tabf_uniforme,2)-1
    k1_uni = round(Nfft*tabf_uniforme(i)/Fe)+1;
    k2_uni = round(Nfft*tabf_uniforme(i+1)/Fe)+1;
    features(i)=sum(S(k1_uni:k2_uni-1));
    end
else 
    for i=1:size(tabf_dyadique,2)-1
    k1_dya = round(Nfft*tabf_dyadique(i)/Fe)+1;
    k2_dya = round(Nfft*tabf_dyadique(i+1)/Fe)+1;
    features(i)=sum(S(k1_dya:k2_dya-1));
    end
end 
    
end



