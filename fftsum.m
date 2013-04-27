function [fft,f_scl]=fftsum(sig,Fe)

    mean_=mean(sig,2);
    for i=1:size(sig,1)
        sig(i,:)=sig(i,:)-mean_(i);
    end
    
    [f_scl_long,ffts,lgth_1]=analyseFFT(sig',Fe);
    
    f_scl=f_scl_long(1:floor(lgth_1/2));
    
    fft=sum(abs(ffts(1:floor(lgth_1/2),:)),2);
    fft=fft/max(fft);
end