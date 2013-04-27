
figure;
hold on
for k=1:100
%     input(k,:)=MakeWaveletPacket(2,0,k,'Daubechies',4,2048);
%     [scal,fft]=fft(input');
    y=fft(MakeWaveletPacket(2,0,k,'Daubechies',4,2048));
    plot(abs(y));
end