close all
class = class_mvt;
N_start = round(4*1024);
N_end = round(6*1024)-1;
filtering = 0;
if filtering
    srate = Fe;
    N = 4;  %Ordre du filtre
    R = 20;    %Ripple du filtre
    W1 = (2*1)/srate;   %Fréquence de coupure basse (1Hz)
    W2 = (2*0.5)/srate;  %Fréquence de coupure haute (10Hz)
    
    Wp = [W1 W2];
    [B,A] = butter(N,W2);
%     [B,A] = cheby2(N,R,Wp); %Filtre de chebychev
%     [B,A] = cheby2(N,R,W2);
    for i=chan
        for j=1:size(class(1).signal,1)
            class(1).signal(j,:,i) = filter(B,A,class(1).signal(j,:,i));
        end
        for j=1:size(class(2).signal,1)
            class(2).signal(j,:,i) = filter(B,A,class(2).signal(j,:,i));
        end
    end
end

class(1).mean = squeeze(mean(class(1).signal,1));
class(2).mean = squeeze(mean(class(2).signal,1));
figure
hold on
plot(class(1).mean(:,6))
plot(class(2).mean(:,6),'r')
hold off
% figure
% cwt(class(1).mean(1:7000,6),500:512,'sym4','abslvl');
% figure
% cwt(class(2).mean(1:7000,6),500:512,'sym4','abslvl');

figure
N_freq=2^11;

for i=1:60
    [s1,f1] = pmusic(class(1).signal(i,N_start:N_end,6),50,N_freq,Fe);
    [s2,f2] = pmusic(class(2).signal(i,N_start:N_end,6),50,N_freq,Fe);
    hold on
    plot(f1,20*log(s1),'b')
    plot(f2,20*log(s2),'r')
    hold off
end

% f = ((1:N_freq)-1);
% figure
% for i=1:60
%     s1 = abs(fft(class(1).signal(i,N_start:N_end,6),N_freq));
%     s2 = abs(fft(class(2).signal(i,N_start:N_end,6),N_freq));
%     hold on
%     plot(f,20*log(s1),'b')
%     plot(f,20*log(s2),'r')
%     hold off
% end
% for i=1:size(class(1).signal,1)
%     figure
%     hold on
%     plot(class(1).signal(i,:,6))
%     plot(class(2).signal(i,:,6),'r')
%     hold off
%     pause;
% end