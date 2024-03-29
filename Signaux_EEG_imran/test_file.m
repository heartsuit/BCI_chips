% Représentation des moyennes signaux
close all
clc

t_start = 0.0;
t_end = 0.65;
N_start = round((2+t_start)*1024);
N_end = round((2+t_end)*1024);
N_tronc = N_start:N_end;
t_tronc = t_start:(t_end-t_start)/(length(N_tronc)-1):t_end;

chan = [1 2 3 4 5 6];

% couleur = ['c','r','g','k','b','m'];
couleur = ['b','b','b','b','b','b'];
for k = 5
    eval (['class = class_ep_' num2str(k)]);
    
    srate = Fe;
    N = 4;  %Ordre du filtre
    R = 20;    %Ripple du filtre
    W1 = (2*1)/srate;   %Fréquence de coupure basse (1Hz)
    W2 = (2*10)/srate;  %Fréquence de coupure haute (10Hz)
    
    Wp = [W1 W2];
    
    [B,A] = butter(N,Wp);
    %     [B,A] = butter(N,W1,'high');
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
    figure
    w=0;
    for i = 1:8
        figure(i)
        axes('FontSize',14);
    end
    for i = [1 2 3 6]
        class(1).signal_tronc = squeeze(class(1).signal(:,N_tronc,i));
        class(2).signal_tronc = squeeze(class(2).signal(:,N_tronc,i));
        
        w = w+1;
        figure(w)
%         
%         subplot(4,2,w)
        
        hold on
%         title(['\fontsize{16}' chan_label(i) ' sans ErrP'])
%         ylabel('\fontsize{16}sans ErrP')
        plot(t_tronc,mean(class(1).signal_tronc,1),'b--','LineWidth',3);
        plot(t_tronc,mean(class(1).signal_tronc,1)+std(class(1).signal_tronc,0,1),['b' '--'],'linewidth',1);
        plot(t_tronc,mean(class(1).signal_tronc,1)-std(class(1).signal_tronc,0,1),['b' '--'],'linewidth',1);
        
        
        plot(t_tronc,mean(class(2).signal_tronc,1),'r','LineWidth',3);
        plot(t_tronc,mean(class(2).signal_tronc,1)+std(class(1).signal_tronc,0,1),['r' '-'],'linewidth',1);
        plot(t_tronc,mean(class(2).signal_tronc,1)-std(class(1).signal_tronc,0,1),['r' '-'],'linewidth',1);
        axis([t_start t_end -0.8 0.8]);
        set(gca,'YTickLabel',[])
        hold off
        xlabel('\fontsize{16}temps en s')
        set(gcf,'position',[520 378 1000 420]);
        set(gcf,'PaperPositionMode','auto');
        set(gca,'YTick',[])
        print(gcf,'-depsc','-loose',['sub_' num2str(k) '_both_' chan_label{i}]);
%         w = w+1;
%         figure(w)
%         
%         subplot(4,2,w)
%         hold on
%         title(['\fontsize{16}' chan_label(i) ' avec ErrP' ])
%         ylabel('\fontsize{16}avec ErrP')
%         plot(t_tronc*1000,mean(class(2).signal_tronc,1),couleur(k),'LineWidth',2);
%         plot(t_tronc*1000,mean(class(2).signal_tronc,1)+std(class(1).signal_tronc,0,1),[couleur(k) '-']);
%         plot(t_tronc*1000,mean(class(2).signal_tronc,1)-std(class(1).signal_tronc,0,1),[couleur(k) '-']);
%         axis([t_start*1000 t_end*1000 -0.8 0.8]);
%         
%         hold off
%         xlabel('\fontsize{16}temps en ms')
%         set(gcf,'position',[520 378 1000 420]);
%         set(gcf,'PaperPositionMode','auto');
%         print(gcf,'-depsc','-loose',['sub_' num2str(k) '_faux_' chan_label{i}]);
    end
end