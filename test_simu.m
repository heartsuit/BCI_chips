close all;

Fe=1024;
nb_class=5;
time=2;
J=6;
nb_trials=120;
sigma=0.5;
posi=0;
[y_learn,signal_learn,f,f_cl,StN]= simu_creation(nb_trials,J,time,nb_class,Fe,sigma,posi);

plot(signal_learn(1,:,1));
keyboard;
hold on

plot(signal_learn(2,:,1),'color','g');
plot(signal_learn(3,:,1),'color','b');
plot(y_learn,'color','r');

a=0;
[fr,z,nf]=analyseFFT(signal_learn([1,26,55,79,115],:,1)',Fe);

figure;
plot(f);

figure;
plot(fr,abs(z(1 :nf,1)),'color','g');
hold on;
plot(fr,abs(z(1 :nf,2)),'color','b');
plot(fr,abs(z(1 :nf,3)),'color','r');
plot(fr,abs(z(1 :nf,4)),'color','c');
plot(fr,abs(z(1 :nf,5)),'color','y');

zi=sum(abs(z),2);
figure;
plot(fr,abs(zi),'color','m');


a=gcf;

 clk=clock;
   
    saveas(gcf-1,['rec_graph/distrib_bb_' num2str(clk(2)) '-' num2str(clk(3)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
    saveas(gcf,['rec_graph/correl_circle_bb_' num2str(clk(2)) '-' num2str(clk(3)) '_' num2str(clk(4)) 'h' num2str(clk(5)) ], 'fig');
