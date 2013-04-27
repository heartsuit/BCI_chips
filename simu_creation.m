function [y_learn,signal_learn,f,f_class,StN,StNN]= simu_creation(nb_trials,J,time,nb_class,Fe,sigma,posi,perfect)
%Cr�ation de signaux simul�s avec un pic fr�quentiel diff�rent pour chaque
%classe. Un bruit blanc a �t� ajout� au tout au final. Le but est de
%v�rifier que le crit�re de s�lection de la mielleure base fonctionne en
%s'assurant que le crit�re fait bien choisir des bases suffisamment fine
%pour �viter d'avoir deux classes d�crite par une m�me valeur d'un
%descripteur...
%
%il n'y a qu'un chan dans notre simu
%
%nb_trials: nombre d'essai qu'on veut avoir au final
%J:profondeur de l'arbre
%nb_samples:nombre de points de l'�chantillon
%nb_class:nombre de classe � discerner
%Fe: fr�quence d'�chantillonnage du signal
%sigma: influence du bruit
%posi: endroit d'o� l'on veut faire d�marrer les classes (compris entre 0
%et 2^J-le nombre de classes)
dt=1/Fe;
t=0:dt:time-1/Fe;
%envelop
gauss=amgauss(length(t),length(t)/2,length(t)*3/8);
signal_learn(:,:,1) = zeros(nb_trials,length(t));
y_learn=zeros(nb_trials,1);
f=zeros(nb_trials,1);
f_class=zeros(nb_class,1);
StN=zeros(nb_trials,1);
StNN=zeros(nb_trials,1);
l=1;

for i=1:nb_class
    f_class(i)=(posi+i+0.5)*Fe/2^(J+1);
    
    
    if i==nb_class
        
        y_learn((i-1)*floor(nb_trials/nb_class)+1:end)=i;
    
        for j=1:length(y_learn((i-1)*floor(nb_trials/nb_class)+1:end))  
            %boucle sp�ciale pour la fin du vecteur pour ne pas n�gliger le 
            %reste dans la division du nombre d'essai par le nombre de classes

            if perfect
               sig=MakeWaveletPacket(J,i+posi,j,'Daubechies',4,length(t)); 
               f((i-1)*floor(nb_trials/nb_class)+j)=f_class(i); %trololol
            else

                g=(posi+(i+0.5)+(0+0.25*randn(1,1)))*Fe/2^(J+1); %fr�quence variant aurour de la fr�quence caract�risant la classe 
                
                f((i-1)*floor(nb_trials/nb_class)+j)=g;

                sig=gauss'.*sin(2*pi*g*t);
            end
            
            noise=sigma*randn(1,length(t));

            signal_learn((i-1)*floor(nb_trials/nb_class)+j,:,1) =...
                sig+noise;

                                  
        end
    else
        y_learn((i-1)*floor(nb_trials/nb_class)+1:i*floor(nb_trials/nb_class))=i;
        for j=1:floor(nb_trials/nb_class)

            if perfect
               sig=MakeWaveletPacket(J,i+posi,j,'Daubechies',4,length(t));
               f((i-1)*floor(nb_trials/nb_class)+j)=(posi+i+0.5)*Fe/2^(J+1);
            else

%                 g=(posi+(i+0.5)+(0+0.25*randn(1,1)))*Fe/2^(J+1); %fr�quence variant aurour de la fr�quence caract�risant la classe 
                
                g=(posi+i+0.8*rand(1,1))*Fe/2^(J+1);                


                f((i-1)*floor(nb_trials/nb_class)+j)=g;

                sig=gauss'.*sin(2*pi*g*t);
            end
            
            noise=sigma*randn(1,length(t));

            signal_learn((i-1)*floor(nb_trials/nb_class)+j,:,1) =...
                sig+noise;
        end
    end
    StN(l)=20*log10(std(sig)/std(noise));
    l=l+1; 
end


%%%%%%%%%%%%%%%%% PLOTTING SOME SIMULATED SIGNALS %%%%%%%%%%%%%%%%%%%%%%%%%
figure;
plot(signal_learn(10,:,1));
hold on;
% plot(signal_learn(19,:,1),'color','r');
% keyboard;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end