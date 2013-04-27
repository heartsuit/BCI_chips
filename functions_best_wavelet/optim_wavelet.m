function optim_one_wavelet(x_learn,y_learn,c_q)

pas = 0.2; % discretisation step
a = -pi+pas:pas:pi;
lh_filtre = 4;
critere_min=c_q;

for i=1:length(a)
   hh = construction_filtreh(lh_filtre,a(i));
   nbz = nbturn(hh); % turncount of mother wavelet (hh)
    if nbz < 30
        % Computing of the marginales and the criterion
        MxK = marginal_dwt_popu(xK,hh,lh_filtre,deb); % compute the marginals of channel K
        nbf = size(MxK,1);
        Rx((K-1)*nbf+1:nbf*K,:,:) = MxK; % update the global marginals
        critere = crossval_learn_SVM_OVRs (Rx, nbloc, svmparam);

        if (critere < critere_max)
            critere_min=critere;
            i_opt = i;
        end;
    end;
end
wav_opt_param=[-pi+pas*(i_opt+1)];

% "optimal" filter
[h_opt,par] = construction_filtreh(lh_filtre,[a(i_opt)]);
end
    