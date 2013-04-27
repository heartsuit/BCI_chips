function [h_opt, g_opt, critere_min] = optim_one_wavelet(x_learn,y_learn,chans,ind_des,measure)


%measure = 'pce' ou 'fisher'

pas = 0.2; % discretisation step
a = -pi+pas:pas:pi;
lh_filtre = 4;
critere_min=inf;

nb_samples = size(x_learn,2);
nb_trials_all = length(y_learn);
i_opt = 1;

for ii=1:length(a)
    hh = construction_filtreh(lh_filtre,a(ii));
    gh = hh(end:-1:1);
    for jj=1:length(gh)
        gh(jj) = (-1)^jj * gh(jj);
    end
    nbz = nbturn(hh); % turncount of mother wavelet (hh)
    if nbz < 30
        % Computing of the marginales and the criterion
        N_dec = floor(log2(nb_samples));
        for i=1:nb_trials_all
            for k=chans
                [dwt_coeff(i,:,k), l(i,:,k)] = wavedec_simplif(x_learn(i,:,k),'per',N_dec,hh,gh,0,0);
                marginals(i,:,k) = calc_marg(dwt_coeff(i,:,k),l(i,:,k),ind_des);
            end
        end
        features = reshape(marginals,nb_trials_all,[]);
        if strcmp(measure,'pce')
            %Calcul du crière du pce
            [y_estim, y_real] = xval_procedure(features,y_learn,10);
            c_q = sum(y_estim ~= y_real)/length(y_estim);
        elseif strcmp(measure,'fisher')
            %Clacul du critère de fisher
            c_q = fct_calc_fisher(features,y_learn);
        end
        if (c_q < critere_min)
            critere_min=c_q;
            i_opt = ii;
        end
    end
end

% "optimal" filter
if exist('i_opt','var')
    h_opt = construction_filtreh(lh_filtre,a(i_opt));
    g_opt = h_opt(end:-1:1);
    for jj=1:length(g_opt)
        g_opt(jj) = (-1)^jj * g_opt(jj);
    end
end
    