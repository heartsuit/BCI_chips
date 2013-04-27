%   CSP Function

%   Coded by James Ethridge and William Weaver

function [f,class] = CSP(class,chan)
    
    nb_classes = size(class,2);
    if nargin < 2
        chan = 1:size(class(1).signal,3);
    else
        for i=1:nb_classes
            class(i).signal = class(i).signal(:,:,chan);
        end
    end
    if (nb_classes ~= 2)
        disp('Must have 2 classes for CSP!')
    end
    
    
    nb_chan = length(chan);
    
    
%     Rmoy = zeros(nb_classes,nb_chan,nb_chan);
    %finding the covariance of each class and composite covariance
    Rsum = 0;
    for i = 1:nb_classes
        nb_trials = size(class(i).signal,1);
        class(i).R = zeros(nb_trials,nb_chan,nb_chan);
        for j = 1:nb_trials
            E = squeeze(class(i).signal(j,:,:));
            E = E';
            class(i).R(j,:,:) = (E*E');
            class(i).R(j,:,:) = squeeze(class(i).R(j,:,:))/trace(squeeze(class(i).R(j,:,:)));
        end
        class(i).Rmoy = squeeze(mean(class(i).R,1));
        %Ramoser equation (2)
        Rsum = Rsum+class(i).Rmoy;
    end
    
    
    %   Find Eigenvalues and Eigenvectors of RC
    %   Sort eigenvalues in descending order
    [EVecsum,EValsum] = eig(Rsum);
    [EValsum,ind] = sort(diag(EValsum),'descend');
    EVecsum = EVecsum(:,ind);
    
    %   Find Whitening Transformation Matrix - Ramoser Equation (3)
        W = sqrt(inv(diag(EValsum))) * EVecsum';
    
    S = zeros(nb_classes,nb_chan,nb_chan);
    for i = 1:nb_classes
        class(i).S = W * class(i).Rmoy * W'; %       Whiten Data Using Whiting Transform - Ramoser Equation (4)
    end
    
    
    
    % Ramoser equation (5)
   % [U{1},Psi{1}] = eig(S{1});
   % [U{2},Psi{2}] = eig(S{2});
    
    %generalized eigenvectors/values
    [B,D] = eig(class(2).S);
    % Simultanous diagonalization
			% Should be equivalent to [B,D]=eig(S{1});
    
    %verify algorithim
    %disp('test1:Psi{1}+Psi{2}=I')
    %Psi{1}+Psi{2}
    
    %sort ascending by default
    %[Psi{1},ind] = sort(diag(Psi{1})); U{1} = U{1}(:,ind);
    %[Psi{2},ind] = sort(diag(Psi{2})); U{2} = U{2}(:,ind);
    [max_D,ind]=max(diag(D));
    B=B(:,ind);
    %Resulting Projection Matrix-these are the spatial filter coefficients
    f = B'*W;
    for i=1:nb_classes
        nb_trials = size(class(i).signal,1);
        for j = 1:nb_trials
            class(i).signal_SpatFilter(j,:) = f*squeeze(class(i).signal(j,:,:))';
        end
    end
end
