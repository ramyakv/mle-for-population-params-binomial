clear all; close all; clc

% Population size
N = 1e6;

%tlist = [10, 50, 100, 500, 1000];
tlist = [2, 4, 6, 8, 10, 50, 100, 250, 500, 1000];
len_t = length(tlist);

maxExpers = 10;

for t_ind = 1 : len_t
    % Observations per individual
    t = tlist(t_ind);

    % grid for recovery
    %m = 100;
    m = 1000;
    
    for exper = 1 : maxExpers
	t
	exper
        % generate distribution over p's
        %% Uniform
        % pdist = rand(1, N);

        %% Truncated Normal
        % mu = 0.5;
        % sig = 0.1;
        % llim = zeros(1, N);
        % ulim = ones(1, N);
        % pdist = mu + sig .* truncatedNormal( (llim - mu)./sig , (ulim - mu)./sig )';

        %% Single spike
        pdist = 0.5.*ones(1, N);

        %% 3-spike
        % points = [1/4, 1/2, 3/4];
        % masses = [1/3, 1/3, 1/3];
        % pdist = pointMassMixture(points, masses, N);

        % double check to ensure range [0, 1]
        pdist(pdist < 0) = 0;
        pdist(pdist > 1) = 1;

        %% observations
        if (t == 1)
            headseen = rand(t, N) < pdist;
        else
        %     observations = rand(t, N) < kron(ones(t, 1), pdist);
        %     headseen = sum(observations);
            headseen = 0;
            for tosses = 1 : t
                observations = rand(1, N) < pdist;
                headseen = headseen + observations;
            end
        end
        clear observations

        Hobs = -ones(1, t+1);
        for i = 1 : t+1
            numheadseen = find(headseen == i-1);
            Hobs(i) = length(numheadseen)./N;
            clear numheadseen
        end

        tp = t;

        %[phat, hhat, gridvec, optval, status] = learnpopulationparamsBinomialwithCountFractions( (Hobs(1:tp+1))', tp, m);

        [phatmle, hhatmle, gridvecmle, optvalmle, statusmle] = learnpopulationparamsBinomialMLE( (Hobs(1:tp+1))', tp, m);

        hgrid = gridvecmle - 1/(2.*m);
        hgrid = [hgrid, 1 + 1/(2.*m)];
        hpdist = histcounts(pdist, hgrid, 'Normalization', 'probability');
        cdfhpdist = cumsum(hpdist);

        pempest = headseen./t;
        hemp = histcounts(pempest, hgrid, 'Normalization', 'probability');
        cdfemp = cumsum(hemp);

        cdfmle = cumsum(phatmle);

        w1l1mle = norm(cdfmle - cdfhpdist', 1); 
        emd(exper, t_ind) = w1l1mle/m;

        w1l1emp = norm(cdfemp - cdfhpdist, 1); 
        emdemp(exper, t_ind) = w1l1emp/m;
    
    end

end
