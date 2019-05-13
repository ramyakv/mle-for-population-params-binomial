clear all; close all; clc

Nlist = [1e6];
len_n = length(Nlist);

tlist = [2, 3 , 4, 5, 6, 7, 8, 9, 10, 11, 12];
len_t = length(tlist);

maxExpers = 10;

for n_ind = 1 : len_n
    
    % Population size
    N = Nlist(n_ind);

    for t_ind = 1 : len_t
        % Observations per individual
        t = tlist(t_ind);

        % grid for recovery
        m = 1000;

        for exper = 1 : maxExpers
            N
            t
            exper
            % generate distribution over p's

            %% Uniform on [0, 1]
            pdist = rand(1, N);

            %% Truncated Gaussian:
            %mu = 0.5;
            %sig = 0.1;
            %llim = zeros(1, N);
            %ulim = ones(1, N);
            %pdist = mu + sig .* truncatedNormal( (llim - mu)./sig , (ulim - mu)./sig )';

            %% Single spike at 1/2
            %pdist = 0.5.*ones(1, N);

            %% 3-spike
            %points = [1/4, 2/4, 3/4];
            %masses = [1/3, 1/3, 1/3];
            %pdist = pointMassMixture(points, masses, N);

            % double check to ensure pdist is on [0, 1]
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

            %% Observed Fingerprint
            Hobs = -ones(1, t+1);
            for i = 1 : t+1
                numheadseen = find(headseen == i-1);
                Hobs(i) = length(numheadseen)./N;
                clear numheadseen
            end

            %% MLE
            [phatmle, hhatmle, gridvecmle, optvalmle, statusmle] = learnpopulationparamsBinomialMLE( (Hobs(1:t+1))', t, m);

            %% Estimated moments from observations
            betaest = zeros(1, t);
            iter = 1;
            for k = 1 : t
                num = 0;
                for i = k : t
                    %pheads = length(find(headseen == i));
                    pheads = Hobs(i+1)*N;
                    %num = num + pheads .* nchoosek(i, k);
                    num = num + pheads .* exp( gammaln(i+1) - gammaln(k+1) - gammaln(i - k + 1) );
                end
                %den = nchoosek(t, k);
                den = exp( gammaln(t+1) - gammaln(k+1) - gammaln(t - k + 1) );
                betaest(iter) = (1/N) .* (num/den); 
                clear num
                clear den
            iter = iter + 1;
            end

            %% Estimate using moment matching
            w = ones(1, t);
            [phatM, bhat, gridvec, optval, status] = learnpopulationparamsBinomialwithMoments(betaest, t, w, m);

            %% W1 dist measured using l_1 of CDF

            % CDF of pdist
            hgrid = gridvecmle - 1/(2.*m);
            hgrid = [hgrid, 1 + 1/(2.*m)];
            hpdist = histcounts(pdist, hgrid, 'Normalization', 'probability');
            cdfhpdist = cumsum(hpdist);

            % Empirical cdf
            pempest = headseen./t;
            hemp = histcounts(pempest, hgrid, 'Normalization', 'probability');
            cdfemp = cumsum(hemp);

            % CDF of mle output
            cdfmle = cumsum(phatmle);
            cdfM = cumsum(phatM);

            % EMD for mle
            w1l1mle = norm(cdfmle - cdfhpdist', 1); 
            emd{n_ind}(exper, t_ind) = w1l1mle/m;

            % EMD for moment matching 
            w1l1M = norm(cdfM - cdfhpdist', 1); 
            emdM{n_ind}(exper, t_ind) = w1l1M/m;

            % EMD for empirical
            w1l1emp = norm(cdfemp - cdfhpdist, 1); 
            emdemp{n_ind}(exper, t_ind) = w1l1emp/m;
        
        end

    end
end
