%% Election Data
% electionDataMatrix
% headseen = electionData(:, 2);
% N = 3116; % size(electionData, 1);
% t = 8; 

% %% Flights Data
flightsDataMatrix
headseen = flightsData(:, 2);
N = size(flightsData, 1);
t = 10;

Hobs = -ones(1, t+1);
for i = 1 : t+1
    numheadseen = find(headseen == i-1);
    Hobs(i) = length(numheadseen)./N;
    clear numheadseen
end

m = 1000;
%% MLE
tic
[phatmle, hhatmle, gridvecmle, optvalmle, statusmle] = learnpopulationparamsBinomialMLE( (Hobs(1:t+1))', t, m);
toc()
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
hgrid = gridvecmle - 1/(2.*m);
hgrid = [hgrid, 1 + 1/(2.*m)];

% Empirical cdf
pempest = headseen./t;
hemp = histcounts(pempest, hgrid, 'Normalization', 'probability');
cdfemp = cumsum(hemp);

% CDF of mle output
cdfmle = cumsum(phatmle);
cdfM = cumsum(phatM);

figure
plot([0:m]./m, cdfmle, 'b-s', 'LineWidth', 2)
hold on
plot([0:m]./m, cdfM, 'g-.x','LineWidth', 2)
plot([0:m]./m, cdfemp, 'r--d','LineWidth', 2)
legend('MLE', 'TVK17', 'Empirical')

% fingerprints from empirical and moment estimates
gridvec = gridvecmle;
for kk = 1 : t+1
    hhatemp(kk) = exp( gammaln(t+1) - gammaln(kk) - gammaln(t - (kk - 1) + 1) ) * ((gridvec.^(kk-1) .* (1 - gridvec).^(t-kk+1)) * hemp');
    hhatM(kk) = exp( gammaln(t+1) - gammaln(kk) - gammaln(t - (kk - 1) + 1) ) * ((gridvec.^(kk-1) .* (1 - gridvec).^(t-kk+1)) * phatM);
end

% Moments from empirical and moment estimates
betaMLE = -ones(1, t);
betaEmp = -ones(1, t);
for i = 1 : t
betaMLE(i) = gridvecmle.^i * phatmle;
betaEmp(i) = gridvecmle.^i * hemp';
end

figure
plot([1:t], betaMLE, 'b-s', 'Linewidth', 2)
hold on
plot([1:t], bhat, 'g-.x', 'Linewidth', 2)
plot([1:t], betaEmp, 'r--d', 'Linewidth', 2)
plot([1:t], betaest, 'k-x', 'Linewidth', 2)
legend('MLE', 'TVK17', 'Empirical', 'Obs')