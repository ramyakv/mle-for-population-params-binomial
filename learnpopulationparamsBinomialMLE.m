% Algorithm to estimate the population parameters for binomial
function [q, hhat, gridvec, optval, status] = learnpopulationparamsBinomialMLE(h, t, m)
%function [q, hhat, gridvec, optval, status] = learnpopulationparamsBinomialMLE(h, t, m, gridvec)
% N is the number of individuals and t is the observations per individual
% Input: 
% h: count fractions.
% t: number of tosses.
% m: number of points in the grid
%
% Output: 
% hhat: fit found by the linear program
% q: probability mass function on the gridvector that fits the observed
% moments

% grid vector for estimating the pmf
gridvec = (1/m).* [0:m];

%gridvec(1) = 0.01;
%gridvec(m+1) = 0.99;

%gridvec = (1/m) .* [0.1:m];
M = length(gridvec);

cvx_begin

variable q(M) nonnegative
variable hhat(t+1)

minimize( - h' * log(hhat) )

subject to

for kk = 1 : t+1
    %hhat(kk) == nchoosek(t, kk-1) * ((gridvec.^(kk-1) .* (1 - gridvec).^(t-kk+1)) * q);
    hhat(kk) == exp( gammaln(t+1) - gammaln(kk) - gammaln(t - (kk - 1) + 1) ) * ((gridvec.^(kk-1) .* (1 - gridvec).^(t-kk+1)) * q);
end
sum(q) == 1
q >= 0

cvx_end

optval = cvx_optval;
status = cvx_status;


% nchoosek(i, kj)
% exp( gammaln(i+1) - gammaln(kj + 1) - gammaln(i - kj + 1) )