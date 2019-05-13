% Algorithm to estimate the population parameters for binomial
function [q, b, gridvec, optval, status] = learnpopulationparamsBinomialwithMoments(beta, s, w, m)
% N is the number of individuals and t is the observations per individual
% Input: 
% beta: true moments
% s: Number of momets to be used for estimation
% w: weights given to each moment, a 1-by-s 
% m: number of points in the grid
%
% Output: 
% b: fit for moments found by linear programming
% q: probability mass function on the gridvector that fits the observed
% moments

% grid vector for estimating the pmf
gridvec = (1/m).* [0:m];

% linear program to estimate probability distribution to best match the
% moments


cvx_begin

variable q(m+1) nonnegative
variable b(s)

%minimize( w * abs(b - beta') )
%minimize( w * (b - beta').^2 )
minimize(  square_pos(norm( w' .* (b - beta'), 2)) )

subject to

for kk = 1 : s
    b(kk) == gridvec.^kk * q;
end
sum(q) == 1
q >= 0

cvx_end

optval = cvx_optval;
status = cvx_status;
