
% MEAN_SQ_ERR  Calculates the mean squared error of the non-zero entries of b.
%  AUTHOR: Timothy Sipkens, 2021-10-21

function SIGMA = mean_sq_err(A, x, b)

sqErr = (A * x - b) .^ 2; % squared errors, normalized by expected error in b
SIGMA = mean(sqErr(b~=0)); % average square error for cases where b~= 0

end
%=========================================================================%
