
% EXP_DIST_LPR  A helper function for 'exp_dist' to compute prior covariance.
%  
%  LPR = invert.exp_dist_lpr(LD, VEC) computes the prior covariance using a
%  length scale of LD and a vector of sizes on which the distribution is to
%  be reconstructed, VEC. 
%  
%  LPR = invert.exp_dist_lpr(..., BC) adds an input to specify boundary
%  conditions. Currently, BC = 0 pushes for 0s that boundary.
%  
%  ------------------------------------------------------------------------
%  
%  Author: Timothy Sipkens, 2019-12-11

function [Lpr, D, Gpr] = exp_dist_lpr(ld, vec, bc)

%-- Parse inputs ------%
if ~exist('bc', 'var'); bc = []; end
if isempty(bc); bc = 1; end


%-- Compute distances between elements -----------------------------------%
[vec_a, vec_b] = ndgrid(vec,vec);  % grid of distances

dr = log10(vec_a) - log10(vec_b);
D = abs(dr ./ ld); % distance, abs replaces sqrt(d.^2)


%-- Compute prior covariance matrix --------------------------------------%
Gpr = exp(-D);

Gpr_inv = pinv(Gpr);
[Lpr, ~] = chol(Gpr_inv);
clear Gpr_inv; % to save memory

% NOTE: Will cause incompatibility between Lpr and Gpr.
if bc == 0  % then push for 0 at boundary
    Lpr(1, :) = 0;
    Lpr(1, 1) = Lpr(2, 2);
    Lpr(end, :) = 0;
    Lpr(end, end) = Lpr(2, 2);
end

% Remove entries where relative distance to current 
% pixel is greater than some constant. 
Lpr(D > 1.75) = 0;
Lpr = sparse(Lpr);

end

