
% GEN_PMA  Bridging function used to evaluate particle mass analyer (PMA) transfer function.
% 
%  INPUTS:
%   sp          Structure defining various setpoint parameters 
%               (e.g. m_star, V). Use 'get_setpoint' method to generate 
%               this structure.
%   m           Particle mass
%   d           Particle mobility diameter
%   z           Integer charge state
%   prop        Device properties (e.g. classifier length) (optional)
%   opt         Alphanumeric code for transfer function evaluation method
%               (optional)
%
%  OUTPUTS:
%   Lambda      Transfer function
%   prop        CPMA device settings
%   f_z         Charge fractions
%   Lambda_z    Transfer function with individual charge contributions
%   qbar        Average charge on particles (not on transmitted particles)
% 
%  ------------------------------------------------------------------------
% 
%  AUTHOR: Timothy Sipkens, 2021-10-06


function [Lambda, prop, f_z, Lambda_z, qbar] = gen_pma(sp, m, d, z, prop, opt, varargin)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('opt','var'); opt = []; end

% By default, use Taylor series solution baout rc (Case 1C) without diffusion.
% See Sipkens et al., Aerosol Sci. Technol. (2019) for more information.
if isempty(opt); opt = '1C'; end

% If not given, import default properties of PMA, 
% as selected by prop_pma function.
if ~exist('prop','var'); prop = []; end
if isempty(prop); prop = prop_pma; end

if ~exist('z','var'); z = []; end
if isempty(z); z = (1:4)'; end
%-------------------------------------------------------------------------%


% Compute charge state.
[f_z, qbar] = charger(d, z, [], varargin{:}); % get fraction charged for d

Lambda_z = tfer_pma(sp, m, d, z, prop, opt);  % PMA transfer function

Lambda_z = Lambda_z .* permute(f_z, [3, 2, 1]);  % incorporate charge fraction

Lambda = sum(Lambda_z, 3);  % sum over multiple charge states

if nargout > 3
    Lambda_z = permute(Lambda_z, [1, 3, 2]);
end

end

