
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
% 
%  ------------------------------------------------------------------------
% 
%  AUTHOR: Timothy Sipkens, 2021-10-06


function [Lambda, prop] = gen_pma(sp, m, d, z, prop, opt, varargin)

addpath tfer_pma; % add mat-tfer-pma package to MATLAB path


%-- Parse inputs ---------------------------------------------------------%
if ~exist('opt','var'); opt = []; end

% By default, use Taylor series solution baout rc (Case 1C) without diffusion.
% See Sipkens et al., Aerosol Sci. Technol. (2019) for more information.
if isempty(opt); opt = '1C'; end

% If not given, import default properties of PMA, 
% as selected by prop_pma function.
if ~exist('prop','var'); prop = []; end
if isempty(prop); prop = kernel.prop_pma; end

if ~exist('z','var'); z = []; end
if isempty(z); z = (1:4)'; end
%-------------------------------------------------------------------------%


% Compute charge state.
f_z = kernel.tfer_charge(d .* 1e-9, z, [], varargin{:}); % get fraction charged for d

% Assign tfer_pma function to use.
fun = str2func(['tfer_',opt]); % call relevant function from submodule

% For first charge state.
Lambda = f_z(1,:) .* ...
    fun(sp, m' .* 1e-18, d' .* 1e-9, z(1), prop); % CPMA transfer function

% Add additional charge states.
for ii=2:length(z)
    Lambda = Lambda + f_z(ii,:) .* ...
        fun(sp, m' .* 1e-18, d' .* 1e-9, z(ii), prop);
end

end
