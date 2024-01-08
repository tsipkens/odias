
% GEN_SMPS  Evaluates the transfer function of a differential mobility analyzer.
% 
%  INPUTS:
%   d_star              Particle diameter, measurement set point for DMA [m]
%   d                   Particle diameter, points in integral, can be vector [m]
% 
%   varargin{1} = 
%       prop	DMA properties, struct, generated using prop_DMA function (optional)
%   
%   varargin{2} = 
%       opts.diffusion  Indicates whether to include diffusion, boolean (optional)
%           .solver     Indicates the method by which diffusion is calculated (optional)
%           .param      String indicated which parameter set to use (see prop_DMA.m)
%
%  OUTPUTS:
%   Omega           Transfer function
%  
%  AUTHOR:	Timothy Sipkens, 2020-03-09

function [Omega, f_z, Omega_z, qbar] = gen_smps(d_star, d, z, argin_dma, argin_z)


%-- Parse inputs -----------------%
if ~exist('z', 'var'); z = []; end
if isempty(z); z = (1:4)'; end

if ~exist('argin_dma', 'var'); argin_dma = {}; end
if ~iscell(argin_dma); argin_dma = {argin_dma}; end

if ~exist('argin_z', 'var'); argin_z = {}; end
if ~iscell(argin_z); argin_z = {argin_z}; end
%---------------------------------%


%== Evaluate particle charging fractions =================================%
[f_z, qbar] = charge(d .* 1e-9, z, [], argin_z{:}); % get fraction charged for d


%== Evaluate DMA transfer function =======================================%
Omega_z = tfer_dma( ...  % evalute transfer function
    d_star' .* 1e-9, d .* 1e-9, ...
    z, argin_dma{:});

Omega_z = Omega_z .* permute(f_z, [3, 2, 1]);  % incorporate charge fraction

Omega = sum(Omega_z, 3);  % sum over multiple charge states

if nargout > 2
    Omega_z = permute(Omega_z, [1, 3, 2]);
end

disp(' ');


end



