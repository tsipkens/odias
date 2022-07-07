
% GEN_DMA  Evaluates the transfer function of a differential mobility analyzer.
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

function [Omega, f_z, Omega_z, qbar] = gen_dma(d_star, d, z, argin_dma, argin_z)


%-- Parse inputs -----------------%
if ~exist('z', 'var'); z = []; end
if isempty(z); z = (1:4)'; end

if ~exist('argin_dma', 'var') argin_dma = {}; end
if ~iscell(argin_dma); argin_dma = {argin_dma}; end

if ~exist('argin_z', 'var') argin_z = {}; end
if ~iscell(argin_z); argin_dma = {argin_z}; end
%---------------------------------%


n_b = length(d_star);
n_i = length(d);
n_z = length(z);

%== Evaluate particle charging fractions =================================%
[f_z, qbar] = kernel.tfer_charge(d .* 1e-9, z, [], argin_z{:}); % get fraction charged for d


%== Evaluate DMA transfer function =======================================%
disp('Computing DMA kernel:');
Omega = sparse(n_b, n_i);

if nargout > 2
    Omega_z = zeros([size(Omega), length(z)]);
end

tools.textbar([0, n_z]);
for ii=1:n_z  % loop through charge states
    Omega_ii = f_z(ii,:) .* kernel.tfer_dma( ...
        d_star' .* 1e-9, d .* 1e-9, ...
        z(ii), argin_dma{:});
    
    % Remove numerical noise in kernel.
    Omega_ii(Omega_ii<(1e-7.*max(max(Omega_ii)))) = 0;
    
    % Compile with other charge states.
    Omega = Omega + Omega_ii;
    
    if nargout > 2
        Omega_z(:, :, ii) = Omega_ii;
    end
    
    tools.textbar([ii, n_z]);
end
if nargout > 2
    Omega_z = permute(Omega_z, [1, 3, 2]);
end

disp(' ');


end



