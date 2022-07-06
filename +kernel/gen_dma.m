
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

function [Omega] = gen_dma(d_star, d, z, argin_dma, varargin)


%-- Parse inputs -----------------%
if ~exist('z', 'var'); z = []; end
if isempty(z); z = (1:4)'; end

if ~exist('argin_dma', 'var') argin_dma = {}; end
if ~iscell(argin_dma); argin_dma = {argin_dma}; end

% the rest of the arguments are passed to kernel.tfer_charge(...)
%---------------------------------%


n_b = length(d_star);
n_i = length(d);

%== Evaluate particle charging fractions =================================%
f_z = sparse(kernel.tfer_charge(d .* 1e-9, z, [], varargin{:})); % get fraction charged for d
n_z = length(z);


%== Evaluate DMA transfer function =======================================%
disp('Computing DMA kernel:');
Omega = sparse(n_b,n_i);
tools.textbar([0, n_b]);
for kk=1:n_z  % loop through charge states
    Omega_z = kernel.tfer_dma( ...
        d_star' .* 1e-9, ...
        d .* 1e-9, ...
        z(kk), argin_dma{:});
    
    % Remove numerical noise in kernel.
    Omega_z(Omega_z<(1e-7.*max(max(Omega_z)))) = 0;
    
    % Compile with other charge states.
    Omega = Omega + f_z(kk,:) .* sparse(Omega_z);
    
    tools.textbar([kk, n_z]);
end

disp(' ');


end



