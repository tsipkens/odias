
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
%   Zp_tilde        Non-dimensional electrical mobility, vector
%  
%  AUTHOR:	Timothy Sipkens, 2020-03-09

function [Omega] = gen_dma(d_star, d, varargin)

n_b = length(d_star);
n_i = length(d);

%== Evaluate particle charging fractions =================================%
z_vec = (1:3)';
f_z = sparse(kernel.tfer_charge(d.*1e-9, z_vec)); % get fraction charged for d
n_z = length(z_vec);


%== Evaluate DMA transfer function =======================================%
disp('Computing DMA kernel:');
Omega = sparse(n_b,n_i);
tools.textbar([0, n_b]);
for kk=1:n_z
    Omega_z = kernel.tfer_dma(...
        d_star' .* 1e-9,...
        d .* 1e-9,...
        z_vec(kk),varargin{:});
    
    % Remove numerical noise in kernel.
    Omega_z(Omega_z<(1e-7.*max(max(Omega_z)))) = 0;
    
    % Compile with other charge states.
    Omega = Omega + f_z(kk,:) .* sparse(Omega_z);
    
    tools.textbar([kk, n_z]);
end

disp(' ');


end
