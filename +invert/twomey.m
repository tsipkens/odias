
% TWOMEY  Performs inversion using the iterative Twomey approach.
% 
%  X = invert.twomey(A,B) inverts the system A*X = B using Twomey 
%  inversion, where A is a model matrix or kernel and B is a data vector.
% 
%  X = invert.twomey(A,B,XI) applies Twomey inversion using an initial 
%  guess of XI. If excluded, a vector of ones is used. 
% 
%  X = invert.twomey(A,B,XI,ITER) applies Twomey inversion using ITER 
%  number of iterations. If excluded, 100 iterations is used. 
% 
%  X = invert.twomey(A,B,XI,ITER,F_SIG) adds a flag for whether or not 
%  to check for convergence of the algorithm (mean square error < 1).
% 
%  X = invert.twomey(...,F_BAR) applies Twomey inversion, 
%  as above, with F_BAR flagging whether a progress bar is 
%  output to the console. 
% 
%  AUTHOR:  Timothy Sipkens, 2018-11-21

function [x] = twomey(A, b, xi, iter, f_sigma, f_bar)

%-- Parse inputs ---------------------------------------------------------%
if ~exist('iter', 'var'); iter = []; end
if isempty(iter); iter = 100; end % default of 100 iterations

if ~exist('xi', 'var')
    xi = ones(size(A,2), 1); % intiate vector of ones if xi is not specified
end
x = xi;

% whether to display text-based progress bar
if ~exist('f_bar', 'var'); f_bar = []; end
if isempty(f_bar); f_bar = 0; end

% controls exit conditions
if ~exist('f_sigma', 'var'); f_sigma = []; end
if isempty(f_sigma); f_sigma = 1; end
%-------------------------------------------------------------------------%


%-- Start evaluation -----------------------------------------------------%
s = 1 ./ max(A, [], 2); % factor to scale data and model matrix
s_min = max(max(A, [], 2)) * 1e-10;
s(max(A, [], 2) < s_min) = 0;
A = bsxfun(@times, s, A); % scale model matrix
b = max(b .* s, 0); % remove any negative data, scale data

if f_bar; disp('Twomey progress:'); tools.textbar([0,iter]); end
    % output textbar for progress

factor = 1; % factor, allows one to decrease step size in Twomey


%-- Perform Twomey iterations --------------------------------------------%
for kk=1:iter % perform multiple Twomey passes
    for ii=1:length(b) % loop through data vector for single Twomey pass
        if b(ii)~=0 % do not consider null data
            y = A(ii,:) * x;
            if y~=0 % prevents errors in dividing by zero
                X = b(ii) ./ y; % calculate weight factor
                C = 1 + factor .* (X-1) .* A(ii,:)';
                x = C .* x;
            end
        end
    end
    
    if f_bar; tools.textbar([kk,iter]); end % outputs progresss

    % Exit conditions if sigma is specified (e.g. Twomey-Markowski)
    if f_sigma
        if tools.mean_sq_err(A, x, b)<0.6
            if f_bar
                tools.textbar([iter,iter]);
                disp([' Algorithm converged at ', num2str(kk), ' iterations']);
            end
            break;
        end
    end
end


end
