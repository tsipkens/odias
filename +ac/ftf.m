
% FTF  Full transfer function approach, including size distribution. 
%  Not an average charge method, rather acting as a baseline.
%  
%  AUTHOR: Timothy Sipkens, 2022-05-27

function o = ftf(m_star, kernel, p)

tools.textheader('Running FAC');  % add header to console

zmax = 300;
zbar = zeros(1, length(m_star));    % to store average charge state
fz0 = zeros(zmax, length(m_star));  % to store charge distribution

for ii=1:length(m_star)
    
    % First compute the mobility diameter. 
    d0 = (m_star(ii) .* (1:zmax) .* 1e-18 ./ prop.m0) .^ (1 / prop.Dm);

    % OPTIONAL: Convert to charging equivalent diameter.
    if f_deq  % check input flag
        di = working.dm2deq(di);
    end
    
    fz = kernel.tfer_charge(d0 .* 1e-9, 1:400, [], model, opt);
    fz0(:, ii) = diag(fz,0)';
    
    zbar(ii) = sum((1:zmax) .* diag(fz,0)') ./ sum(diag(fz,0));
    
    % Update progress bar for overall m_star loop.
    disp(' ');
    disp('Overall progress:');
    tools.textbar([0, length(m_star)]);
    tools.textbar([ii, length(m_star)]);
    
    %{
    % Debugging figure.
    figure(gcf); clf;
    addpath cmap;
    cmap_sweep(length(d0), inferno(length(d0)));
    plot(d0, fz);
    hold on;
    plot(d0, fz0(:,ii), 'ko');
    hold off;
    set(gca, 'XScale', 'log');
    %}
end

mbar = zbar' .* m_star;

tools.textheader();  % mark as complete

end
