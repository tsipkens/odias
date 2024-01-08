
% IMPORT_DMA_ALT  Import a raw CPC file containing SMPS data. 
%  
%  AUTHOR: Timothy Sipkens, 2023-09-22
%=========================================================================%

function [data, d_star, prop_dma] = import_dma_alt(fn)

in = readcell(fn, 'Delimiter', '\t');

idx_data1 = find(strcmp(in(:, 1), 'Diameter Midpoint (nm)')) + 1;
idx_data2 = find(strcmp(in(:, 1), 'Scan Time (s)')) - 1;

%== Get data and d_star ========================%
data = in(idx_data1:idx_data2, :);
data = replace(data, ',', '.');
data = str2double(data);

d_star = data(:, 1);
data = data(:, 2:end);


%== Consult header =============================%
in0 = readcell(fn);

idx_data0 = find(strcmp(in0(:, 1), 'Diameter Midpoint (nm)'));
in0 = in; % (1:idx_data0, :);

opts = struct();
opts.params = 'custom';

idx = find(contains(in0(:,1), 'Inner Radius'));
opts.prop.R1 = str2double(replace(in0{idx, 2}, ',', '.')) ./ 100;

idx = find(contains(in0(:,1), 'Outer Radius'));
opts.prop.R2 = str2double(replace(in0{idx, 2}, ',', '.')) ./ 100;

idx = find(contains(in0(:,1), 'Characteristic Length'));
opts.prop.L = str2double(replace(in0{idx, 2}, ',', '.')) ./ 100;

idx = find(contains(in0(:,1), 'Sample Flow'));
Q1 = str2double(replace(in0{idx, 2}, ',', '.'));
Q2 = str2double(replace(in0{idx, 4}, ',', '.'));
opts.prop.Q_s = Q1/60/1000;  % sample flow [m^3/s]
opts.prop.Q_a = Q1/60/1000;  % exhaust flow [m^3/s]
opts.prop.Q_m = Q2/60/1000;  % sheath flow [m^3/s]
opts.prop.Q_c = Q2/60/1000;  % aerosol flow [m^3/s]

opts.prop.T = 298;
opts.prop.p = 1;

prop_dma = kernel.prop_dma(opts);

end
