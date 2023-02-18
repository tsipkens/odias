

function meancharge0 = collkernel2charge(collkernel, ntvec)

% Infer npmax from size of data.
npmax = length(collkernel) - 1;

% Normalized collision kenrel.
collkernelratio = collkernel ./ collkernel(1);

% Generate necessary dummy variables.
dummy1 = ones(npmax, npmax+1);  % by default no contribution below
dummy2 = zeros(npmax, 1);
for np=1:npmax
    for nj = 0:np
        dummy1(np, nj+1) = diff_product(collkernelratio, nj, np);
    end
    dummy2(np) = product_term(collkernelratio, np);
end

% Main loop to compute charge.
meancharge0 = zeros(1, length(ntvec));
for nn=1:length(ntvec)
    nt = ntvec(nn);
    meancharge = 0;
    fq(1) = exp(-collkernel(1) * nt);
    sscheck = fq(1);
    
    for np = 1:npmax
	    fq(np + 1) = sum(exp(-collkernel((0:np)+1) .* nt) ./ dummy1(np, (0:np)+1));
        fq(np + 1) = fq(np + 1) * dummy2(np);
        
        if isnan(fq(np + 1)); fq(np + 1) = 0; end

        sscheck = sscheck + fq(np + 1);
        meancharge = meancharge + np * fq(np + 1);
    end

    meancharge0(nn) = meancharge;
end

end


function out = diff_product(x, nj, np)

fl = (0:np) ~= nj;
out = prod(x(fl) - x(nj + 1));

end


function out = product_term(x, np)

out = prod(x(2:np));

end