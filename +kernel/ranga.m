

% RANGA A function to determine the average charge on particles.
%   ORIGINAL AUTHOR: Ranganathan Gopalakrishnan
%   MODIFIED BY: Timothy Sipkens, 2023-02-08

clear;

% physical constants
k = 1.381e-23;
e = 1.602e-19;
N_av = 6.023e23;
eps = 8.85e-12;  % J/K, C, #/mole
R_univ = 8.314e0;  % dimless, J/kg/mole
C1 = 25.836;
C2 = 11.211;
C3 = 3.502;
C4 = 7.211;

% gas properties
p = 1.01325e5; % Pa
T = 2.98e2; % K
mwt = 2.9e-2; % kg/mole
mu = 1.81e-5; % Pa.s
alpha = 1.257e0; beta = 4e-1; gamma=1.1e0;

% ion properties
MODEL = 'SIPKENS'
ni = 1;

% particle properties
PA_nd=3.14159e0; % dimless
Rs_nd=1e0; % dimless
dielec=13.5;  % dielectric constant
georatio=10;
STARTnt=1.0d12; ENDnt=1.0d13; TAG='06P76D12';

	
if strcmp(MODEL, 'MOHSEN')
    
    Mwt_ion = 0.1
    zp_ion  = 1.9d-4 % kg/mole, m^2/s/V

elseif strcmp(MODEL, 'VOHRAA')
    
    mwt_ion = 0.05
    zp_ion  = 1.9d-4 % kg/mole, m^2/s/V

elseif strcmp(MODEL, 'SIPKENS')
    
    mwt_ion = 0.109
    zp_ion  = 1.4d-4 % kg/mole, m^2/s/V
	
end
	
		
% Mobility equivalent diameter
d_me_array = logspace(log10(4), log10(2e3), 120);
d_me_array = d_me_array([1,40,65]);

meancharge_d =[];
for i = 1:length(d_me_array)
    
    d_me = d_me_array(i) * 1e-9;
    
    disp('Running:')
    disp(['d_me = ', num2str(d_me * 1e9)]);
    
    % calculation of gas properties
    R_gas = R_univ/mwt; % J/kg/K
    pho = p/(R_gas*T); % kg/m^3
    m_gas = mwt/N_av; % mass of gas molecule, kg
    mts = (8*k*T/pi/m_gas)^0.5; % m/s
    mfp = mu/(0.499*pho*mts); % hard sphere mean free path, m
    
    % calculation of primary particle radius of aggregate with the same mobility diameter
    R_me = d_me/2;
	    
    rp = R_me;
    PA = PA_nd*(rp^2); % m^2
    Rs = Rs_nd*(rp); % m
    
    % calculation of ion properties
    m_i = mwt_ion/N_av;
    f_i = ni*e/zp_ion;
    
    % Adjust npmax depending on the particle size. 
    % This could cause discontinuities depending on the conditions.
    if (d_me*1e9 <= 10)
	    npmax = 5;
    elseif (d_me*1e9 <= 50)
	    npmax = 10;
    elseif (d_me*1e9 <= 400)
	    npmax = 20;
    elseif (d_me*1e9 <= 800)
	    npmax = 50;
    elseif (d_me*1e9 <= 1200)
	    npmax = 100;
    else
	    npmax = 300;
    end
    
    psiEarray = [];
    psiIarray = [];
    etafarray = [];
    etacarray = [];
    KnDarray = [];
    collkernel = [];
    for np=0:npmax
	    psiE = -np*ni*(e^2)/(4*pi*eps*k*T*Rs);
	    psiI =  ((ni*e)^2)/(4*pi*dielec*eps*k*T*Rs);
	    
	    psiEarray(np+1) = psiE;
	    psiIarray(np+1) = psiI;
    
        if ((psiE == 0) && (psiI == 0))
		    eta_c = 1;
        else
		    SI = simpsons_cont(0,0.99999e0,1000,psiE,psiI);
		    eta_c = 1.0/SI;
        end
    
        if ((psiE == 0) && (psiI == 0))
		    eta_f = 1e0;

        elseif (psiI == 0)
		    eta_f = exp(psiE);
    
        else
            if (psiE < 0)
                eta_f = eta_fm_repulsive(psiE, psiI);
            else
                eta_f = eta_fm_attractive(psiE, psiI);
            end
    
        end
    
	    etafarray(np + 1) = eta_f;
	    etacarray(np + 1) = eta_c;
    
        L_KnD = PA*eta_f/(pi*Rs*eta_c);
        L_H3 = PA*PA*eta_f*eta_f/(pi*pi*Rs*eta_c);
        
        KnD=(m_i*k*T)^0.5/(f_i*L_KnD);
        H = (4*pi*(KnD^2) + C1*(KnD^3) + ((8*pi)^0.5)*C2*(KnD^4) )/(1 + C3*KnD + C4*(KnD^2) + C2*(KnD^3) );
        KnDarray(np + 1) = KnD;
        collkernel(np + 1) = H*L_H3*f_i/m_i;
	    
    end

    collkernelratio = collkernel ./ collkernel(1);
    
    
    nt = STARTnt;
    
    meancharge_nt =[];
    while (nt <= ENDnt)
        meancharge=0e0;
        ssfrac(1) = exp(-collkernel(1) * nt);
        sscheck = ssfrac(1);
	    
        for np=1:npmax
	        
	        ssfrac(np + 1) = 0;
        
            for nj = 0:np
	            
		        dummy1 = diff_product(collkernelratio, npmax, nj, np);
		        ssfrac(np + 1) = ssfrac(np + 1) + exp(-collkernel(nj + 1)*nt)/dummy1;
            
            end
	        dummy2 = product_term(collkernelratio, npmax, np);
	        ssfrac(np + 1) = ssfrac(np + 1)*dummy2;
	        % write(*,*) np, ssfrac(np + 1)
		        
	        sscheck = sscheck + ssfrac(np + 1);
	        meancharge = meancharge + np * ssfrac(np + 1);
        end
        
        nt = nt*georatio;

        meancharge_nt = [meancharge_nt, meancharge];
    
    end
    meancharge_d = [meancharge_d; meancharge_nt]

end



function out = diff_product(collkernelratio,~,nj,np)

fl = (0:np) ~= nj;
out = prod(collkernelratio(fl) - collkernelratio(nj + 1));

end


function out = product_term(collkernelratio, ~, np)

out = prod(collkernelratio(2:np));

end


function out = Cc(Kn,alpha,beta,gamma)

out = 1 + Kn*(alpha+beta*exp(-gamma/Kn));

end


function out = eta_fm_repulsive(psiE, psiI)

dv = 1e-4;
dr = 1e-3;

rmin = log10(1);

if (psiI == 0)
	rc = 1;
else
	rc = newton_raphson_maxima_repulsive(psiE, psiI);
end

phi_c = pot(rc, psiE, psiI);
rmax=log10(1d2);

vmin=log10(phi_c^0.5);
vmax=log10(1d2); % represents infinity

out = 0;
v = vmin;
v_prev = vmin - dv;
while (v <= vmax)
	b = 1e99;
	r = rmax;
	linv = (10^v);
	while (r >= rmin)
		linr = 10 ^ r;
		b = min(...
            (linr*linr*(1-(1/linv/linv*pot(linr,psiE,psiI))))^0.5,...
            b);
		r = r - dr;
	end
	term = 2*linv*linv*linv*exp(-linv*linv)*b*b*...
        (10^v - 10^v_prev);
	out = out + term;

    v_prev = v;
    v = v + dv;
end

end



function r2 = newton_raphson_maxima_repulsive(psiE, psiI)

r1 = 100;
residue = 1;
while (residue > 1e-6)  % iterate until convergence
	r2 = r1 - f(r1, psiE, psiI) / df(r1,psiE, psiI);
	residue = abs((r2 - r1) / r1);
	r1 = r2;
end

end

function out = f(x, psiE, psiI)

out = psiE*(x^5) - 2*psiE*(x^3) + 2*psiI*x*x + psiE*x - psiI;

end

function out = df(x, psiE, psiI)

out = 5*psiE*(x^4) - 6*psiE*x*x + 4*psiI*x + psiE;

end

function out = pot(r, psiE, psiI)

if (psiI == 0e0)
    out = -psiE/r;
else
    out = -psiE/r - psiI/2/r/r/(r*r-1);
end

end




function out = eta_fm_attractive(psiE,psiI)

dv = 1e-3; dr = 5e-4;

rmin = 1 + dr;
rmax = 100;
vmin = -5;
vmax = 2;

rvec = rmax:-dr:rmin;
po = potential(rvec, psiE, psiI);

vvec = (vmin:dv:vmax)';
linv = 10 .^ vvec;

b = zeros(length(vvec), 1);
for vv=1:length(vvec)
    b2 = rvec.^2 .* (1 - (po ./ (linv(vv)^2)));
    b2 = b2(b2 >= 0);
    b2 = min(sqrt(b2));

    if isempty(b2); b2 = 1; end
    b(vv) = max(b2, 1);
end

term = 2 .* linv .^ 3 .* exp(-(linv .^ 2)) .* b .^ 2 .* ...
    (linv - 10 .^ (vvec - dv));

out = sum(term);

end



function out = potential(r,psiE,psiI)

out = -psiE./r - psiI./2./r./r./(r.*r-1);

end




function si = simpsons_cont(a,b,n,psiE,psiI)

ti_n = trapezoidal_cont(a,b,n,psiE,psiI);
ti_n_1 = trapezoidal_cont(a,b,n-1,psiE,psiI);

si = ti_n + 1.0/3.0*(ti_n - ti_n_1);

end

function ti = trapezoidal_cont(a,b,n,psiE,psiI)

delx=(b-a)/n;

ti = 0.0;
for i = 2:n
    ti=ti+func(a+i*delx,psiE,psiI);
end

ti=0.5*delx*(func(a,psiE,psiI) + 2*ti + func(b,psiE,psiI));

end

function out = func(x,psiE,psiI)

term = -psiE*x-psiI/2*(x^4)/(1-x*x);
out  = exp(term);

end

