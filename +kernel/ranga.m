

% RANGA A function to determine the average charge on particles.
%   ORIGINAL AUTHOR: Ranganathan Gopalakrishnan
%   MODIFIED BY: Timothy Sipkens, 2023-02-08

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
p = 1.01325e5; T=2.98e2; mwt=2.9e-2; mu=1.81e-5; % Pa, K, kg/mole, Pa.s
alpha = 1.257e0; beta = 4e-1; gamma=1.1e0;

% ion properties
MODEL = 'SIPKENS'
ni=1;

% particle properties
PA_nd=3.14159e0; Rs_nd=1e0; % dimless, dimless
dielec=13.5; georatio=10;
STARTnt=1.0d12; ENDnt=1.0d13; TAG='06P76D12';
NPMAX=10

	
if strcmp(MODEL, 'MOHSEN')
    
    Mwt_ion=0.1
    zp_ion=1.9d-4 % kg/mole, m^2/s/V

elseif strcmp(MODEL, 'VOHRAA')
    
    mwt_ion=0.05
    zp_ion=1.9d-4 % kg/mole, m^2/s/V

elseif strcmp(MODEL, 'SIPKENS')
    
    mwt_ion=0.109
    zp_ion=1.4d-4 % kg/mole, m^2/s/V
	
end
	
		
% Mobility equivalent diameter
d_me_array = logspace(log10(4), log10(2e3), 120);
d_me_array = d_me_array([1,40,50,100]);

meancharge_d =[];
for i = 1:length(d_me_array)
    
    d_me = d_me_array(i) * 1e-9;
    
    disp('Running:')
    disp(['d_me = ', num2str(d_me * 1e9)]);
    
    % calculation of gas properties
    R_gas=R_univ/mwt; % J/kg/K
    pho=p/(R_gas*T); % kg/m^3
    m_gas=mwt/N_av; % mass of gas molecule, kg
    mts=(8*k*T/pi/m_gas)^0.5; % m/s
    mfp=mu/(0.499*pho*mts); % hard sphere mean free path, m
    
    % calculation of primary particle radius of aggregate with the same mobility diameter
    R_me=d_me/2;
	    
    rp = R_me;
    PA=PA_nd*(rp^2); % m^2
    Rs=Rs_nd*(rp); % m
    
    % calculation of ion properties
    m_i=mwt_ion/N_av;
    f_i=ni*e/zp_ion;
    
    if (d_me*1e9 <= 10)
	    NPMAX = 5;
    elseif (d_me*1e9 <= 50)
	    NPMAX = 10;
    elseif (d_me*1e9 <= 400)
	    NPMAX = 20;
    elseif (d_me*1e9 <= 800)
	    NPMAX = 50;
    elseif (d_me*1e9 <= 1200)
	    NPMAX = 100;
    else
	    NPMAX = 300;
    end
    
    psiEarray = [];
    psiIarray = [];
    etafarray = [];
    etacarray = [];
    kndarray = [];
    collkernel = [];
    for np=0:NPMAX
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
                eta_fm_repulsive;
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
        kndarray(np + 1) = KnD;
        collkernel(np + 1) = H*L_H3*f_i/m_i;
	    
    end
    
    
    for np=0:NPMAX
	    
	    collkernelratio(np + 1) = collkernel(np + 1)/collkernel(1);
			    
    end
    
    
    nt = STARTnt;
    
    meancharge_nt =[];
    while (nt <= ENDnt)
        meancharge=0e0;
        ssfrac(1) = exp(-collkernel(1) * nt);
        sscheck = ssfrac(1);
	    
        for np=1:NPMAX
	        
	        ssfrac(np + 1) = 0;
        
            for nj = 0:np
	            
		        dummy1 = diff_product(collkernelratio, NPMAX, nj, np);
		        ssfrac(np + 1) = ssfrac(np + 1) + exp(-collkernel(nj + 1)*nt)/dummy1;
            
            end
	        dummy2 = product_term(collkernelratio, NPMAX, np);
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

out = 1;
for ni=1:np-1
	out = out*collkernelratio(ni + 1);
end

end


function out = Cc(Kn,alpha,beta,gamma)

out = 1+Kn*(alpha+beta*exp(-gamma/Kn));

end

function out = eta_fm_attractive(psiE,psiI)

dv=1e-3; dr=5e-4;

out=0e0;

rmin=1+dr;
rmax=1e2;
vmin=-5e0;
vmax=2e0;

v=vmin;
v_prev=vmin-dv;
while (v <= vmax)
	b=1e99;
	r=rmax;
	j=0;
	linv=(10^v);
    while (r >= rmin)
        b2 = r * r * (1-(1/linv/linv*potential(r,psiE,psiI)));
        
        if (b2 >= 0e0)
	        b=min(b2^0.5,b);
	        j=1;
        end
        
        r = r-dr;
    end

	if (j == 0)
		b = 1;
	else
	    b = max(b, 1);
	end
			
	term = 2*linv*linv*linv*exp(-linv*linv)*b*b* ...
        (10 ^ v - 10 ^ v_prev);
	
	out = out + term;
    v_prev = v;
    v = v + dv;
end

end



function out = potential(r,psiE,psiI)

out = -psiE/r-psiI/2/r/r/(r*r-1);

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

ti=0.5*delx*(func(a,psiE,psiI)+2*ti+func(b,psiE,psiI));

end

function out = func(x,psiE,psiI)

term=-psiE*x-psiI/2*(x^4)/(1-x*x);
out = exp(term);

end

