function [ionic_current, J_Na, J_K,  J_p, J_l] = ionicupdatevth(m,n,h,p,Z1_K,Z1_Na,V_l,V_1,d, V_rest,th)
E = V_1 + V_rest;
mave = m;
nave = n;
have = h;
pave = p;

powerofnave = 2;
powerofmave = 2; 

g_l = 30.3e-6; % S per cm^2
    P_Na = 0.008; P_p = 0.00054; P_K = 0.0012; % cm/s

scaling = 1e-7; % it is the ratio of the axon diameter to the tract diameter.
% scaling should be very small. if it is 1e-5, we see some noisy action
% potentials. but very noisy. so we reduce it further.
J_K = P_K.*nave.^powerofnave.*Z1_K; % A/cm^2
J_Na = P_Na.*mave.^powerofmave.*have.*Z1_Na;
J_p = P_p.*pave.^2.*Z1_Na;
J_l = g_l.*(V_1 - V_l); % Note sub ell
J_FH = J_K + J_Na + J_p + J_l;
ionic_current = (1/(cos(th)^2)).*pi.*d.*J_FH;% - noiseon*scaling*randn(1); % A/cm