function [ionic_current, J_Na, J_K,  J_p, J_l] = ionicupdatetheta(m,n,h,p,Z1_K,Z1_Na,V_l,V_1,d, V_rest,theta)
E = V_1 + V_rest;
mave = m;
nave = n;
have = h;
pave = p;

powerofnave = 2;
powerofmave = 2; 

g_l = 30.3e-6; % S per cm^2
    P_Na = 0.008; P_p = 0.00054; P_K = 0.0012; % cm/s


J_K = P_K.*nave.^powerofnave.*Z1_K; % A/cm^2
J_Na = P_Na.*mave.^powerofmave.*have.*Z1_Na;
J_p = P_p.*pave.^2.*Z1_Na;
J_l = g_l.*(V_1 - V_l); % Note sub ell
J_FH = J_K + J_Na + J_p + J_l;
ionic_current = (1/cosd(theta))*pi.*d.*J_FH; % A/cm

% include randomness at the ionic current level
scale = 0;% suppress the noise. .8e-5; %1e-5 is a good value for scale --> show EEG like response.
ionic_current = ionic_current + scale*randn(1,1);



        