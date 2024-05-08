function [ionic_current, J_Na, J_K,  J_p, J_l] = ionicupdatethetag(m,n,h,p,Z1_K,Z1_Na,V_l,V_1,d, V_rest,theta,g,gf,time)
E = V_1 + V_rest;
mave = m;
nave = n;
have = h;
pave = p;
gcostheta = (1+g)*cosd(theta);
powerofnave = 2;
powerofmave = 2; 

g_l = 30.3e-6; % S per cm^2
    P_Na = 0.008; P_p = 0.00054; P_K = 0.0012; % cm/s

h=g;
J_K = P_K.*nave.^powerofnave.*Z1_K; % A/cm^2
J_Na = P_Na.*mave.^powerofmave.*have.*Z1_Na;
J_p = P_p.*pave.^2.*Z1_Na;
J_l = g_l.*(V_1 - V_l); % Note sub ell
J_FH = J_K + J_Na + J_p + J_l;
ionic_current = (1/gcostheta)*pi.*d.*J_FH;% A/cm
ionic_current = ((1+h)/(1+h*cos(2*pi*gf*time)))*ionic_current;


        