function [alpham, alphan, alphah, alphap, betam, betan, betah, betap] = actinit(E_rest)
%% The following follow Dr. Snider's values
Aalpham = 3.6e5; %e6; % per s
Balpham = 22e-3; 
Calpham = 3e-3; % volts 

Abetam = 4e5; %0.4e6; 
Bbetam = 13e-3; 
Cbetam = 20e-3;

Aalphah = 1e5; %e6; 
Balphah = -10e-3; 
Calphah = 6e-3;

Abetah = 4.5*10^3; %0.005*10^(exponent);%See F-H 1969 Equation 19. 4.5*10^(exponent); %e3; % per s
Bbetah = 45e-3;%32e-3; %45e-3; 
Cbetah = 10e-3;% 10e-3;  % volts 

Aalphan = 2e4; %e6; 
Balphan = 35e-3; 
Calphan = 10e-3;

Abetan = 5e4; %e6; 
Bbetan = 10e-3; 
Cbetan = 10e-3;

Aalphap = 6e3; %e6; 
Balphap = 40e-3; 
Calphap = 10e-3;

Abetap = 9e4; %e6; 
Bbetap = -25e-3; 
Cbetap = 20e-3;


alpham = Aalpham*(E_rest - Balpham)/(1 - exp((Balpham-E_rest)/Calpham));
 alphan = Aalphan*(E_rest - Balphan)/(1 - exp((Balphan-E_rest)/Calphan));
 alphap = Aalphap*(E_rest - Balphap)/(1 - exp((Balphap-E_rest)/Calphap));
 alphah = Aalphah*(Balphah-E_rest)/(1 - exp((E_rest - Balphah)/Calphah));
 betam = Abetam*(Bbetam-E_rest)/(1 - exp((E_rest - Bbetam)/Cbetam));
 betan = Abetan*(Bbetan-E_rest)/(1 - exp((E_rest - Bbetan)/Cbetan));
 betap = Abetap*(Bbetap-E_rest)/(1 - exp((E_rest - Bbetap)/Cbetap));
 betah = (Abetah./(1 + exp((Bbetah - E_rest)./Cbetah))); % per s