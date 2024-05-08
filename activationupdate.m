function [mnew, nnew, hnew, pnew] = activationupdate(m, n, h, p, V_1, Q, delT)
%% This file corrects the betah formula and bases it upon Frankenhaeuser 1960 rather
%% than on Reutskiy. 

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

% Aalpham = 0.36*10^(exponent); %e6; % per s
% Balpham = 22e-3; 
% Calpham = 3e-3; % volts 
% 
% Abetam = 0.4*10^(exponent); %0.4e6; 
% Bbetam = 13e-3; 
% Cbetam = 20e-3;
% 
% Aalphah = 0.1*10^(exponent); %e6; 
% Balphah = -10e-3; 
% Calphah = 6e-3;
% 
% Abetah = 4.5*10^3; %0.005*10^(exponent);%See F-H 1969 Equation 19. 4.5*10^(exponent); %e3; % per s
% Bbetah = 45e-3;%32e-3; %45e-3; 
% Cbetah = 10e-3;% 10e-3;  % volts 
% 
% Aalphan = 0.02*10^(exponent); %e6; 
% Balphan = 35e-3; 
% Calphan = 10e-3;
% 
% Abetan = 0.05*10^(exponent); %e6; 
% Bbetan = 10e-3; 
% Cbetan = 10e-3;
% 
% Aalphap = 0.006*10^(exponent); %e6; 
% Balphap = 40e-3; 
% Calphap = 10e-3;
% 
% Abetap = 0.09*10^(exponent); %e6; 
% Bbetap = -25e-3; 
% Cbetap = 20e-3;

            
alpham = Aalpham.*(V_1 - Balpham)...
    ./(1 - exp((Balpham - V_1)./Calpham));
alphan = Aalphan.*(V_1 - Balphan)...
    ./(1 - exp((Balphan - V_1)./Calphan));
alphap = Aalphap.*(V_1 - Balphap)...
    ./(1 - exp((Balphap - V_1)./Calphap));
alphah = Aalphah.*(Balphah-V_1)...
    ./(1 - exp((V_1 - Balphah)./Calphah));

betam = Abetam.*(Bbetam-V_1)...
    ./(1 - exp((V_1 - Bbetam)./Cbetam));
betan = Abetan.*(Bbetan-V_1)...
    ./(1 - exp((V_1 - Bbetan)./Cbetan));
betap = Abetap.*(Bbetap-V_1)...
    ./(1 - exp((V_1 - Bbetap)./Cbetap));
betah = (Abetah./(1 + exp((Bbetah - V_1)./Cbetah))); % per s
% betah = Abetah.*(V_1 - Bbetah)...
%   ./(1 - exp((Bbetah - V_1)./Cbetah));

expdecay_m = exp(-Q.*(alpham+betam).*delT);

m_next = m.*expdecay_m + (1-expdecay_m).*alpham./(alpham+betam);
expdecay_n = exp(-Q.*(alphan+betan).*delT);
n_next = n.*expdecay_n + (1-expdecay_n).*alphan./(alphan+betan);
expdecay_h = exp(-Q.*(alphah+betah).*delT);
h_next = h.*expdecay_h + (1-expdecay_h).*alphah./(alphah+betah);
expdecay_p = exp(-Q.*(alphap+betap).*delT);
p_next = p.*expdecay_p + (1-expdecay_p).*alphap./(alphap+betap);


% % Reutskiy Formula

% m_next = m.*((2 - delT.*Q.*(alpham + betam))./(2 + delT.*Q.*(alpham + betam))) + ((2.*delT.*alpham.*Q)./(2 + delT.*Q.*(alpham + betam)));
% n_next = n.*((2 - delT.*Q.*(alphan + betan))./(2 + delT.*Q.*(alphan + betan))) + ((2.*delT.*alphan.*Q)./(2 + delT.*Q.*(alphan + betan)));
% h_next = h.*((2 - delT.*Q.*(alphah + betah))./(2 + delT.*Q.*(alphah + betah))) + ((2.*delT.*alphah.*Q)./(2 + delT.*Q.*(alphah + betah)));
% p_next = p.*((2 - delT.*Q.*(alphap + betap))./(2 + delT.*Q.*(alphap + betap))) + ((2.*delT.*alphap.*Q)./(2 + delT.*Q.*(alphap + betap)));
% 
mnew = m_next;
nnew = n_next;
hnew = h_next;
pnew = p_next;