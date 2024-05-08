%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Author: Aman Chawla     %%
%%% Date: September 1, 2015 %%
%%% Rev: October 21, 2015%%
%%% Rev: November 1, 2015 (double stimulation) %%
%%% Rev: Nov 3, 2015 (single stimulation) %%
%%% Rev: Dec 27, 2015 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [outputText, axon] = packagedCNNDec24_2016v6(N, d, ell, L, delT, timesteps, endlength, temperature, Ranvier_array, r0, cm, gm, cnd, gnd, rf, InjectionDuration, injectnumberarray, stimposarray, delaydeltarray, Eonoff, Estrength, K, theta)

bundle1 = bUndle; % create bundle1 as a bundle array
bundle1.N = N; %6; % Number of axons in bundle
bundle1.d = d; %1e-3;
bundle1.ell = ell; %2.5e-4;
bundle1.L = L; %.2;
bundle1.delX = bundle1.L/10;
bundle1.delT = delT; %2e-6;
bundle1.timesteps = timesteps; %300; 
bundle1.duration = bundle1.timesteps*bundle1.delT;
bundle1.endlength = endlength; %200; 
bundle1.temperature = temperature; %24;% temperature in degrees Celsius
bundle1.Ranvier = Ranvier_array; %20:20:180;
LTminus1 = bundle1.timesteps;
bundle1.K = K;
LYY = LTminus1;   % Both symbols are used in the code
bundle1.r0 = r0; %1e5; % interfiber medium resistance
bundle1.rf = rf; %100/(pi*(bundle1.d/2)^2); % Ohm/cm (Waxman); see Snider's page

alpha_0 = bundle1.r0/bundle1.rf ; 
alpha_N = alpha_0/(1 + bundle1.N*alpha_0);

bundle1.cm = cm; %1.87e-11; % F/cm (both)
bundle1.gm = gm; %5.6e-9; % per Ohm-cm (both)
bundle1.cnd = cnd; %3.14e-9; % F/cm (from Waxman - Reutskiy has wrong dimensions)
bundle1.gnd = gnd; %0;  % Not used in Reutskiy, but included in cde here

egm = exp(-bundle1.gm*bundle1.delT/bundle1.cm);
egnd = exp(-bundle1.gnd*bundle1.delT/bundle1.cnd);

%InjectionDuration = 1e-6; % sec. (Reutskiy page 4)
bundle1.injectnumber = injectnumberarray(1); %2; % axon to be injected in the bundle.

% Specified the electric field based interaction gain
gain = 0;%5e-11; %1 is high, 0.001 is low , etc. 

Q = 3^((bundle1.temperature-20)/10); 

bundle1.F = 96485; %Coulombs per mole
bundle1.R = 8.3145; %Joules per Kelvin - mole
T_ab = bundle1.temperature+273.15; %Kelvin temperature for 25 degrees Celsius
gamma = bundle1.F./(bundle1.R*T_ab); % per volt

% Specify the ionic concentrations inside and outside
bundle1.Nao = 114.5e-6; % moles per cm^3
bundle1.Nai = 13.74e-6; 
bundle1.Ko = 2.5e-6; 
bundle1.Ki = 120e-6;

% Sodium Nernst potential
NaNernstPot  = (1/gamma)*log(bundle1.Nao/bundle1.Nai);



% The above, so far, are the common parameters. 

 factor = 4/cosd(theta);%input('What value of injection intensity factor would you like to use? Suggest value is 1. ')
    InjectionIntensity = factor*(-1e-9/bundle1.ell); % A/cm (Reutskiy page 4)
    CoulombsInjected = InjectionIntensity*bundle1.ell*InjectionDuration;
    E_rest = -70e-3; % volts
    V_l = -0.026*3.5; % volts. adjusted so that the initial ionic current (unstimulated condition) is 0
    g_l = 30.3e-6; % S per cm^2
    P_Na = 0.008; P_p = 0.00054; P_K = 0.0012; % cm/s
    % compute alphas and betas starting values
[alpham, alphan, alphah, alphap, betam, betan, betah, betap] = actinit(E_rest);
 % create axon objects. axon is an array of N axons.
    axon = aXon;
    
     E_rest =-70e-3;
 
 for i = 1:bundle1.N % axon index
    
   
    axon(i).number = i;
    axon(i).stimpos = stimposarray(1);
    
    
    
    
%     % allocate space
    axon(i).J_ex = zeros(bundle1.endlength, LTminus1);

    
     
    
    
    
    % put starting values
  
    axon(i).E(1:bundle1.endlength,1) = E_rest.*ones(bundle1.endlength,1); % since V = E - E_r
    axon(i).V(1:bundle1.endlength,1) =  0.*ones(bundle1.endlength,1);
    
   
    axon(i).m(1:bundle1.endlength,1) = (alpham/(alpham+betam));%.*ones(1:bundle1.endlength,1);
    axon(i).n(1:bundle1.endlength,1) = (alphan/(alphan+betan));%.*ones(1:bundle1.endlength,1);
    axon(i).h(1:bundle1.endlength,1) = (alphah/(alphah+betah));%.*ones(1:bundle1.endlength,1);
    axon(i).p(1:bundle1.endlength,1) = (alphap/(alphap+betap));%.*ones(1:bundle1.endlength,1);
 
%      [axon(i).mcoeff(1:bundle1.endlength,1), axon(i).ncoeff(1:bundle1.endlength,1),  axon(i).hcoeff(1:bundle1.endlength,1),axon(i).pcoeff(1:bundle1.endlength,1)] = activationupdate(axon(i).m(1:bundle1.endlength,1), axon(i).n(1:bundle1.endlength,1), axon(i).h(1:bundle1.endlength,1), axon(i).p(1:bundle1.endlength,1), axon(i).Vcoeff(1:bundle1.endlength,1), Q, bundle1.delT);
   


    axon(i).Zfactor(1:bundle1.endlength,1)=exp(-gamma.*axon(i).E(1:bundle1.endlength,1));
    axon(i).Z_Na(1:bundle1.endlength,1) = gamma*bundle1.F.*axon(i).E(1:bundle1.endlength,1).*((bundle1.Nao.*axon(i).Zfactor(1:bundle1.endlength,1) - bundle1.Nai)./(1-axon(i).Zfactor(1:bundle1.endlength,1)));
    axon(i).Z_K(1:bundle1.endlength,1) = gamma*bundle1.F.*axon(i).E(1:bundle1.endlength,1).*((bundle1.Ko.*axon(i).Zfactor(1:bundle1.endlength,1) - bundle1.Ki)./(1-axon(i).Zfactor(1:bundle1.endlength,1)));
   
            
    [axon(i).J_ion(1:bundle1.endlength,1), axon(i).J_Na(1:bundle1.endlength,1), axon(i).J_K(1:bundle1.endlength,1), axon(i).J_p(1:bundle1.endlength,1), axon(i).J_l(1:bundle1.endlength,1)] = ionicupdatetheta(axon(i).m(1:bundle1.endlength,1),axon(i).n(1:bundle1.endlength,1),axon(i).h(1:bundle1.endlength,1),axon(i).p(1:bundle1.endlength,1),axon(i).Z_K(1:bundle1.endlength,1),axon(i).Z_Na(1:bundle1.endlength,1),V_l,axon(i).E(1:bundle1.endlength,1),bundle1.d, E_rest, theta);
               
       
       
%     axon(i).Ecoeff(1:bundle1.endlength,1)  = axon(i).Vcoeff(1:bundle1.endlength,1) + E_rest;
%      
%     axon(i).Zfactorcoeff(1:bundle1.endlength,1)=exp(-gamma.*axon(i).Ecoeff(1:bundle1.endlength,1));
%     axon(i).Z_Nacoeff(1:bundle1.endlength,1) = gamma*bundle1.F.*(axon(i).Ecoeff(1:bundle1.endlength,1)).*((bundle1.Nao.*axon(i).Zfactor(1:bundle1.endlength,1) - bundle1.Nai)./(1-axon(i).Zfactor(1:bundle1.endlength,1)));
%     axon(i).Z_Kcoeff(1:bundle1.endlength,1) = gamma*bundle1.F.*(axon(i).Ecoeff(1:bundle1.endlength,1)).*((bundle1.Ko.*axon(i).Zfactor(1:bundle1.endlength,1) - bundle1.Ki)./(1-axon(i).Zfactor(1:bundle1.endlength,1)));
%    
%             
%     [axon(i).J_ioncoeff(1:bundle1.endlength,1), axon(i).J_Nacoeff(1:bundle1.endlength,1), axon(i).J_Kcoeff(1:bundle1.endlength,1), axon(i).J_pcoeff(1:bundle1.endlength,1), axon(i).J_lcoeff(1:bundle1.endlength,1)] = ionicupdate(axon(i).mcoeff(1:bundle1.endlength,1),axon(i).ncoeff(1:bundle1.endlength,1),axon(i).hcoeff(1:bundle1.endlength,1),axon(i).pcoeff(1:bundle1.endlength,1),axon(i).Z_Kcoeff(1:bundle1.endlength,1),axon(i).Z_Nacoeff(1:bundle1.endlength,1),V_l,(axon(i).Vcoeff(1:bundle1.endlength,1)+E_rest),bundle1.d, E_rest);

 end
     
 % Prepare for current injection
% injectnumberarray = [injectnumber1, injectnumber2, injectnumber3, injectnumber4];
% stimposarray = [stimpos1, stimpos2, stimpos3, stimpos4];
% delaydeltarray = [delaydelta, delaydelta2, delaydelta3];
for my_count = 1:bundle1.N
    axon(injectnumberarray(my_count)).stimpos = stimposarray(my_count);
 % inject current
    if my_count >=2
%         my_count
    axon(injectnumberarray(my_count)).J_ex(axon(injectnumberarray(my_count)).stimpos, [10+delaydeltarray(my_count-1):delaydeltarray(my_count-1)+floor(InjectionDuration/bundle1.delT)+10,1060+delaydeltarray(my_count-1):delaydeltarray(my_count-1)+floor(InjectionDuration/bundle1.delT)+1060]) = factor*InjectionIntensity;
    else
    axon(injectnumberarray(my_count)).J_ex(axon(injectnumberarray(my_count)).stimpos, 10:floor(InjectionDuration/bundle1.delT)+10) = factor*InjectionIntensity;
    end
end

 % Prepare to apply field
        V_app = 12; % Volts DC
        R_T1 = -1*V_app/(Estrength/bundle1.ell); % Ohms was 2.5e-11
        R_12 = -1*V_app/(Estrength/bundle1.ell); % Ohms
        R_23 = -1*V_app/(Estrength/bundle1.ell); % Ohms
        R_34 = -1*V_app/(Estrength/bundle1.ell); % Ohms
        
%         % Apply field
%         axon(injectnumber1).J_ex = axon(injectnumber1).J_ex + Eonoff.*(V_app/R_T1);
%         axon(injectnumber2).J_ex = axon(injectnumber2).J_ex + Eonoff.*(V_app/R_12);
%         axon(injectnumber3).J_ex = axon(injectnumber3).J_ex + Eonoff.*(V_app/R_23);
%         axon(injectnumber4).J_ex = axon(injectnumber4).J_ex + Eonoff.*(V_app/R_34);



         
         % BLOCK tridiagonal matrix for internodal segments

         a = (alpha_N-1)/(2*bundle1.rf*(bundle1.delX)^2);
         b = bundle1.cm/bundle1.delT + (1-alpha_N)/(bundle1.rf*(bundle1.delX)^2) + bundle1.gm/2;
         c = a;
         ee = alpha_N/(2*bundle1.rf*(bundle1.delX)^2);
         % create the Amd subblock
         
     B = diag(b.*ones(bundle1.N, 1)) + ((-2*ee).*ones(bundle1.N, bundle1.N)-diag((-2*ee).*ones(bundle1.N,1)));
% repair B
for i = 1:N
    for j = 1:N
        if i==j
            B(i,j) = (B(i,j)- (bundle1.cm/bundle1.delT) - (1/(bundle1.rf*(bundle1.delX)^2)) - (bundle1.gm/2))*K(i,j) + (bundle1.cm/bundle1.delT + (1/(bundle1.rf*(bundle1.delX)^2)) + bundle1.gm/2);
        else
            B(i,j) = B(i,j)*K(i,j);
        end
    end
end

         
         % create the Asubd subblock
         
              A = diag(a.*ones(bundle1.N, 1)) + ((ee).*ones(bundle1.N, bundle1.N)-diag((ee).*ones(bundle1.N,1)));

             %repair A
             for i = 1:N
                for j = 1:N
                    if i==j
                        A(i,j) = (A(i,j) + (1/(2*bundle1.rf*(bundle1.delX)^2)))*K(i,j) - (1/(2*bundle1.rf*(bundle1.delX)^2));
                    else
                        A(i,j) = A(i,j)*K(i,j);
                    end
                end
             end


         
         % create the Asupd subblock
         
              C = diag(c.*ones(bundle1.N, 1)) + ((ee).*ones(bundle1.N, bundle1.N)-diag((ee).*ones(bundle1.N,1)));
  
              
              % repair C
              
              for i = 1:N
                  for j = 1:N
                      if i == j
                        C(i,j) = (C(i,j) + (1/(2*bundle1.rf*(bundle1.delX)^2)))*K(i,j) - (1/(2*bundle1.rf*(bundle1.delX)^2));
                      else 
                          C(i,j) = C(i,j)*K(i,j);
                      end
                  end
              end
              
         % construct the block tridiagonal matrix
         leftmatrix = blktridiag(B, A, C, bundle1.endlength);
         
  
        
        % Repair Ranvier-node parameters
    
       
        
         b_nd = bundle1.cnd/bundle1.delT + (1-alpha_N)/(bundle1.rf*(bundle1.delX)^2);
        a_nd = a;
         c_nd = c;
         e_nd = ee;
       
    
         % create the Amd subblock
         
            B_nd = diag(b_nd.*ones(bundle1.N, 1)) + ((-2*e_nd).*ones(bundle1.N, bundle1.N)-diag((-2*e_nd).*ones(bundle1.N,1)));

         % Repair B_nd
         for i = 1:N
             for j = 1:N
                 if i==j
                     B_nd(i,j) = (B_nd(i,j) - (bundle1.cnd/bundle1.delT) -(1/(bundle1.rf*(bundle1.delX)^2)))*K(i,j) + ((bundle1.cnd/bundle1.delT) + (1/(bundle1.rf*(bundle1.delX)^2)));
                 else
                     B_nd(i,j) = B_nd(i,j)*K(i,j);
                 end
             end
         end
         
         % create the Asubd subblock
         A_nd = A;
         
         % create the Asupd subblock
         
         C_nd = C;
         
            for i = 2:bundle1.endlength % The first position should not be a node
                if ismember(i, bundle1.Ranvier)
                    leftmatrix(((i-1)*bundle1.N)+1:i*bundle1.N,((i-1)*bundle1.N)-bundle1.N+1:((i-1)*bundle1.N)) = A_nd;
                    leftmatrix(((i-1)*bundle1.N)+1:i*bundle1.N,((i-1)*bundle1.N)+1:i*bundle1.N) = B_nd;
                    leftmatrix(((i-1)*bundle1.N)+1:i*bundle1.N,(i*bundle1.N)+1:((i+1)*bundle1.N)) = C_nd;
                else % do nothing
                end
            end
           
            
            bundle1.a = a;
            bundle1.b = b; 
            bundle1.c = c;
            bundle1.ee = ee;
            bundle1.b_nd = b_nd;


            outputLines = {};  % Initialize cell array to store text lines
   
    % Add your processing code here
   

% Begin time stepping
   for time = 1:(LTminus1-1)  %Recall initial value of V_1 is set to zero.
       t=time; % Both notations are used 
      % disp(t)
        outputLines{end+1} = sprintf('Processing input: %f', t);
   % Step (i) of Reutskiy algorithm
           for i = 1:bundle1.endlength
        for j=1:bundle1.N
            if ismember(i,bundle1.Ranvier)
       % Update the ionic current
           
            axon(j).E(i,t) = axon(j).V(i,t) + E_rest;
            axon(j).Zfactor(i,t)=exp(-gamma*axon(j).E(i,t));
            axon(j).Z_Na(i,t) = gamma*bundle1.F*axon(j).E(i,t)*((bundle1.Nao*axon(j).Zfactor(i,t) - bundle1.Nai)/(1-axon(j).Zfactor(i,t)));
            
            axon(j).Z_K(i,t) = gamma*bundle1.F*axon(j).E(i,t)*((bundle1.Ko*axon(j).Zfactor(i,t) - bundle1.Ki)/(1-axon(j).Zfactor(i,t)));
            axnum = j;
            nodenum = i;
            presenttime = t;
            [axon(j).J_ion(i,t), axon(j).J_Na(i,t), axon(j).J_K(i,t), axon(j).J_p(i,t), axon(j).J_l(i,t)]  = ionicupdatetheta(axon(j).m(i,t),axon(j).n(i,t),axon(j).h(i,t),axon(j).p(i,t),axon(j).Z_K(i,t),axon(j).Z_Na(i,t),V_l,(axon(j).V(i,t)+E_rest),bundle1.d, E_rest, theta);
          axnum = j;
%            [axon(axnum).J_ioncoeff(i,time), axon(axnum).J_Nacoeff(i,time), axon(axnum).J_Kcoeff(i,time), axon(axnum).J_pcoeff(i,time), axon(axnum).J_lcoeff(i,time)] = ionicupdate(axon(axnum).mcoeff(i,time),axon(axnum).ncoeff(i,time),axon(axnum).hcoeff(i,time),axon(axnum).pcoeff(i,time),axon(axnum).Z_Kcoeff(i,time),axon(axnum).Z_Nacoeff(i,time),V_l,(axon(axnum).Ecoeff(i,time)),bundle1.d, E_rest);
% Include the current injections from fields.
% if axnum==1
%     otheraxons  = 2;
% else 
%     otheraxons = 1; 
% end

% cycle through every other node on axons other than the present one. 
% from each such node, gather its field at a distance , specifically, at
% the distance which intervenes between that node and the present node. 
if gain>0
for maxon = setdiff(1:N, [axnum])
    for node = i %[i-20,i,i+20] %ignoring some `other' nodes, for saving time
  if (node>1)&&(node<300)
%       disp('hi6')
            if t ==1
          
            result=mainporeAug172020v2(axon(maxon).m(node,t),axon(maxon).n(node,t),axon(maxon).h(node,t),axon(maxon).p(node,t),t,0,0,0,0);
            oldresult = result; 
            fieldzero = 1e-40;
            myfield(maxon,node) = fieldzero;
            else
            result=mainporeAug172020(axon(maxon).m(node,t),axon(maxon).n(node,t),axon(maxon).h(node,t),axon(maxon).p(node,t),t,oldresult,axon(maxon).m(node,t-1),axon(maxon).n(node,t-1),axon(maxon).h(node,t-1),axon(maxon).p(node,t-1));
          %  [c,h] = contourf(real(result{1,1}), real(result{1,2}),real(result{1,3}), 20);
%             oldsourcelowv = result{1,4}; % this is the open channel configuration used 
%             % to compute the result i.e., the field shown. in the next
%             % iteration, this will be merely modified slightly to save on
%             % computational time. 
           oldresult = result;
%             colorbar
%t1=clabel(c,h); % acquire the text and line objects created to label the plot. From these extract the weakest field value
myfield(maxon, node) = result{1,4}; %min(min(t1)); %assuming t1 is 2 dimensional
%        caxis('auto');
%          %caxis([200 1000000])
% %caxis([0 realmax])
%        axis([2.5e-6 5.5e-6 0.1e-6 3.1e-6]);
%        xlabel('\rho (meters)');
% 
%        ylabel('z (meters)');
%        title(num2str(t));
%         drawnow
            end 
  end
    end
end
% slopek = 0.2e4; %slope of i-v relationship in mA/(mV times cm^2)
% memthick = 10e-9; % membrane thickness in meters
% mynetcurrent = sum(sum(myfield))*slopek*memthick;
%          
%  axon(j).J_ion(i,t)=axon(j).J_ion(i,t)+mynetcurrent;   
 
 slopek = 1; %slope of K i-v relationship in mA/(mV times cm^2)
memthick = 1; % membrane thickness in meters
%gain = 1;
mynetcurrent =  gain.*myfield(maxon, node);%*slopek*memthick; assuming conductivity of sea water is 5S/m
        % size(axon(axnum).J_ioncoeff(i,time))
        % if at the present time the net nodal current is outward,
        % mynetcurrent should be subtracted
        % else if at the present time the net ionic current is inward,
        % my netcurrent should be added.
        if axon(axnum).J_ion(i,time) >0 % inward
 axon(axnum).J_ion(i,time)=axon(axnum).J_ion(i,time)+mynetcurrent; 
        elseif axon(axnum).J_ion(i,time)<=0
            axon(axnum).J_ion(i,time) = axon(axnum).J_ion(i,time)-mynetcurrent;
        end
        end
            
            end
        end
         
           end     
                         
            % Create the right vector
             dee = bundle1.cm/bundle1.delT - (1-alpha_N)/(bundle1.rf*((bundle1.delX)^2)) - (bundle1.gm/2);
            dee_nd = bundle1.cnd/bundle1.delT + (alpha_N - 1)/(bundle1.rf*(bundle1.delX)^2);
            Dd = rightcreateReutskiy1v3(axon, bundle1, dee, dee_nd, t);
            
           
             % Step (ii) of Reutskiy algorithm.
      axon = ReutskiyStep2or7(leftmatrix, Dd, bundle1, axon, time,0);
       
       % Step (iii) of Reutskiy algorithm
       
      
       % Compute the interpolated vectors of the potentials. 
       for axnum = 1:bundle1.N
       axon(axnum).Vcoeff(:,time) = 0.5.*(axon(axnum).V(:,time+1) + axon(axnum).V(:,time)); 
       end
       

    
       % Step (iv) of Reutskiy algorithm
        for i = 1:bundle1.endlength
        for j=1:bundle1.N
            
            % execute the Ranvier nodal code
            %  Update the concentrations m, n, h, p 
      
            [mnew, nnew, hnew, pnew] = activationupdate(axon(j).m(i,t), axon(j).n(i,t), axon(j).h(i,t), axon(j).p(i,t), (axon(j).Vcoeff(i,t)), Q, bundle1.delT);

            axon(j).m(i,t+1) = mnew;
            axon(j).n(i,t+1) = nnew;
            axon(j).h(i,t+1) = hnew;
            axon(j).p(i,t+1) = pnew;
            end % Finish with Ranvier nodes
       % Skip myelin nodes
        end
        
     
    
        % Step (v) of Reutskiy algorithm
      
        
       
       % Compute the interpolated vectors of the activation variables. 
       for axnum = 1:bundle1.N
       axon(axnum).mcoeff(:,time) = 0.5.*(axon(axnum).m(:,time+1) + axon(axnum).m(:,time)); 
       axon(axnum).ncoeff(:,time) = 0.5.*(axon(axnum).n(:,time+1) + axon(axnum).n(:,time)); 
       axon(axnum).hcoeff(:,time) = 0.5.*(axon(axnum).h(:,time+1) + axon(axnum).h(:,time)); 
       axon(axnum).pcoeff(:,time) = 0.5.*(axon(axnum).p(:,time+1) + axon(axnum).p(:,time)); 

       end
        
       
       % Step (vi) of Reutskiy algorithm
       
       
       
         for i = 1:bundle1.endlength
        for j=1:bundle1.N
           
       % Update the ionic-coeff current
           
            axon(j).Ecoeff(i,t) = axon(j).Vcoeff(i,t) + E_rest;
            axon(j).Zfactorcoeff(i,t)=exp(-gamma*axon(j).Ecoeff(i,t));
            axon(j).Z_Nacoeff(i,t) = gamma*bundle1.F*axon(j).Ecoeff(i,t)*((bundle1.Nao*axon(j).Zfactorcoeff(i,t) - bundle1.Nai)/(1-axon(j).Zfactorcoeff(i,t)));
            
            axon(j).Z_Kcoeff(i,t) = gamma*bundle1.F*axon(j).Ecoeff(i,t)*((bundle1.Ko*axon(j).Zfactorcoeff(i,t) - bundle1.Ki)/(1-axon(j).Zfactorcoeff(i,t)));
            
            axnum = j;
            if ismember(i, bundle1.Ranvier)
            [axon(axnum).J_ioncoeff(i,time), axon(axnum).J_Nacoeff(i,time), axon(axnum).J_Kcoeff(i,time), axon(axnum).J_pcoeff(i,time), axon(axnum).J_lcoeff(i,time)] = ionicupdatetheta(axon(axnum).mcoeff(i,time),axon(axnum).ncoeff(i,time),axon(axnum).hcoeff(i,time),axon(axnum).pcoeff(i,time),axon(axnum).Z_Kcoeff(i,time),axon(axnum).Z_Nacoeff(i,time),V_l,(axon(axnum).Ecoeff(i,time)),bundle1.d, E_rest, theta);
% Include the current injections from fields.

% if axnum==1
%     otheraxons  = 2;
% else 
%     otheraxons = 1; 
% end
% cycle through every other node on axons other than the present one. 
% from each such node, gather its field at a distance , specifically, at
% the distance which intervenes between that node and the present node. 
if gain > 0
for maxon = setdiff(1:N, [axnum])
   for node = i %[i-20,i,i+20] %ignoring some `other' nodes, for saving time
  if (node>1)&&(node<300)
%       disp('hi5')
            if t ==1
          
            result=mainporeAug172020v2(axon(maxon).m(node,t),axon(maxon).n(node,t),axon(maxon).h(node,t),axon(maxon).p(node,t),t,0,0,0,0);
            oldresult = result; 
            fieldzero = 1e-40;
            myfield(maxon,node) = fieldzero;
            slopek = 1; %slope of K i-v relationship in mA/(mV times cm^2)
memthick = 1; % membrane thickness in meters
%gain = 1;
mynetcurrent =  gain.*myfield(maxon, node);%*slopek*memthick; assuming conductivity of sea water is 5S/m
        % size(axon(axnum).J_ioncoeff(i,time))
        if axon(axnum).J_ioncoeff(i,time) >0 %inward
 axon(axnum).J_ioncoeff(i,time)=axon(axnum).J_ioncoeff(i,time)+mynetcurrent; 
        elseif axon(axnum).J_ioncoeff(i,time)<=0
            axon(axnum).J_ioncoeff(i,time) = axon(axnum).J_ioncoeff(i,time)-mynetcurrent;
        end
            else
            result=mainporeAug172020(axon(maxon).m(node,t),axon(maxon).n(node,t),axon(maxon).h(node,t),axon(maxon).p(node,t),t,oldresult,axon(maxon).m(node,t-1),axon(maxon).n(node,t-1),axon(maxon).h(node,t-1),axon(maxon).p(node,t-1));
            %[c,h] = contourf(real(result{1,1}), real(result{1,2}),real(result{1,3}), 20);
%             oldsourcelowv = result{1,4}; % this is the open channel configuration used 
%             % to compute the result i.e., the field shown. in the next
%             % iteration, this will be merely modified slightly to save on
%             % computational time. 
           oldresult = result;
%             colorbar
%t1=clabel(c,h); % acquire the text and line objects created to label the plot. From these extract the weakest field value
myfield(maxon, node) = result{1,4}; %assuming t1 is 2 dimensional
%        caxis('auto');
%          %caxis([200 1000000])
% %caxis([0 realmax])
%        axis([2.5e-6 5.5e-6 0.1e-6 3.1e-6]);
%        xlabel('\rho (meters)');
% 
%        ylabel('z (meters)');
%        title(num2str(t));
%         drawnow
             
            
            slopek = 1; %slope of K i-v relationship in mA/(mV times cm^2)
memthick = 1; % membrane thickness in meters
%gain = 1;
mynetcurrent =  gain.*myfield(maxon, node);%*slopek*memthick; assuming conductivity of sea water is 5S/m
        % size(axon(axnum).J_ioncoeff(i,time))
        if axon(axnum).J_ioncoeff(i,time)>0 %inward .Refer to Figure 6 in Frankenhaueser-Huxley J. Physiol. (1964), 171, pp. 302-315
 axon(axnum).J_ioncoeff(i,time)=axon(axnum).J_ioncoeff(i,time)+ mynetcurrent; 
        elseif axon(axnum).J_ioncoeff(i,time) <= 0
            axon(axnum).J_ioncoeff(i,time) = axon(axnum).J_ioncoeff(i,time) - mynetcurrent;
        end
 end
  end
   end
end
end
            end
        end
         end

         
        
         
         
                
            % Create the right vector
            
    
             Ddcoeff = rightcreateReutskiy6v3(axon, bundle1, dee, dee_nd, t);
             
             
             % Step (vii) of Reutskiy algorithm
             axon = ReutskiyStep2or7(leftmatrix, Ddcoeff, bundle1, axon, time,1);
             
             
       
       % step (viii) of Reutskiy algorithm
        % Compute the interpolated vectors of the potentials. 
       for axnum = 1:bundle1.N
       axon(axnum).Vcoeff(:,time) = 0.5.*(axon(axnum).V(:,time+1) + axon(axnum).V(:,time)); 
       end
       
       % Step (ix) of Reutskiy algorithm
       
       for i = 1:bundle1.endlength
        for j=1:bundle1.N
            if ismember(i,bundle1.Ranvier)
            % execute the Ranvier nodal code
            %  Update the concentrations m, n, h, p 
            [axon(j).m(i,t+1), axon(j).n(i,t+1), axon(j).h(i,t+1), axon(j).p(i,t+1)] = activationupdate(axon(j).m(i,t), axon(j).n(i,t), axon(j).h(i,t), axon(j).p(i,t), (axon(j).Vcoeff(i,t)), Q, bundle1.delT);
            end % Finish with Ranvier nodes
    end    % Skip myelin nodes
       end
       
      % plottimes = [30, 100, 200, 300, 400, 500, 600, 800, 1000];
      
       
%             figure(1);
%             nodetoplot = 30;
% %             subplot(4,2,1)
% %             xlim([1,LTminus1])
% %             plot(1:t, axon(j).m(20, 1:t),'b*')
% %             xlabel('Time in 0.5 microsecond steps');
% %             ylabel('m(t)');
% %             drawnow
% %             subplot(4,2,2)
% %             plot(1:t, axon(j).n(20, 1:t), 'r*')
% %             xlabel('Time in 0.5 microsecond steps');
% %             ylabel('n(t)');
% %             drawnow
% %             subplot(4,2,3)
% %             plot(1:t, axon(j).h(20,1:t), 'g*');
% %             xlabel('Time in 0.5 microsecond steps');
% %             ylabel('h(t)');
% %             drawnow
% %             subplot(4,2,4)
% %             plot(1:t, axon(j).p(20,1:t), 'k*');
% %              xlabel('Time in 0.5 microsecond steps');
% %             ylabel('p(t)');
% %             drawnow
%             subplot(2,4,1)
%             plot(1:t, axon(1).V(nodetoplot,1:t), 'r*-');
%              xlabel('Time in 0.5 microsecond steps');
%             ylabel('V(t)');
%             drawnow
%             subplot(2,4,2)
%             plot(1:t, axon(1).m(nodetoplot,1:t), 'b*-');
%              xlabel('Time in 0.5 microsecond steps');
%             ylabel('m(t)');
%             drawnow
%             subplot(2,4,3)
%             plot(1:t, axon(1).h(nodetoplot,1:t), 'k*-');
%              xlabel('Time in 0.5 microsecond steps');
%             ylabel('h(t)');
%             drawnow
%             subplot(2,4,4)
%             plot(1:t, axon(1).n(nodetoplot,1:t), 'c*-');
%              xlabel('Time in 0.5 microsecond steps');
%             ylabel('n(t)');
%             title(strcat('You are looking at node ', num2str(nodetoplot)));
%             drawnow
%               subplot(2,4,5)
%             plot(1:t, axon(2).V(nodetoplot+10,1:t), 'r*-');
%              xlabel('Time in 0.5 microsecond steps');
%             ylabel('V2(t)');
%             drawnow
%             subplot(2,4,6)
%             plot(1:t, axon(2).m(nodetoplot+10,1:t), 'b*-');
%              xlabel('Time in 0.5 microsecond steps');
%             ylabel('m2(t)');
%             drawnow
%             subplot(2,4,7)
%             plot(1:t, axon(2).h(nodetoplot+10,1:t), 'k*-');
%              xlabel('Time in 0.5 microsecond steps');
%             ylabel('h2(t)');
%             drawnow
%             subplot(2,4,8)
%             plot(1:t, axon(2).n(nodetoplot+10,1:t), 'c*-');
%              xlabel('Time in 0.5 microsecond steps');
%             ylabel('n2(t)');
%             title(strcat('You are looking at node',num2str(nodetoplot+10)));
%             drawnow
%             subplot(4,2,6)
%                plot(1:t, axon(j).J_Na(20,1:t), 'c*');
%              xlabel('Time in 0.5 microsecond steps');
%             ylabel('Sodium Current J_{Na(t)}');
%             drawnow
%             subplot(4,2,7)
%             if t ==1
%           
%             result=mainporeAug172020v2(axon(j).m(20,t),axon(j).n(20,t),axon(j).h(20,t),axon(j).p(20,t),t,0,0,0,0)
%             oldresult = result; 
%             fieldzero = 1e-40;
%             myfield(1,1) = fieldzero;
%             else
%             result=mainporeAug172020(axon(j).m(20,t),axon(j).n(20,t),axon(j).h(20,t),axon(j).p(20,t),t,oldresult,axon(j).m(20,t-1),axon(j).n(20,t-1),axon(j).h(20,t-1),axon(j).p(20,t-1))
%             [c,h] = contourf(real(result{1,1}), real(result{1,2}),real(result{1,3}), 20);
% %             oldsourcelowv = result{1,4}; % this is the open channel configuration used 
% %             % to compute the result i.e., the field shown. in the next
% %             % iteration, this will be merely modified slightly to save on
% %             % computational time. 
%            oldresult = result;
% %             colorbar
% t1=clabel(c,h); % acquire the text and line objects created to label the plot. From these extract the weakest field value
% myfield(1,t) = min(min(t1)); %assuming t1 is 2 dimensional
% %        caxis('auto');
% %          %caxis([200 1000000])
% % %caxis([0 realmax])
% %        axis([2.5e-6 5.5e-6 0.1e-6 3.1e-6]);
% %        xlabel('\rho (meters)');
% % 
% %        ylabel('z (meters)');
% %        title(num2str(t));
% %         drawnow
   end

    outputLines{end+1} = 'Processing completed.';
    outputText = strjoin(outputLines, '\n');  % Join all lines into a single string with newlines
            
%    figure(10)
% for i = 1:endlength
%     for plotnumber = 1:N
%         subplot(N,1,plotnumber)
%     if ismember(i,Ranvier_array)
%        
% plot([1:LTminus1].*delT, axon(plotnumber).V(i, 1:LTminus1))
% axis([[1, LTminus1]*delT,-.1,.15]);
% hold on;
% grid on;
%     end
%     end
% end
% 
% title('Voltage vs Time at 1st six Ranvier nodes of axon 1 and 2');
% 
% grid on;
      
   end