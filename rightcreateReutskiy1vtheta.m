function Dd = rightcreateReutskiy1vtheta(axon, bundle1,time,theta)
egm = 1; %exp(-bundle1.gm*bundle1.delT/bundle1.cm);
egnd =1;% exp(-bundle1.gnd*bundle1.delT/bundle1.cnd);
K = bundle1.K;
a = bundle1.a; % a = (alpha_N-1)/(2*bundle1.rf*(bundle1.delX)^2);
a_nd = a;
e_nd = bundle1.ee; 
ee = bundle1.ee; % ee = alpha_N/(2*bundle1.rf*(bundle1.delX)^2);
           

             for i=1:bundle1.endlength
                % first work on the nodal positions
                if ismember(i,bundle1.Ranvier)
                    for l=1:bundle1.N
                    sumdiff1(l) = 0;
                    
                    % Repair sumdiff1
                    for n=1:bundle1.N
                        
                        
                        sumdiff1(l) = sumdiff1(l) + K(l,n)*(axon(n).V(i+1,time) - 2*axon(n).V(i,time) + axon(n).V(i-1,time));
                    end
                    bundle1.rf(1,l) = bundle1.rf/theta(1,l);
            alpha_0(1,l) = bundle1.r0/bundle1.rf(1,l) ; 
alpha_N(1,l) = alpha_0(1,l)/(1 + bundle1.N*alpha_0(1,l));

 a_nd(1,l) = (alpha_N(1,l)-1)/(2*bundle1.rf(1,l)*(bundle1.delX)^2);
                   % Repair a_nd
                   a_nd(1,l) = (a_nd(1,l) + (1/(2*bundle1.rf(1,l)*(bundle1.delX)^2)))*K(l,l) - (1/(2*bundle1.rf(1,l)*(bundle1.delX)^2)); 
                   % Repair dee_nd
                   
            dee_nd(1,l) = bundle1.cnd/bundle1.delT + (alpha_N(1,l) - 1)/(bundle1.rf(1,l)*(bundle1.delX)^2);
                   dee_nd(1,l) = (dee_nd(1,l) - bundle1.cnd/bundle1.delT + (1/(bundle1.rf(1,l)*(bundle1.delX)^2)))*K(l,l) + bundle1.cnd/bundle1.delT - 1/(bundle1.rf(1,l)*(bundle1.delX)^2);
                   ee = alpha_N/(2*bundle1.rf*(bundle1.delX)^2);
                   
                        Dd((i-1)*bundle1.N+l,1) = egnd*(-a_nd(1,l)*axon(l).V(i-1, time) + dee_nd(1,l)*axon(l).V(i,time) - a_nd(1,l)*axon(l).V(i+1,time) - e_nd*(sumdiff1(l)-K(l,l)*(axon(l).V(i+1,time) - 2*axon(l).V(i,time) + axon(l).V(i-1,time))))+ (-axon(l).J_ex(i,time)+axon(l).J_ion(i,time) - egnd*((axon(l).J_ex(i,time) - axon(l).J_ion(i,time))));
                   
                    end
                    
                else % internodal positions
                    for l=1:bundle1.N
                    sumdiff2(l) = 0;
                    if i == 1 % implement the left b.c. (symmetric) on all axons
                        for n=1:bundle1.N
                        sumdiff2(l) = sumdiff2(l) + K(l,n)*(axon(n).V(i+1,time) - 1*axon(n).V(i,time));
                    end
                    % Repair dee and a 
                      dee = bundle1.cm/bundle1.delT - (1-alpha_N)/(bundle1.rf*((bundle1.delX)^2)) - (bundle1.gm/2);
            dee_nd = bundle1.cnd/bundle1.delT + (alpha_N - 1)/(bundle1.rf*(bundle1.delX)^2);
                    dee =(dee-bundle1.cm/bundle1.delT + (1/(bundle1.rf*(bundle1.delX)^2)) + bundle1.gm/2)*K(l,l) + bundle1.cm/bundle1.delT - (1/(bundle1.rf*(bundle1.delX)^2)) - bundle1.gm/2;
                   a = (a_nd + (1/(2*bundle1.rf*(bundle1.delX)^2)))*K(l,l) - (1/(2*bundle1.rf*(bundle1.delX)^2)); 
                   
                    Dd((i-1)*bundle1.N+l,1) = egm*((-a + dee)*axon(l).V(i,time) - a*axon(l).V(i+1,time) - ee*(sumdiff2(l)-K(l,l)*(axon(l).V(i+1,time) - axon(l).V(i,time) )));
                       
                    
                      sumdiff3(l) = 0;
                    elseif i == bundle1.endlength % implement the right (zero) b.c. on all axons
                        for n=1:bundle1.N
                        sumdiff3(l) = sumdiff3(l) +K(l,n)*(-2*axon(n).V(i,time) + axon(n).V(i-1,time));
                        end
                    % Repair dee and a 
                      dee = bundle1.cm/bundle1.delT - (1-alpha_N)/(bundle1.rf*((bundle1.delX)^2)) - (bundle1.gm/2);
            dee_nd = bundle1.cnd/bundle1.delT + (alpha_N - 1)/(bundle1.rf*(bundle1.delX)^2);
                    dee =(dee-bundle1.cm/bundle1.delT + (1/(bundle1.rf*(bundle1.delX)^2)) + bundle1.gm/2)*K(l,l) + bundle1.cm/bundle1.delT - (1/(bundle1.rf*(bundle1.delX)^2)) - bundle1.gm/2;
                   a = (a_nd + (1/(2*bundle1.rf*(bundle1.delX)^2)))*K(l,l) - (1/(2*bundle1.rf*(bundle1.delX)^2)); 
                
                    
                        Dd((i-1)*bundle1.N+l,1) = egm*(-a*axon(l).V(i-1, time) + dee*axon(l).V(i,time) - ee*(sumdiff3(l)-K(l,l)*(- 2*axon(l).V(i,time) + axon(l).V(i-1,time))));
                    
                    else
                        sumdiff(l) = 0;
                        % implement the intermediate positions
                    for n=1:bundle1.N
                        sumdiff(l) = sumdiff(l) + K(l,n)*(axon(n).V(i+1,time) - 2*axon(n).V(i,time) + axon(n).V(i-1,time));
                    end
                    
                    % Repair dee and a 
                      dee = bundle1.cm/bundle1.delT - (1-alpha_N)/(bundle1.rf*((bundle1.delX)^2)) - (bundle1.gm/2);
            dee_nd = bundle1.cnd/bundle1.delT + (alpha_N - 1)/(bundle1.rf*(bundle1.delX)^2);
                    dee =(dee-bundle1.cm/bundle1.delT + (1/(bundle1.rf*(bundle1.delX)^2)) + bundle1.gm/2)*K(l,l) + bundle1.cm/bundle1.delT - (1/(bundle1.rf*(bundle1.delX)^2)) - bundle1.gm/2;
                   a = (a_nd + (1/(2*bundle1.rf*(bundle1.delX)^2)))*K(l,l) - (1/(2*bundle1.rf*(bundle1.delX)^2)); 
                   
                    
                        Dd((i-1)*bundle1.N+l,1) =  egm*(-a*axon(l).V(i-1, time) + dee*axon(l).V(i,time) - a*axon(l).V(i+1,time) - ee*(sumdiff(l)-K(l,l)*(axon(l).V(i+1,time) - 2*axon(l).V(i,time) + axon(l).V(i-1,time))));
                    
                   
                    end
                    end
                end
             end

                
              

             

