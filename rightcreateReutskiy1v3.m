function Dd = rightcreateReutskiy1v3(axon, bundle1, dee, dee_nd, time)
egm = exp(-bundle1.gm*bundle1.delT/bundle1.cm);
egnd = exp(-bundle1.gnd*bundle1.delT/bundle1.cnd);
% egnd2 = exp(-bundle1.gnd*bundle1.delT/cnd2);

a = bundle1.a;
a_nd = a;
e_nd = bundle1.ee; 
ee = bundle1.ee;

             for i=1:bundle1.endlength
                % first work on the nodal positions
                if ismember(i,bundle1.Ranvier)
                    sumdiff1 = 0;
                    for n=1:bundle1.N
                        sumdiff1 = sumdiff1 + axon(n).V(i+1,time) - 2*axon(n).V(i,time) + axon(n).V(i-1,time);
                    end
                    
                    for l=1:bundle1.N
                        Dd((i-1)*bundle1.N+l,1) = egnd*(-a_nd*axon(l).V(i-1, time) + dee_nd*axon(l).V(i,time) - a_nd*axon(l).V(i+1,time) - e_nd*(sumdiff1-(axon(l).V(i+1,time) - 2*axon(l).V(i,time) + axon(l).V(i-1,time))))+ (-axon(l).J_ex(i,time)+axon(l).J_ion(i,time) - egnd*((axon(l).J_ex(i,time) - axon(l).J_ion(i,time))));
                    end
%                     % type2node positions
%                 elseif ismember(i,[40,50])
%                     sumdiff3 = 0;
%                     for n=1:bundle1.N
%                         sumdiff3 = sumdiff3 + axon(n).V(i+1,time) - 2*axon(n).V(i,time) + axon(n).V(i-1,time);
%                     end
%                     
%                     for l=1:bundle1.N
%                         Dd((i-1)*bundle1.N+l,1) = egnd2*(-a_nd*axon(l).V(i-1, time) + dee_nd2*axon(l).V(i,time) - a_nd*axon(l).V(i+1,time) - e_nd*(sumdiff1-(axon(l).V(i+1,time) - 2*axon(l).V(i,time) + axon(l).V(i-1,time))));
%                     end
                else % internodal positions
                    sumdiff2 = 0;
                    if i == 1 % implement the left b.c. (symmetric) on all axons
                        for n=1:bundle1.N
                        sumdiff2 = sumdiff2 + axon(n).V(i+1,time) - 1*axon(n).V(i,time);
                    end
                    for l=1:bundle1.N
                        Dd((i-1)*bundle1.N+l,1) = egm*((-a + dee)*axon(l).V(i,time) - a*axon(l).V(i+1,time) - ee*(sumdiff2-(axon(l).V(i+1,time) - axon(l).V(i,time) )));
                       
                    end
                      sumdiff3 = 0;
                    elseif i == bundle1.endlength % implement the right (zero) b.c. on all axons
                        for n=1:bundle1.N
                        sumdiff3 = sumdiff3 - 2*axon(n).V(i,time) + axon(n).V(i-1,time);
                    end
                    for l=1:bundle1.N
                        Dd((i-1)*bundle1.N+l,1) = egm*(-a*axon(l).V(i-1, time) + dee*axon(l).V(i,time) - ee*(sumdiff3-(- 2*axon(l).V(i,time) + axon(l).V(i-1,time))));
                    end
                    else
                        sumdiff = 0;
                        % implement the intermediate positions
                    for n=1:bundle1.N
                        sumdiff = sumdiff + axon(n).V(i+1,time) - 2*axon(n).V(i,time) + axon(n).V(i-1,time);
                    end
                    for l=1:bundle1.N
                        Dd((i-1)*bundle1.N+l,1) =  egm*(-a*axon(l).V(i-1, time) + dee*axon(l).V(i,time) - a*axon(l).V(i+1,time) - ee*(sumdiff-(axon(l).V(i+1,time) - 2*axon(l).V(i,time) + axon(l).V(i-1,time))));
                    end
                   
                    end
                end
             end 

             

