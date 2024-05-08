load avtogether1.mat
axonav = axon
clear axon
load aalone1.mat
axona = axon
clear axon
load valone1.mat
axonv = axon
clear axon

N=30;
endlength = 250;
timesteps = 500;
LTminus1 = timesteps;
nodalspacing = 20;
Ranvier_array = 20:nodalspacing:endlength-20;
delT = 2e-6;


figure(2)

    for plotnumber = 1:N
        for i = 1:endlength
         subplot(N,1,plotnumber)
                outputs = screation(1);

    if ismember(i,Ranvier_array)
        datamya = axona(plotnumber).V(i,1:LTminus1);
         datamyv = axonv(plotnumber).V(i,1:LTminus1);
          datamyav = axonav(plotnumber).V(i,1:LTminus1);
%         spikeposonlya  = find(datamya>0.1);
%         spikevala = datamya(spikeposonlya);
%         yvalsa = embed(spikeposonlya,zeros(1,LTminus1),spikevala);
        
plot([1:LTminus1].*delT, datamyav-(datamya+datamyv), outputs);
axis([[1, LTminus1]*delT,-.2,.2]);
hold on;
grid on;
    end
    end
    end
