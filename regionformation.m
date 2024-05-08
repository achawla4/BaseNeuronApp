function [A, V] = regionformation(wmatrix,alpha,beta,N)

%Wmatrix ==> containing the unadultered raw distances
%alpha,beta --> alpha is region A radius and beta is region V radius
%N --> W matrix size

l = 1;
wrowA = wmatrix(l,:);
[I,A, distA] = find(wrowA<alpha)

% A now contains axon numbers that will be called auditory

% pick an axon not in A
m=0;
while m==0
    for searchcount = 1:N
        if ~ismember(searchcount,A)
            m=searchcount;
        end
    end
end
disp(m)
% m is not in A

wrowV = wmatrix(m,:);

[I2, V,distV]= find(wrowV<beta)

C = intersect(A,V)
A = setdiff(A, C);
V = setdiff(V, C);
[]

