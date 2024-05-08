function wmatrix = interaxonal(m,N)
% m is the input sparse matrix
% wmatrix is the output full matrix storing interaxonal distance values
initwmatrix = zeros(N,N);
% normalize m before use.
maxentry = max(m);
newm = m;
[is,js,vals]=find(newm);
myfull = full(newm);
for i = 1:N
    for j = 1:N
        % find the location of the ith axon
        iindex=find(vals==i);
        % find the location of the jth axon
        jindex = find(vals==j);
        % Compute the distance between ith and jth axons
        distance = sqrt((is(iindex)- is(jindex))^2+(js(iindex)-js(jindex))^2);
        % store distance as W(i,j)
        initwmatrix(i,j)=distance;
    end
end
wmatrix = initwmatrix./(sqrt(2)*N);