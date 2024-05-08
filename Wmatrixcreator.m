function wmatrix = Wmatrixcreator(N,M)
% This file creates a random Wmatrix of N rows and N cols
% These N axons are placed randomly in a grid of size M by M


% M = 100;
xp = randperm(M);
yp = randperm(M);
% N = 10;
is = xp(1:N);
js = yp(1:N);
kp = randperm(N);
S = sparse(is,js,kp)
wmatrix = interaxonal(S,N);
% Those are actual distances in wmatrix
% However, the larger the distance, the smaller the contribution of that
% axon to the present axon
% Thus we must invert each entry in wmatrix, after adding 1 (keeping in
% mind that the self-distance is 0)
wmatrix = 1+wmatrix;
wmatrix = 1./wmatrix
% normalize the wmatrix, i.e. largest distance must be 1.