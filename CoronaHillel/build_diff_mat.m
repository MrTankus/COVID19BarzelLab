function A = build_diff_mat(locs,M_imm)
D = pdist(locs); % distances between all pairs
D = squareform(D); % to matrix
N = size(D,1);
f= @(d) d.^-2;
A = f(D);
A(1:N+1:end) = 0;
A = A./sum(A).*M_imm'; % sum of row i is # immigrants from the city

% sum of column j is not known: problem.
