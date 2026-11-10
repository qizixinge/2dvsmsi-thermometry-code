function T2 = rev_t2fit(S)
% RDNLS grid fit of RDNLS.m (0:0.05:400 ms, n = 1..3, 8 ms spacing) without
% the edge-band magnitude mask, for voxel lists. S: nvox x 3 magnitudes.
t2srch = reshape(0:0.05:400,1,[]);
n = (1:3)';
E = exp(-n*8./t2srch);                 % 3 x D
den = sum(E.^2,1);                     % 1 x D
T2 = zeros(size(S,1),1);
chunk = 2000;
for i0 = 1:chunk:size(S,1)
    idx = i0:min(i0+chunk-1,size(S,1));
    s = S(idx,:);                      % m x 3
    J = sum(s.^2,2) - (s*E).^2./den;   % m x D
    [~,id] = min(J,[],2);
    T2(idx) = t2srch(id);
end
end
