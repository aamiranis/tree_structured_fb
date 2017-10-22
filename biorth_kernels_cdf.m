function [h0,g0] = biorth_kernels_cdf(K1,K2)

% K1, K2 must be even for now (will relax this later)

coeffs = poly([-1*ones(1,K1) repmat([-1+1i -1-1i],[1 K2])]);

K = length(coeffs)-1;
temp = [zeros(1,K-1) fliplr(coeffs) zeros(1,K-1)];

C = zeros(2*K,K);
for i = 1:2*K
    C(i,:) = temp( (3*K-1)-(K-1)-(i-1) : (3*K-1)-(i-1) );
end
A = C(2:2:end,:); % matlab ordering of coeffs is reverse of the usual
b = zeros(K,1);
b(end) = 1; % matlab ordering of coeffs is reverse of the usual

r = (A \ b)';
p = conv(coeffs,r);

% % Spectral factorization - maximal balancing and orthogonality

roots_r = roots(r);
num_roots_r = length(roots_r);
all_perms = nchoosek(1:num_roots_r, floor(num_roots_r/2));
valid_perms = [];
for i = 1:size(all_perms,1)
    indices = all_perms(i,:);
    if (abs(imag(sum(roots_r(indices)))) < eps)
        valid_perms = [valid_perms; indices]; 
    end
end

num_r = 10;
num_theta = 18;
best_perm = 1;
max_orth = 0;
for i = 1:size(valid_perms,1)
    indices = valid_perms(i,:);
    indices_comp = setdiff(1:num_roots_r, indices);
    h0 = poly([ -1*ones(1,K1/2)...
                repmat([-1+1i -1-1i],[1 K2/2])...
                roots_r(indices)' ]);
    g0 = poly([ -1*ones(1,K1/2)...
                repmat([-1+1i -1-1i],[1 K2/2])...
                roots_r(indices_comp)' ]);
    a = conv(h0,g0);
    factor = p(1)/a(1);
    h0 = sqrt(abs(factor)) * h0;
    g0 = sign(factor) * sqrt(abs(factor)) * g0;
    h0 = sign(polyval(h0,1)) * h0;
    g0 = sign(polyval(g0,1)) * g0;

    % evaluate orthogonality
    Lam = (1/num_r:1/num_r:1)'*exp(1i*2*pi*(0:num_theta-1)/num_theta);
    Lam = Lam(:);
    H1 = polyval(g0,-Lam);
    H0 = polyval(h0,Lam);
    
%     Pow = (H1.*conj(H1) + H0.*conj(H0)).^(0.5);
    Pow = abs(H1.*H1 + H0.*H0).^(0.5);
    % approx. orthogonality measure from Narang & Ortega '13
    % range: [0 1], 1 -> most orthogonal, 0 -> least orthogonal
    orth = 1 - abs(max(Pow) - min(Pow))/abs(max(Pow) + min(Pow));
    
    if (orth > max_orth)
        best_perm = i;
        max_orth = orth;
%         PR = H1.*fliplr(H0) + H0.*fliplr(H1);
    end
end

indices = valid_perms(best_perm,:);
indices_comp = setdiff(1:num_roots_r, indices);
h0 = poly([ -1*ones(1,K1/2)...
            repmat([-1+1i -1-1i],[1 K2/2])...
            roots_r(indices)' ]);
g0 = poly([ -1*ones(1,K1/2)...
            repmat([-1+1i -1-1i],[1 K2/2])...
            roots_r(indices_comp)' ]);
a = conv(h0,g0);
factor = p(1)/a(1);
h0 = sqrt(abs(factor)) * h0;
g0 = sign(factor) * sqrt(abs(factor)) * g0;
h0 = sign(polyval(h0,1)) * h0;
g0 = sign(polyval(g0,1)) * g0;