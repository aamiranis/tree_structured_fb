rng('default');

addpath(genpath('exportfig'));

%% Graph definition

N_list = [25 25 25 25];
i = [0 cumsum(N_list)];
j = [cumsum(N_list) N_list(1)];

M = length(N_list);
N = sum(N_list);

A = zeros(N);
for k = 1:M
    N_curr = N_list(k);
    N_next = N_list(mod(k,M)+1);
%     P = rand(N_curr,N_next);
    P = eye(N_curr);
    Q = double(rand(N_curr,N_next) < 0.2);
    Q = Q - diag(diag(Q));
    P = P + Q;
    P = diag(sum(P,2).^-1)*P;
    A(i(k)+1:i(k+1),mod(j(k),N)+1:j(k+1)) = P;
end

Coordinates = zeros(N,2);
for k = 1:M
%     Coordinates(i(k)+1:i(k+1),1) = cos((k-1)/M*2*pi) + 0.75*(rand(N_list(k),1) - 0.5);
%     Coordinates(i(k)+1:i(k+1),2) = sin((k-1)/M*2*pi) + 0.75*(rand(N_list(k),1) - 0.5);
    % arrange points on a rectangular grid
    K = N_list(k);
    xx = zeros(K,1);
    yy = zeros(K,1);
    yy_val = 0;
    for n = 0:K-1
        xx_val = mod(n,floor(sqrt(K)));
        if (xx_val == 0), yy_val = yy_val + 1; end
        xx(n+1) = xx_val + 1;
        yy(n+1) = yy_val;
    end
    Coordinates(i(k)+1:i(k+1),1) = xx/floor(sqrt(K)) + 1.5*cos((k-1)/M*2*pi);
    Coordinates(i(k)+1:i(k+1),2) = yy/floor(sqrt(K)) + 1.5*sin((k-1)/M*2*pi);
end

[V,D] = eig(A);
lambda = diag(D);

% Plot graph
figure;
gplot(A,Coordinates,'k- ');
hold on;
scatter(Coordinates(:,1), Coordinates(:,2), 80, [1 0 0], 'Filled');
axis equal; axis off;
export_fig('plots/filtexp_4bc_graph.pdf','-transparent');

% Plot eigenvalues
figure;
viscircles([0 0], 1, 'LineStyle', '--');
hold on;
scatter(real(lambda), imag(lambda), 40, [0 0 0], 'filled');
axis equal; axis off;
export_fig('plots/filtexp_4bc_eig_values.pdf','-transparent');

%% Filtering experiment

f = zeros(N,1);
f(i(1)+1:i(2)) = rand(N_list(1),1);

% % Ideal filtering

h0_ideal = abs(angle(lambda)) <= pi/4;
h1_ideal = abs(angle(lambda)) > pi/4 & abs(angle(lambda)) <= pi/2;
h2_ideal = abs(angle(lambda)) >= 3*pi/4;
h3_ideal = abs(angle(lambda)) > pi/2 & abs(angle(lambda)) < 3*pi/4;

f0_ideal = V * (h0_ideal .* (V^-1 * f));
f1_ideal = V * (h1_ideal .* (V^-1 * f));
f2_ideal = V * (h2_ideal .* (V^-1 * f));
f3_ideal = V * (h3_ideal .* (V^-1 * f)); 

% % Polynomial filtering

[h0_poly,g0_poly] = biorth_kernels_cdf(2,2);
f0_poly = polyvalm(h0_poly, A) * polyvalm(h0_poly, A*A) * f;
f1_poly = polyvalm(h0_poly, A) * polyvalm(g0_poly,-A*A) * f;
f2_poly = polyvalm(g0_poly,-A) * polyvalm(h0_poly, A*A) * f;
f3_poly = polyvalm(g0_poly,-A) * polyvalm(g0_poly,-A*A) * f;

% % Plotting

%% Plot results

% Original signal
figure;
gplot(A,Coordinates,'k- ');
hold on;
axis equal; axis off;
scatter(Coordinates(:,1),Coordinates(:,2),80,real(f),'Filled'); colorbar;
caxis([-1 1]);
export_fig('plots/filtexp_4bc_signal.pdf','-transparent');

% Ideal filter responses
figure;
gplot(A,Coordinates,'k- ');
hold on;
axis equal; axis off;
scatter(Coordinates(:,1),Coordinates(:,2),80,real(f0_ideal),'Filled'); colorbar;
caxis([-1 1]);
export_fig('plots/filtexp_4bc_ideal_0.pdf','-transparent');

figure;
gplot(A,Coordinates,'k- ');
hold on;
axis equal; axis off;
scatter(Coordinates(:,1),Coordinates(:,2),80,real(f1_ideal),'Filled'); colorbar;
caxis([-1 1]);
export_fig('plots/filtexp_4bc_ideal_1.pdf','-transparent');

figure;
gplot(A,Coordinates,'k- ');
hold on;
axis equal; axis off;
scatter(Coordinates(:,1),Coordinates(:,2),80,real(f2_ideal),'Filled'); colorbar;
caxis([-1 1]);
export_fig('plots/filtexp_4bc_ideal_2.pdf','-transparent');

figure;
gplot(A,Coordinates,'k- ');
hold on;
axis equal; axis off;
scatter(Coordinates(:,1),Coordinates(:,2),80,real(f3_ideal),'Filled'); colorbar;
caxis([-1 1]);
export_fig('plots/filtexp_4bc_ideal_3.pdf','-transparent');

% Polynomial filter responses
figure;
gplot(A,Coordinates,'k- ');
hold on;
axis equal; axis off;
scatter(Coordinates(:,1),Coordinates(:,2),80,real(f0_poly),'Filled'); colorbar;
caxis([-1 1]);
export_fig('plots/filtexp_4bc_poly_0.pdf','-transparent');

figure;
gplot(A,Coordinates,'k- ');
hold on;
axis equal; axis off;
scatter(Coordinates(:,1),Coordinates(:,2),80,real(f1_poly),'Filled'); colorbar;
caxis([-1 1]);
export_fig('plots/filtexp_4bc_poly_1.pdf','-transparent');

figure;
gplot(A,Coordinates,'k- ');
hold on;
axis equal; axis off;
scatter(Coordinates(:,1),Coordinates(:,2),80,real(f2_poly),'Filled'); colorbar;
caxis([-1 1]);
export_fig('plots/filtexp_4bc_poly_2.pdf','-transparent');

figure;
gplot(A,Coordinates,'k- ');
hold on;
axis equal; axis off;
scatter(Coordinates(:,1),Coordinates(:,2),80,real(f3_poly),'Filled'); colorbar;
caxis([-1 1]);
export_fig('plots/filtexp_4bc_poly_3.pdf','-transparent');

close all;