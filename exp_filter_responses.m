% Filter responses using the CDF maximally-flat design
% Roots are placed at -1, -1+i, -1-i

[h0,g0] = biorth_kernels_cdf(2,2);

%% Kernels h0 and h1

[xx,yy] = meshgrid(-1:0.02:1,-1:0.02:1);
z = xx + 1i*yy;
xc = cos(2*pi*(0:1/50:1));
yc = sin(2*pi*(0:1/50:1));

% Plot h0
f = polyval(h0,z);
f(abs(z) > 1) = NaN + 1i*NaN;

figure;
contourf(xx,yy,abs(f),20);
xlabel('Re(\lambda)','FontSize',14);
ylabel('Im(\lambda)','FontSize',14);
set(gca,'YTick',[-1 -0.5 0 0.5 1]);
axis equal;
caxis([0 2]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
title('|h_0(\lambda)|','FontSize',14);
set(gca,'FontSize',14);
export_fig('plots/response_mag_h0_lambda.pdf','-transparent');

figure;
contourf(xx,yy,angle(f),20);
xlabel('Re(\lambda)','FontSize',14);
ylabel('Im(\lambda)','FontSize',14);
set(gca,'YTick',[-1 -0.5 0 0.5 1]);
axis equal;
caxis([-pi pi]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
title('\angle(h_0(\lambda)) (radians)','FontSize',14);
set(gca,'FontSize',14);
export_fig('plots/response_phase_h0_lambda.pdf','-transparent');


% Plot h1
f = polyval(g0,-z);
f(abs(z) > 1) = NaN + 1i*NaN;

figure;
contourf(xx,yy,abs(f),20);
xlabel('Re(\lambda)','FontSize',14);
ylabel('Im(\lambda)','FontSize',14);
set(gca,'YTick',[-1 -0.5 0 0.5 1]);
axis equal;
caxis([0 2]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
title('|h_1(\lambda)|','FontSize',14);
set(gca,'FontSize',14);
export_fig('plots/response_mag_h1_lambda.pdf','-transparent');

figure;
contourf(xx,yy,angle(f),20);
xlabel('Re(\lambda)','FontSize',14);
ylabel('Im(\lambda)','FontSize',14);
set(gca,'YTick',[-1 -0.5 0 0.5 1]);
axis equal;
caxis([-pi pi]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
title('\angle(h_1(\lambda)) (radians)','FontSize',14);
set(gca,'FontSize',14);
export_fig('plots/response_phase_h1_lambda.pdf','-transparent');

%% Responses of four channels of two-stage filterbank

[xx,yy] = meshgrid(-1:0.02:1,-1:0.02:1);
z = xx + 1i*yy;
xc = cos(2*pi*(0:1/50:1));
yc = sin(2*pi*(0:1/50:1));

% Plot h0(\lambda) * h0(\lambda^2)

f = polyval(h0,z) .* polyval(h0,z.*z);
f(abs(z) > 1) = NaN + 1i*NaN;

figure;
contourf(xx,yy,abs(f),20);
axis equal; axis off;
caxis([0 3]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
set(gca,'FontSize',16);
export_fig('plots/response_mag_ch0.pdf','-transparent');

figure;
contourf(xx,yy,angle(f),20);
axis equal; axis off;
caxis([-pi pi]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
set(gca,'FontSize',16);
export_fig('plots/response_phase_ch0.pdf','-transparent');


% Plot h0(\lambda) * h1(\lambda^2)

f = polyval(h0,z) .* polyval(g0,-z.*z);
f(abs(z) > 1) = NaN + 1i*NaN;

figure;
contourf(xx,yy,abs(f),20);
axis equal; axis off;
caxis([0 3]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
set(gca,'FontSize',16);
export_fig('plots/response_mag_ch1.pdf','-transparent');

figure;
contourf(xx,yy,angle(f),20);
axis equal; axis off;
caxis([-pi pi]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
set(gca,'FontSize',16);
export_fig('plots/response_phase_ch1.pdf','-transparent');


% Plot h1(\lambda) * h0(\lambda^2)

f = polyval(g0,-z) .* polyval(h0,z.*z);
f(abs(z) > 1) = NaN + 1i*NaN;

figure;
contourf(xx,yy,abs(f),20);
axis equal; axis off;
caxis([0 3]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
set(gca,'FontSize',16);
export_fig('plots/response_mag_ch2.pdf','-transparent');

figure;
contourf(xx,yy,angle(f),20);
axis equal; axis off;
caxis([-pi pi]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
set(gca,'FontSize',16);
export_fig('plots/response_phase_ch2.pdf','-transparent');


% Plot h1(\lambda) * h1(\lambda^2)

f = polyval(g0,-z) .* polyval(g0,-z.*z);
f(abs(z) > 1) = NaN + 1i*NaN;

figure;
contourf(xx,yy,abs(f),20);
axis equal; axis off;
caxis([0 3]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
set(gca,'FontSize',16);
export_fig('plots/response_mag_ch3.pdf','-transparent');

figure;
contourf(xx,yy,angle(f),20);
axis equal; axis off;
caxis([-pi pi]);
colorbar;
hold on;
plot(xc,yc,'k- ','LineWidth',2);
set(gca,'FontSize',16);
export_fig('plots/response_phase_ch3.pdf','-transparent');

close all;