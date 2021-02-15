% Generate Figure 5 from "Prandtl’s extended mixing length model applied to the
% two-dimensional turbulent classical far wake" by Hutchinson et al.

% Plot shear stresses of various closure models alongside the extracted data..

% Nick Hale, Jan 2021

close all

% Load the extracted data:
load('wygdata11');
% Airfoil
S1 = .103; X1 = data1(:,1); Y1 = data1(:,2); str = 'airfoil';
% % Solid strip
S2 = .072; X2 = data2(:,1); Y2 = data2(:,2); str = 'solid';

uv1 = @(xi,F) -S1*xi.*F;
uv2 = @(xi,F) -S2*xi.*F;

% Choose beta:
bet = 0.0; 
% Marker size in plotting
MS = 7; 

% Set options for EPML code:
do_norm = true; 
do_iter = true;

% We're going to use a slightly varied colour order:
cols = get(gca, 'colororder');

% Plot the data fit:
xx = linspace(-3, 3, 1001);
wyg = @(x) exp(-0.637*x.^2-0.056*x.^4);
xi = xx;
F = wyg(xx);
h3 = plot(xi, uv1(xi, F), 'k', 'LineWidth', 3); hold on
plot(xi, uv2(xi, F), 'k', 'LineWidth', 3)

% Solve for various K_2 values:
[xi, F, yfun] = EPML(0, bet, do_norm, do_iter);
h5 = plot(xi, uv1(xi,F), 'color', cols(1,:), 'LineWidth', 3);
plot(xi, uv2(xi,F), 'color', cols(1,:), 'LineWidth', 3)

K2 = 0.375;
[xi, F, yfun] = EPML(K2, bet, do_norm, do_iter);
h6 = plot(xi, uv1(xi,F), 'color', cols(3,:), 'LineWidth', 3);
plot(xi, uv2(xi,F), 'color', cols(3,:), 'LineWidth', 3)

% Plot the CEV model solution:
cev = @(x) exp(-x.^2*log(2));
xi = xx;
F = cev(xx);
h4 = plot(xi, uv1(xi, F), ':', 'color', cols(4,:), 'LineWidth', 3);
plot(xi, uv2(xi, F), ':', 'color', cols(4,:), 'LineWidth', 3)

% Plot the data:
h1 = plot(X1, Y1, 'ok', ...
    'markersize', MS, 'markerfacecolor', 'c', 'linewidth', 1.5); hold on
h2 = plot(X2, Y2, 'ok', ...
    'markersize', MS, 'markerfacecolor', 'm', 'linewidth', 1.5); hold on

% Legend and axes labels:
legend([h1 h2 h3 h4 h5 h6], 'Airfoil data', 'Solid strip data', '(5.1)', ...
    'CEV', 'PML', 'EPML ($\tilde K_2 = 0.375$)', ...
'interpreter', 'latex'), shg
xlabel('$\xi_N$', 'interpreter', 'latex')
ylabel('$g_N(\xi_N)$','interpreter', 'latex')
set(gca, 'fontsize', 16)

axis([-3 3 -.06 .06])

% Make it nice and big:
set(gcf,  'position', [235        1081        1450        1081]);

% Print:
print -depsc figure5.eps
print -dpng figure5.png

