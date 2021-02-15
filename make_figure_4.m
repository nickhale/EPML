% Generate Figure 4 from "Prandtl’s extended mixing length model applied to the
% two-dimensional turbulent classical far wake" by Hutchinson et al.

% Plot solutions of various closure models alongside the extracted data..

% Nick Hale, Nov 2020 (Revised Jan 2021)

% Choose beta:
bet = 0.01; 
% Marker size in plotting
MS = 7; 

% Set options for EPML code:
do_norm = true; 
do_iter = true;

% We're going to use a slightly varied colour order:
cols = get(gca, 'colororder');

% Plot the data fit:
xx = linspace(0, 2, 1001);
wyg = @(x) exp(-0.637*x.^2-0.056*x.^4);
h2 = plot(xx, wyg(xx), '-k', 'LineWidth', 3); hold on

% Solve for various K_2 values:
k = 0;
for K2 = [0 .25 .5]
    k = k + 1;
    [x, y, yfun] = EPML(K2, bet, do_norm, do_iter);
    h456{k} = plot(x, y, 'color', cols(k,:), 'LineWidth', 3);
end

% Plot the CEV model solution:
cev = @(x) exp(-x.^2*log(2));
h3 = plot(xx, cev(xx), ':', 'color', cols(4,:), 'LineWidth', 3);

% Load and plot the extracted data:
load('wygdata');
X = data(:,1); X(X<0,1) = -X(X<0,1); % Mirror x data
Y = data(:,2); 
h1 = plot(X, Y, 'ok', ...
    'markersize', MS, 'markerfacecolor', 'c', 'linewidth', 1.5); 

axis([0 2 0 1]), shg

% Legend and axes labels:
legend([h1,h2,h3, h456{:}], ...
    'Data', '(5.1)', 'CEV', 'PML', 'EPML ($\tilde K_2 = 0.25$)', ...
    'EPML ($\tilde K_2 = 0.5$)', 'interpreter', 'latex'), shg
xlabel('$|\xi_N|$', 'interpreter', 'latex')
ylabel('$F_N(\xi_N)$', 'interpreter', 'latex')
set(gca, 'fontsize', 16)

% Make it nice and big:
set(gcf,  'position', [235        1081        1450        1081]);

% Print:
print -depsc figure4.eps
print -dpng figure4.png

