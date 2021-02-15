% Generate Figure 3 from "Prandtl’s extended mixing length model applied to the
% two-dimensional turbulent classical far wake" by Hutchinson et al.

% Import a scanned image From the WYgnanski et al paper and overlay
% extracted data points.

% Nick Hale, Nov 2020

% Import image:
A = imread('wyg11.bmp');
[Am,An] = size(A);
C = ~A(:,1:An/3);
C = C(30:end-100,160:end-10);
[i,j] = find(flipud(C));
% Correct distortions:
ii = 6.125*(j/max(j)-.51)+.00005*i;
jj = .12*(1.01*i/max(i)-.000006*j)-0.06;
idx = abs(ii) > 2.95; ii(idx) = []; jj(idx) = [];
numel(ii)
plot(ii,jj, '.k'), axis on, hold on, shg

% Add the extracted data:
load('wygdata11')
X1 = data1(:,1); Y1 = data1(:,2);
plot(X1, Y1, 'ok', ...
    'markersize', 3, 'markerfacecolor', 'c', 'linewidth', 1)
h1 = plot(NaN, NaN, 'ok', 'markerfacecolor', 'c', 'linewidth', 1);
X2 = data2(:,1); Y2 = data2(:,2);
plot(X2, Y2, 'ok', ...
    'markersize', 3, 'markerfacecolor', 'm', 'linewidth', 1)
h2 = plot(NaN, NaN, 'ok', 'markerfacecolor', 'm', 'linewidth', 1);
hold off

% Add legend and axes labels:
l = legend([h1, h2], 'Airfoil data', 'Solid strip data', 'fontsize', 6);
xlabel('$\xi_N$','interpreter','latex'), ...
ylabel('$g_N(\xi_N)$','interpreter','latex')
axis([-3 3 -.0665 .0665])

set(gca, 'fontsize', 10);

% Save
set(gcf, 'Renderer', 'Painters');
print -depsc figure3b.eps
set(gcf, 'papersize', [12 9])
print(gcf, 'figure3b.pdf', '-dpdf', '-fillpage')
!cp figure3b.pdf ~/Desktop/paper




