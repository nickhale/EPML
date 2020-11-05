% Generate Figure 3 from "Prandtl’s extended mixing length model applied to the
% two-dimensional turbulent classical far wake" by Hutchinson et al.

% Import a scanned image From the WYgnanski et al paper and overlay
% extracted data points.

% Nick Hale, Nov 2020

% Import image:
A = imread('wygnanski.bmp');
[Am,An] = size(A);
C = ~A(:,1:An/3);
C = C(20:end-100,200:end-175);
[i,j] = find(flipud(C));
% Correct distortions:
ii = 7.655*(j/max(j)-.51)+.000015*i;
jj = 1.0455*i/max(i)-.03;
plot(ii,jj, '.k'), axis on, hold on, shg

% Add the extracted data:
load('wygdata')
X = data(:,1); Y = data(:,2);
plot(X, Y, 'ok', ...
    'markersize', 3, 'markerfacecolor', 'c', 'linewidth', 1)
h = plot(NaN, NaN, 'ok', 'markerfacecolor', 'c', 'linewidth', 1);
hold off

% Add legend and axes labels:
legend(h, 'Extracted data')
xlabel('$\xi_N$','interpreter','latex'), ...
ylabel('$F_N(\xi_N)$','interpreter','latex')
axis([-4 4, -.05 1.05]) 

% Save
print -depsc figure3.eps
print -dpng figure3.png
