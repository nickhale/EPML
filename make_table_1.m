% Generate Table 1 from "Prandtl’s extended mixing length model applied to the
% two-dimensional turbulent classical far wake" by Hutchinson et al.

bet = .01;

load('wygdata');
X = data(:,1); Y = data(:,2);
X(X<0,1) = -X(X<0,1); % Mirror x data

wyg = @(x) exp(-0.637*x.^2-0.056*x.^4);
fprintf('wyg: \t\t\t %f\n', norm(wyg(X) - Y, 2)) 

cev = @(xc) exp(-xc.^2*log(2));
fprintf('CEV: \t\t\t %f\n', norm(cev(X) - Y, 2)) 

[x, y, yfun] = EPML(0,bet,do_norm, do_iter);
fprintf('PML  (%2.2f): \t\t %f\n', bet, norm(yfun(X) - Y, 2))

for K2 = [0.1:.1:.5 .375]
    [x, y, yfun] = EPML(K2,bet,do_norm, do_iter);
    fprintf('EPML (%2.2f,%2.2f): \t %f\n', K2, bet, norm(yfun(X) - Y, 2))
end