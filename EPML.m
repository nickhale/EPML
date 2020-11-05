function [x, y, yfun] = EPML(K2, bet, do_norm, do_iter)
%EPML  Solve Prandtl's extended mixing length model with Hermite spectral method 
%   [X, Y] = EPML(K2, B) solve the EPML model with $K_2$ = K2 and $\beta$ = B. 
%   Solutions are returned on a discrete grid of scaled Gauss-Hermite points.
%
%   The solution is 'normalised' so that y(0) = 1 and y(1) = 0.5. This can be 
%   disabled with ... = EPML(K2, B, DO_NORM) with DO_NORM = false.
%
%   By default, norm estimates in the Newton iterates are displayed to
%   check for convergence. This can be disabled with 
%   ... = EPML(K2, B, DO_NORM, DO_ITER) with DO_ITER = false.
%
%   [X, Y, YFUN] = EPML(K2, B) returns also a function handle of an interpolant 
%   to the solution, so that it can be evaluated at arbitrary values.
%
%   The code requires DMSUITE [1] and Chebfun [2] to be in the MATLAB path. 
%   [1] https://www.mathworks.com/matlabcentral/fileexchange/29-dmsuite
%   [2] https://www.chebfun.org

% For more details, see Prandtl’s extended mixing length model applied to the
% two-dimensional turbulent classical far wake" by Hutchinson et al.

% Nick Hale, Nov 2020

% Code summary: The algorithm discretised the DE using a Hermite spectral
% method (DMSUITE). the integral constraint is enforced by replacing a row
% of the Hermite collocation matrix with suitably scaled Guass-Hermite
% quadrature weights (Chebfun). The resulting nonlinear eqations are
% solving via a Newton iteration. A spline interpolant of the result is
% returned, since the global evaluation problem on the Hermite nodes can be
% ill-conditioned when the Hermite scaling, b, is large.

if ( nargin == 0 )
    % Model parameters. These can be changed:
    bet = 0.01;
    K2 = 0;
end
if ( nargin < 3 )
    do_norm = true;  % Normalise solutions so that y(0) = 1, y(1) = 0.5.
end
if ( nargin < 4 )
    do_iter = true;  % Show convergence of Newton iterates
end

% Discretisation parameters
n = 321;             % Discretisation size - can be changed - must be odd
b = 8.35;            % Discretisation scaling - heuristic - change carefully

% Construct differential operators (Hermite diffmats from DMSUITE):
[x, DD] = herdif(n,2,b); D1 = DD(:,:,1); D2 = DD(:,:,2);
zeroidx = (n+1)/2;   % Index of at x = 0 point
% DMSUITE doesn't ensure middle node is at zero.
x(zeroidx) = 0; X = diag(x);
% Quadrature weights (from Chebfun): 
[~,w] = hermpts(n, 'rec'); w = w.*exp((b*x).^2).'/b; 

% Initial guess:
y0 = exp(-0.637*x.^2-0.056*x.^4); % (from Wyg. et al data fit)
% y0 = exp(-log(2)*x.^2); % (from CEV model)
y = 2*y0/(w*y0);     % Normalise to satisfy integral constraint

% SOLVE HERE THE PROBLEM:
if ( do_iter ), disp('Convergence:'), end
for k = 1:10
    
    % Construct operators for current y:
    Dy = D1*y; DY = diag(Dy);   D2y = D2*y;       D2Y = diag(D2y);
    s2 = Dy.^2+K2^2*D2y.^2;     S2 = diag(s2);    DS2 = (DY*D1+K2^2*D2Y*D2);
    s = sqrt(s2);               S = diag(s);

    % Forward operator and Jacobian (multiplied through by s):
    Ny = s.*(x.*y + bet*Dy) + (s2.*Dy);
    Jy = S*(X + bet*D1) + (S2*D1 + DY*DS2);
   
    % Boundary conditions:
    Jy(zeroidx,:) = w; 
    Ny(zeroidx) = w*y - 2; % w*y = 2 since we integrate over [-inf,inf]
    
    % Newton solve:
    dy = -Jy\Ny; 
    y = y + dy;

    % Display convergence and terminate:
    if ( do_iter ), disp(norm(dy, inf)); end
    if ( norm(dy, inf)/max(y) < 1e-10 ), break, end
end

% Normalise/scale solution so that y(0) = 1 and y(1) = 0.5
if ( do_norm ) 
    y = y/y(zeroidx);
    yfun = @(z) spline(x, y, z);
    [~, xidx] = min(abs(y - 0.5)+100*(x<0));
    r = fzero(@(z) yfun(z) - 0.5, x(xidx+[-1 1]));
    if ( isempty(r) ), r = 1; end
    x = x/r;
else
    r = 1;
end

% Spline interpolant for finer plotting:
yfun = @(z) spline(x*r, y, r*z);
% (Evaluation of Hermite-function interpolant is unstable for large x if b >> 1)

end

