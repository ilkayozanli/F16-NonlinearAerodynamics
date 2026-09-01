%RUN_CONTINUATION  Pseudo-arclength continuation of F-16 equilibria.
%
%   The Kolb (2017) method applied to the F-16: track equilibria of the full
%   eight-state system f(x,u) = 0 as the stabilator is varied, classify each
%   by the eigenvalues of the Jacobian, and detect folds and Hopf points.
%
%   Why pseudo-arclength rather than a plain sweep.  A plain sweep solves
%   f(x,delta_h) = 0 independently at each delta_h.  It cannot go round a
%   fold, because at a fold the branch turns back on itself and there is no
%   solution for the next parameter value.  Pseudo-arclength instead treats
%   the parameter as another unknown and steps along the branch by arclength,
%   so a fold is just a point where the branch changes direction.  This is
%   what MATCONT does internally and what Kolb used.
%
%   How it works, per step:
%     1. Compute the tangent to the solution curve as the null vector of the
%        augmented Jacobian [df/dx  df/dlambda].
%     2. Predict: step a distance ds along that tangent.
%     3. Correct: Newton iterate back onto f = 0, with the extra equation
%        that the correction stays perpendicular to the tangent.
%     4. Record eigenvalues; detect a fold when the tangent's parameter
%        component changes sign, and a Hopf when a complex pair crosses the
%        imaginary axis.
%
%   KNOWN LIMITATION, and it belongs in the report.  The aerodynamic model
%   is a lookup table with linear interpolation, so the vector field is
%   continuous but its Jacobian jumps at every table node.  Newton needs a
%   continuous Jacobian.  In practice the branch tracks well over a range of
%   about 15 deg of alpha and then stalls on a node, despite the adaptive
%   step below.  Kolb does not hit this because the F-18 HARV model he uses
%   is analytic.  This is one concrete reason the literature fits tabulated
%   aerodynamic data to smooth global polynomial models before applying
%   continuation methods, which connects directly to the global nonlinear
%   modelling discussed in Section 3 of the report.

clear; close all; clc;

h    = 0;
xcg  = 0.35;
V0   = 152.4;                       % [m/s] anchor condition
SCALE = [100 1 1 1 1 1 1 1 10]';    % arclength scaling: V in 100 m/s, dh in 10 deg
DS0  = 0.01;                        % nominal arclength step
NMAX = 4000;

% Anchor the branch on a converged trim point.
[thtl, dh0, alpha0, res] = f16_trim(V0, h, 0, xcg);
fprintf('Anchor trim: V = %.1f m/s, alpha = %.2f deg, dh = %+.2f deg, res = %.1e\n', ...
    V0, alpha0, dh0, res);

ar = deg2rad(alpha0);
x0 = [V0; ar; 0; 0; 0; 0; 0; ar];
u0 = [thtl; dh0; 0; 0];

fwd = continue_branch(x0, u0, 2, +DS0, NMAX, h, xcg, SCALE);
bwd = continue_branch(x0, u0, 2, -DS0, NMAX, h, xcg, SCALE);

lam = [fliplr(bwd.lam) fwd.lam];
alp = [fliplr(bwd.alpha) fwd.alpha];
mxr = [fliplr(bwd.maxre) fwd.maxre];
tng = [fliplr(bwd.tlam) fwd.tlam];

fprintf('Branch tracked: %d points, delta_h from %+.2f to %+.2f deg, alpha %.1f to %.1f deg\n', ...
    numel(lam), min(lam), max(lam), min(alp), max(alp));

%% ---- Plot -------------------------------------------------------------
figure('Name','Continuation','Position',[80 80 880 580]); hold on;
st = mxr < 0;
plot(lam(st),  alp(st),  '.', 'MarkerSize', 12, 'DisplayName','stable');
plot(lam(~st), alp(~st), 'x', 'MarkerSize', 6, 'LineWidth', 0.8, ...
    'DisplayName','unstable');
xlabel('Stabilator \delta_h (deg)'); ylabel('Equilibrium \alpha (deg)');
title('Continuation of full eight-state equilibria, throttle held fixed');
legend('Location','best'); grid on; box on;

%% ---- Bifurcation detection --------------------------------------------
nf = 0;
for k = 2:numel(lam)-1
    if tng(k-1)*tng(k) < 0 && abs(lam(k)-lam(k-1)) < 1
        nf = nf + 1;
        fprintf('FOLD near delta_h = %+.2f deg, alpha = %.2f deg\n', lam(k), alp(k));
    end
end
if nf == 0
    fprintf(['No fold detected inside the tracked range. The branch was cut ' ...
             'short by\nthe interpolation limitation described in the header, ' ...
             'not necessarily because\nno fold exists further along.\n']);
end

% ------------------------------------------------------------------------
function out = continue_branch(x0, u0, ip, ds0, nmax, h, xcg, SC)
%CONTINUE_BRANCH  Pseudo-arclength continuation in control component ip.

F  = @(w) resid(w, u0, ip, h, xcg, SC);
w  = [x0; u0(ip)]./SC;

Jw = fdjac(F, w);
[~,~,Vt] = svd(Jw);  t = Vt(:,end);
if sign(t(end)) ~= sign(ds0), t = -t; end

ds = abs(ds0);  fails = 0;
out.lam = [];  out.alpha = [];  out.maxre = [];  out.tlam = [];

for step = 1:nmax
    wp = w + ds*t;  wc = wp;  ok = false;
    for it = 1:30
        Jw = fdjac(F, wc);
        A  = [Jw; t'];
        b  = [-F(wc); -(wc - wp)'*t];
        if rcond(A) < 1e-14, break; end
        dw = A\b;
        wc = wc + dw;
        if norm(dw) < 1e-11, ok = true; break; end
    end

    if ~ok || max(abs(F(wc))) > 1e-6
        ds = ds/2;  fails = fails + 1;
        if ds < 1e-5 || fails > 400, return; end
        continue;
    end
    fails = 0;  ds = min(ds*1.3, abs(ds0));

    Jw = fdjac(F, wc);
    [~,~,Vt] = svd(Jw);  tn = Vt(:,end);
    if tn'*t < 0, tn = -tn; end
    w = wc;  t = tn;

    z = w.*SC;  u = u0;  u(ip) = z(end);
    ev = eig(fdjac(@(xx) f16_dynamics(xx, u, h, xcg), z(1:end-1)));

    out.lam(end+1)   = z(end);
    out.alpha(end+1) = rad2deg(z(2));
    out.maxre(end+1) = max(real(ev));
    out.tlam(end+1)  = t(end);

    if z(end) < -25 || z(end) > 25, return; end
end
end

function r = resid(w, u0, ip, h, xcg, SC)
z = w.*SC;  u = u0;  u(ip) = z(end);
r = f16_dynamics(z(1:end-1), u, h, xcg);
end

function J = fdjac(fun, w)
f0 = fun(w);
n  = numel(w);
J  = zeros(numel(f0), n);
for i = 1:n
    d  = 1e-7*max(1, abs(w(i)));
    wp = w; wp(i) = wp(i) + d;
    J(:,i) = (fun(wp) - f0)/d;
end
end
