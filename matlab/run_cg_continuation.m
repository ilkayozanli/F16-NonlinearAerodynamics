%RUN_CG_CONTINUATION  Pseudo-arclength continuation of the deep-stall branch.
%
%   Replaces the parameter sweep of run_cg_bifurcation with a genuine
%   continuation method, applied to the pitch equilibrium condition
%
%       C_m(alpha, x_cg) = 0
%
%   with the stabilator held at its full nose-down limit.
%
%   WHY CONTINUATION RATHER THAN A SWEEP.  A parameter sweep prescribes
%   x_cg and solves for alpha at each step. It cannot follow a branch that
%   turns back on itself, because at the turning point there is no solution
%   for the next value of the parameter. The fold therefore has to be
%   inferred from the appearance and disappearance of solutions, and located
%   by bisection on a yes/no test.
%
%   Pseudo-arclength continuation treats x_cg as a further unknown and steps
%   along the solution curve by arclength rather than by parameter. The fold
%   is then simply the point at which the tangent's parameter component
%   changes sign, and the branch continues through it in a single connected
%   trace. This is the method used by Kolb (2017) and implemented in
%   continuation packages such as MATCONT.
%
%   Each step:
%     1. Tangent to the solution curve from the null space of the Jacobian
%        [dCm/dalpha, dCm/dx_cg], oriented to continue the previous
%        direction.
%     2. Predictor: step a distance ds along the tangent.
%     3. Corrector: Newton iterate back onto C_m = 0, with the extra
%        condition that the correction is orthogonal to the tangent.
%     4. Adapt ds: halve on failure, grow slowly on success. This is
%        essential near the fold, where the curvature is high and a fixed
%        step overshoots.
%
%   SMOOTHING IS REQUIRED. The corrector is a Newton method and needs a
%   continuous Jacobian. With linear interpolation on the raw 5-degree grid
%   the slope jumps at every table node and the corrector stalls. f16_aero is
%   therefore called with smooth = true throughout.
%
%   Expected result: a single connected branch from x_cg = 0.42 down around a
%   fold near x_cg = 0.3071, alpha = 69.8 deg, and back out along the
%   unstable branch. The fold location agrees with the bisection result of
%   run_cg_bifurcation.

clear; close all; clc;

dh    = 25;                 % full nose-down stabilator [deg]
V     = 150;                % [m/s]
dlef  = 0;
SC    = [50; 0.02];         % arclength scaling: alpha in 50 deg, x_cg in 0.02
DS0   = 0.02;               % nominal arclength step
XLIM  = [0.30 0.42];        % stop outside this range of x_cg
ALIM  = [20 92];            % stop outside this range of alpha

%% ---- Starting point on the stable high-incidence branch ----------------
x_start = 0.35;
a_start = fzero(@(a) cm_of(a, dh, dlef, V, x_start), [70 85]);
fprintf('Starting from x_cg = %.3f, alpha = %.2f deg\n', x_start, a_start);

w = [a_start; x_start]./SC;
Jc = jac(w, SC, dh, dlef, V);
t  = [Jc(2); -Jc(1)];  t = t/norm(t);
if t(2) > 0, t = -t; end        % march towards decreasing x_cg first

%% ---- Continuation ------------------------------------------------------
n    = 0;
ds   = DS0;
pts  = zeros(3, 20000);         % [x_cg; alpha; tangent parameter component]

for k = 1:20000
    wp = w + ds*t;
    wc = wp;
    ok = false;
    for it = 1:40
        Jc = jac(wc, SC, dh, dlef, V);
        A  = [Jc(:)'; t(:)'];
        b  = [-resid(wc, SC, dh, dlef, V); -(wc - wp)'*t];
        if rcond(A) < 1e-15, break; end
        dw = A\b;
        wc = wc + dw;
        if norm(dw) < 1e-13, ok = true; break; end
    end

    if ~ok || abs(resid(wc, SC, dh, dlef, V)) > 1e-9
        ds = ds/2;
        if ds < 1e-7, break; end
        continue;
    end
    ds = min(ds*1.25, DS0);

    Jc = jac(wc, SC, dh, dlef, V);
    tn = [Jc(2); -Jc(1)];  tn = tn/norm(tn);
    if tn'*t < 0, tn = -tn; end
    w = wc;  t = tn;

    z = w.*SC;
    n = n + 1;
    pts(:,n) = [z(2); z(1); t(2)];

    if z(2) < XLIM(1) || z(2) > XLIM(2) || z(1) < ALIM(1) || z(1) > ALIM(2)
        break;
    end
end
pts = pts(:,1:n);

fprintf('Traced %d points: x_cg from %.4f to %.4f, alpha from %.1f to %.1f deg\n', ...
    n, min(pts(1,:)), max(pts(1,:)), min(pts(2,:)), max(pts(2,:)));

%% ---- Fold detection ----------------------------------------------------
% A fold is where the tangent's parameter component changes sign: the branch
% stops advancing in x_cg and turns back.
turn = find(pts(3,1:end-1).*pts(3,2:end) < 0) + 1;
fprintf('\nTurning points (folds):\n');
for i = turn
    fprintf('   x_cg = %.5f, alpha = %.2f deg\n', pts(1,i), pts(2,i));
end

%% ---- Stability of each point ------------------------------------------
stab = false(1,n);
for i = 1:n
    h = 0.02;
    s = (cm_of(pts(2,i)+h, dh, dlef, V, pts(1,i)) - ...
         cm_of(pts(2,i)-h, dh, dlef, V, pts(1,i)))/(2*h);
    stab(i) = s < 0;
end

%% ---- Plot --------------------------------------------------------------
figure('Name','C.g. continuation','Position',[100 100 900 560],'Color','w');
hold on;
plot(pts(1,:), pts(2,:), '-', 'Color', [.7 .7 .7], 'LineWidth', 1, ...
    'HandleVisibility','off');
plot(pts(1,stab),  pts(2,stab),  '.', 'MarkerSize', 11, ...
    'DisplayName','stable equilibrium');
plot(pts(1,~stab), pts(2,~stab), '.', 'MarkerSize', 11, ...
    'DisplayName','unstable equilibrium');
for i = turn
    plot(pts(1,i), pts(2,i), 'ko', 'MarkerSize', 9, 'LineWidth', 1.5, ...
        'DisplayName', sprintf('fold at x_{cg} = %.4f', pts(1,i)));
end
grid on; box on;
xlabel('Centre of gravity  x_{cg}  (fraction of mean aerodynamic chord)');
ylabel('Equilibrium angle of attack \alpha (deg)');
title(sprintf(['Pseudo-arclength continuation of the equilibrium branch ' ...
               'at \\delta_h = %+d\\circ'], dh));
legend('Location','best');

% ------------------------------------------------------------------------
function r = resid(w, SC, dh, dlef, V)
z = w.*SC;
r = cm_of(z(1), dh, dlef, V, z(2));
end

function J = jac(w, SC, dh, dlef, V)
e = 1e-7;
J = zeros(1,2);
for i = 1:2
    wp = w; wp(i) = wp(i) + e;
    wm = w; wm(i) = wm(i) - e;
    J(i) = (resid(wp,SC,dh,dlef,V) - resid(wm,SC,dh,dlef,V))/(2*e);
end
end

function Cm = cm_of(alpha, dh, dlef, V, xcg)
% smooth = true is required: the Newton corrector needs a continuous Jacobian.
[~,~,~,~,Cm,~] = f16_aero(alpha, 0, dh, 0, 0, dlef, 0, 0, 0, V, xcg, true);
end
