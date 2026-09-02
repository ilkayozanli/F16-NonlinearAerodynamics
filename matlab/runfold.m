%RUN_FOLD_VERIFICATION  Verify that the c.g. transition is a saddle-node.
%
%   run_cg_bifurcation locates the centre-of-gravity position at which two
%   high-incidence equilibria appear. That establishes WHERE something
%   happens, but not WHAT. Appearance of a pair of equilibria is consistent
%   with a saddle-node bifurcation but does not by itself demonstrate one:
%   a corner in a piecewise-linear interpolant can produce the same
%   appearance without any smooth fold being present.
%
%   This script applies the two tests that distinguish the two cases.
%
%   TEST 1, TANGENCY.  At a saddle-node the two equilibria merge tangentially,
%   so both C_m = 0 and dC_m/dalpha = 0 at the critical point. The slopes at
%   the two roots must therefore approach zero as the parameter approaches
%   its critical value. At an interpolation corner they remain finite.
%
%   TEST 2, SCALING.  The normal form of a saddle-node gives a root
%   separation proportional to the square root of the distance from the
%   critical parameter value. A factor of ten in that distance should give a
%   factor of about 3.16 in separation. No other bifurcation type, and no
%   interpolation artefact, produces this exponent.
%
%   Both tests are run with the smoothed model (smooth = true). This is
%   essential: with linear interpolation on the raw 5-degree grid, the
%   equilibria merge at alpha = 70 deg, which is itself a table node, and the
%   kink there masks the tangency. The critical c.g. is the same either way,
%   which is itself evidence that the transition is a property of the data
%   rather than of the interpolation.
%
%   Expected output: critical x_cg near 0.3071, slopes falling towards zero,
%   and a scaling exponent near 0.5.

clear; close all; clc;

dh   = 25;          % full nose-down stabilator
V    = 150;
dlef = 0;

%% ---- Locate the critical c.g. by bisection -----------------------------
lo = 0.28;  hi = 0.42;
for k = 1:45
    mid = 0.5*(lo+hi);
    if ~isempty(find_roots(mid, dh, dlef, V)), hi = mid; else, lo = mid; end
end
xcrit = hi;
fprintf('Critical centre of gravity: x_cg = %.5f cbar\n\n', xcrit);

%% ---- Test 1 and 2 ------------------------------------------------------
deltas = [3e-2 1e-2 3e-3 1e-3 3e-4 1e-4];
sep    = nan(size(deltas));
s1     = nan(size(deltas));
s2     = nan(size(deltas));

fprintf('  x_cg - x_crit   separation    |dCm/dalpha| at the two roots\n');
fprintf('  ---------------------------------------------------------\n');
for i = 1:numel(deltas)
    x = xcrit + deltas(i);
    r = find_roots(x, dh, dlef, V);
    if numel(r) < 2, continue; end
    r = [min(r) max(r)];
    sep(i) = r(2) - r(1);
    s1(i)  = abs(slope(r(1), x, dh, dlef, V));
    s2(i)  = abs(slope(r(2), x, dh, dlef, V));
    fprintf('    %8.1e    %7.3f deg    %.2e   %.2e\n', deltas(i), sep(i), s1(i), s2(i));
end

ok = ~isnan(sep);
p  = polyfit(log10(deltas(ok)), log10(sep(ok)), 1);
fprintf('\nScaling exponent of separation with distance from critical: %.3f\n', p(1));
fprintf('A saddle-node bifurcation gives 0.5.\n');

if abs(p(1) - 0.5) < 0.12 && s1(find(ok,1,'last')) < 0.3*s1(find(ok,1))
    fprintf(['\nBoth tests are satisfied: the slopes fall towards zero as the\n' ...
             'critical point is approached, and the separation follows a square\n' ...
             'root law. The transition is a saddle-node bifurcation.\n']);
else
    fprintf(['\nAt least one test is not satisfied. Report the transition as a\n' ...
             'fold-type transition consistent with a saddle-node rather than as\n' ...
             'a verified saddle-node bifurcation.\n']);
end

%% ---- Figure ------------------------------------------------------------
figure('Name','Fold verification','Position',[100 100 900 400]);

subplot(1,2,1);
loglog(deltas(ok), sep(ok), 'o-', 'LineWidth', 1.5); hold on;
ref = sqrt(deltas(ok)/deltas(find(ok,1)))*sep(find(ok,1));
loglog(deltas(ok), ref, 'k--', 'LineWidth', 1.2);
grid on; xlabel('x_{cg} - x_{crit}'); ylabel('root separation (deg)');
title('Separation follows a square root law');
legend({'model', 'square root reference'}, 'Location','southeast');

subplot(1,2,2);
loglog(deltas(ok), s1(ok), 'o-', 'LineWidth', 1.5); hold on;
loglog(deltas(ok), s2(ok), 's-', 'LineWidth', 1.5);
grid on; xlabel('x_{cg} - x_{crit}'); ylabel('|dC_m/d\alpha| at the root');
title('Slopes vanish at the merge');
legend({'lower root','upper root'}, 'Location','southeast');

% ------------------------------------------------------------------------
function r = find_roots(xcg, dh, dlef, V)
g = linspace(55, 90, 1401);
Cm = zeros(size(g));
for i = 1:numel(g), Cm(i) = cm_of(g(i), dh, dlef, V, xcg); end
idx = find(Cm(1:end-1).*Cm(2:end) < 0);
r = zeros(1,numel(idx));
for k = 1:numel(idx)
    r(k) = fzero(@(a) cm_of(a, dh, dlef, V, xcg), [g(idx(k)) g(idx(k)+1)]);
end
end

function s = slope(a, xcg, dh, dlef, V)
h = 0.01;
s = (cm_of(a+h,dh,dlef,V,xcg) - cm_of(a-h,dh,dlef,V,xcg))/(2*h);
end

function Cm = cm_of(alpha, dh, dlef, V, xcg)
% smooth = true is required; see the header.
[~,~,~,~,Cm,~] = f16_aero(alpha, 0, dh, 0, 0, dlef, 0, 0, 0, V, xcg, true);
end