%RUN_MODAL_ANALYSIS  Dynamic stability of the F-16 across the alpha range.
%
%   Stability from the eigenvalues of the linearised eight-state dynamics,
%   with modes identified by EIGENVECTOR PARTICIPATION rather than by
%   sorting eigenvalues by size.
%
%   Why that matters.  Sorting by magnitude relabels the modes whenever two
%   of them swap frequency, so a single physical motion gets drawn as one
%   curve at low alpha and a different curve at high alpha, and every swap
%   looks like a jump.  It also miscounts bifurcations, because a relabel
%   looks like a damping sign change.  Kolb (2017) instead reports which
%   state components each eigenvector actually contains and names the mode
%   from that: his Hopf involves Mach, alpha, q and theta so he calls it the
%   phugoid, and his branch point involves beta, p, r and phi so he calls it
%   the spiral mode.  This script does the same.
%
%   Classification, from the normalised eigenvector (airspeed weighted by
%   1/100 so its large numerical size does not dominate):
%     longitudinal share = V, alpha, q, theta      lateral share = beta, p, r, phi
%     longitudinal + oscillatory + V heavy         -> phugoid
%     longitudinal + q and alpha heavy             -> short period
%     lateral + oscillatory                        -> dutch roll
%     lateral + real + roll rate heavy             -> roll subsidence
%     lateral + real + bank heavy                  -> spiral
%
%   EQUILIBRIUM FAMILY.  Level-flight trim cannot exceed alpha = 7 deg,
%   because holding altitude at higher alpha means flying slower than the
%   aircraft can sustain.  The throttle is therefore held fixed and the
%   flight path angle allowed to float, so the aircraft descends and the
%   family reaches alpha = 40 deg.  Same choice Kolb makes.
%
%   SMOOTHING.  f16_dynamics is called with smooth = true.  Eigenvalues
%   depend on the SLOPE of the tables, and plain linear interpolation makes
%   that slope jump at every 5 degree node.  With smoothing off the growth
%   rate jumps from 0.33 to 1.12 across alpha = 10 deg for no physical
%   reason.  See the header of f16_aero.
%
%   Expected results at x_cg = 0.35, idle throttle, sea level:
%     alpha  0 to 15 deg   unstable, real longitudinal root, by design
%     alpha 15 to 32 deg   stable without augmentation
%     alpha ~22.5 deg      roll subsidence and spiral merge into a slow,
%                          heavily damped oscillation (roll-spiral coupling)
%     alpha  ~32 deg       dutch roll crosses to negative damping: HOPF,
%                          period about 5.6 s, which is wing rock
%     alpha  ~34 deg       real root in beta and r goes positive:
%                          directional divergence, nose slice departure

clear; close all; clc;

h    = 0;
xcg  = 0.35;
thtl = 0.0;
V_list = [200:-2.5:100, 99:-1:70, 69.8:-0.2:50];

% Modes 1-3 oscillate, so damping ratio and natural frequency mean something
% for them.  Modes 4-7 are real roots: they have no frequency and no damping
% ratio, so they are reported by their real part instead.  Plotting zeta for
% a real root gives exactly +-1 and is meaningless.
MODES = {'short period','phugoid','dutch roll','roll-spiral oscillation', ...
         'roll subsidence','spiral','pitch divergence','pitch convergence'};
OSC   = [true true true true false false false false];
NM    = numel(MODES);

n = numel(V_list);
alpha_tr = nan(1,n); dh_tr = nan(1,n); gam_tr = nan(1,n);
maxre    = nan(1,n); worst_lat = false(1,n); worst_osc = false(1,n);
WN = nan(NM,n); ZT = nan(NM,n); RE = nan(NM,n);

guess = [0 1 0];

fprintf('Throttle %.2f, sea level, x_cg %.2f\n', thtl, xcg);
fprintf('  V [m/s]  alpha    gamma   delta_h   max Re    worst mode\n');
fprintf('  ---------------------------------------------------------------\n');

for k = 1:n
    V = V_list(k);
    [dh, alpha, gam, res] = f16_trim_glide(V, thtl, h, xcg, guess);
    if res > 1e-7, continue; end

    % Reject branch jumps: the family is smooth, so a large step in the
    % controls means the solver has hopped onto a different branch.
    if k > 1 && ~isnan(dh_tr(find(~isnan(dh_tr),1,'last')))
        last = find(~isnan(dh_tr),1,'last');
        if abs(dh - dh_tr(last)) > 3 || abs(alpha - alpha_tr(last)) > 5
            fprintf('  %6.1f   branch jump detected, sweep stopped\n', V);
            break;
        end
    end
    guess = [dh alpha gam];

    ar = deg2rad(alpha);
    x  = [V; ar; 0; 0; 0; 0; 0; ar + deg2rad(gam)];
    u  = [thtl; dh; 0; 0];

    [Vec, Lam] = eig(numerical_jacobian(@(xx) f16_dynamics(xx, u, h, xcg, true), x));
    ev = diag(Lam);

    alpha_tr(k) = alpha;  dh_tr(k) = dh;  gam_tr(k) = gam;
    [maxre(k), iw] = max(real(ev));

    [cls, islat] = classify_modes(ev, Vec);
    worst_lat(k) = islat(iw);
    worst_osc(k) = abs(imag(ev(iw))) > 1e-6;

    for m = 1:NM
        j = find(cls == m & imag(ev) >= 0, 1);   % one of each conjugate pair
        if isempty(j), continue; end
        RE(m,k) = real(ev(j));
        if OSC(m)
            WN(m,k) = abs(ev(j));
            ZT(m,k) = -real(ev(j))/abs(ev(j));
        end
    end

    fprintf('  %6.1f  %6.2f  %7.2f  %+7.2f  %+8.4f   %s\n', ...
        V, alpha, gam, dh, maxre(k), MODES{cls(iw)});
end

ok = ~isnan(alpha_tr);
a  = alpha_tr(ok);

%% ---- One window, five tabs ---------------------------------------------
fig = figure('Name','F-16 modal stability','Position',[60 50 1000 700], ...
             'NumberTitle','off');
tg = uitabgroup(fig);
co = lines(NM);

% --- tab 1: growth rate, coloured by which mode is worst ---
ax = full_axes(uitab(tg,'Title','Growth rate')); hold(ax,'on');
plot(ax, a, maxre(ok), '-', 'Color', [.6 .6 .6], 'LineWidth', 1);
lonI = ok & ~worst_lat;  latI = ok & worst_lat;
plot(ax, alpha_tr(lonI), maxre(lonI), 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', co(1,:), 'MarkerEdgeColor','none', ...
    'DisplayName','longitudinal mode is worst');
plot(ax, alpha_tr(latI), maxre(latI), 'o', 'MarkerSize', 4, ...
    'MarkerFaceColor', co(3,:), 'MarkerEdgeColor','none', ...
    'DisplayName','lateral mode is worst');
xl = xlim(ax); plot(ax, xl, [0 0], 'k--','LineWidth',1,'HandleVisibility','off');
xlim(ax,xl); grid(ax,'on');
xlabel(ax,'\alpha (deg)'); ylabel(ax,'max Re(\lambda)  (1/s)');
title(ax,'Growth rate, coloured by which family the unstable mode belongs to');
legend(ax,'Location','best');

% --- tab 2: damping, oscillatory modes only ---
ax = full_axes(uitab(tg,'Title','Damping')); hold(ax,'on');
for m = find(OSC)
    if all(isnan(ZT(m,ok))), continue; end
    plot(ax, a, ZT(m,ok), '.', 'MarkerSize', 10, 'Color', co(m,:), ...
        'DisplayName', MODES{m});
end
xl = xlim(ax); plot(ax, xl, [0 0], 'k--','LineWidth',1,'HandleVisibility','off');
xlim(ax,xl); grid(ax,'on');
xlabel(ax,'\alpha (deg)'); ylabel(ax,'\zeta');
title(ax,'Damping ratio by identified mode: crossing below zero is a Hopf');
legend(ax,'Location','best');

% --- tab 3: frequency, oscillatory modes only ---
ax = full_axes(uitab(tg,'Title','Frequencies')); hold(ax,'on');
for m = find(OSC)
    if all(isnan(WN(m,ok))), continue; end
    plot(ax, a, WN(m,ok), '.', 'MarkerSize', 10, 'Color', co(m,:), ...
        'DisplayName', MODES{m});
end
grid(ax,'on'); xlabel(ax,'\alpha (deg)'); ylabel(ax,'\omega_n (rad/s)');
title(ax,'Natural frequency, oscillatory modes only');
legend(ax,'Location','best');

% --- tab 4: real roots, which have no frequency or damping ratio ---
ax = full_axes(uitab(tg,'Title','Real roots')); hold(ax,'on');
for m = find(~OSC)
    if all(isnan(RE(m,ok))), continue; end
    plot(ax, a, RE(m,ok), '.', 'MarkerSize', 10, 'Color', co(m,:), ...
        'DisplayName', MODES{m});
end
xl = xlim(ax); plot(ax, xl, [0 0], 'k--','LineWidth',1,'HandleVisibility','off');
xlim(ax,xl); grid(ax,'on');
xlabel(ax,'\alpha (deg)'); ylabel(ax,'Re(\lambda)  (1/s)');
title(ax,'Non-oscillatory modes: a crossing above zero is a divergence');
legend(ax,'Location','best');

% --- tab 5: root locus ---
ax = full_axes(uitab(tg,'Title','Root locus')); hold(ax,'on');
for m = 1:NM
    re = RE(m,ok);
    if all(isnan(re)), continue; end
    if OSC(m)
        im = WN(m,ok).*sqrt(max(1 - ZT(m,ok).^2, 0));
    else
        im = zeros(size(re));
    end
    plot(ax, re, im, '.', 'MarkerSize', 9, 'Color', co(m,:), 'DisplayName', MODES{m});
    if OSC(m)
        plot(ax, re,-im, '.', 'MarkerSize', 9, 'Color', co(m,:), 'HandleVisibility','off');
    end
end
yl = ylim(ax); plot(ax,[0 0],yl,'k--','LineWidth',1,'HandleVisibility','off');
ylim(ax,yl); grid(ax,'on'); box(ax,'on');
xlabel(ax,'Real part (1/s)'); ylabel(ax,'Imaginary part (rad/s)');
title(ax,'Root locus with modes identified by eigenvector participation');
legend(ax,'Location','best');

% --- tab 6: trim conditions ---
ax = full_axes(uitab(tg,'Title','Trim conditions'));
yyaxis(ax,'left');  plot(ax, a, dh_tr(ok), '.-','LineWidth',1);
ylabel(ax,'\delta_h (deg)');
yyaxis(ax,'right'); plot(ax, a, gam_tr(ok), '.-','LineWidth',1);
ylabel(ax,'\gamma (deg)');
grid(ax,'on'); xlabel(ax,'\alpha (deg)');
title(ax,'Stabilator and flight path angle along the equilibrium family');

%% ---- Report -----------------------------------------------------------
fprintf('\n=== Stability transitions ===\n');
g = maxre(ok);
for i = find(g(1:end-1).*g(2:end) < 0)
    fprintf('Overall stability changes sign between alpha = %.1f and %.1f deg.\n', ...
        a(i), a(i+1));
end

fprintf('\n=== Bifurcations, by identified mode ===\n');
for m = find(OSC)
    z = ZT(m,ok);
    for i = find(~isnan(z(1:end-1)) & ~isnan(z(2:end)) & z(1:end-1).*z(2:end) < 0)
        w = WN(m,i);
        fprintf(['HOPF: %s damping crosses zero between alpha = %.1f and ' ...
                 '%.1f deg,\n      frequency %.2f rad/s, limit cycle period ' ...
                 'about %.1f s.\n'], MODES{m}, a(i), a(i+1), w, 2*pi/w);
    end
end
for m = find(~OSC)
    r = RE(m,ok);
    for i = find(~isnan(r(1:end-1)) & ~isnan(r(2:end)) & r(1:end-1).*r(2:end) < 0)
        fprintf(['DIVERGENCE: %s real root crosses zero between alpha = %.1f ' ...
                 'and %.1f deg.\n'], MODES{m}, a(i), a(i+1));
    end
end

% ------------------------------------------------------------------------
function [cls, islat] = classify_modes(ev, Vec)
%CLASSIFY_MODES  Name each eigenvalue from its eigenvector participation.
%   1 short period              2 phugoid            3 dutch roll
%   4 roll-spiral oscillation    5 roll subsidence    6 spiral
%   7 pitch divergence           8 pitch convergence
%
%   Above about alpha = 22.5 deg the roll subsidence and spiral roots
%   coalesce into a fourth oscillation, slow and heavily damped.  It needs
%   its own category: without one it is silently discarded, because only the
%   first eigenvalue matching each category is recorded, and the frequency
%   and real-root plots go blank over that whole band.
%
%   Modes within a family are separated by SPEED rather than by which state
%   component happens to be largest.  Roll subsidence is the fast real
%   lateral root and spiral the slow one; dutch roll is the fast lateral
%   oscillation and roll-spiral the slow one; short period is the fast
%   longitudinal oscillation and phugoid the slow one.  Comparing individual
%   eigenvector components instead is fragile: bank angle can outweigh roll
%   rate in the roll subsidence eigenvector at moderate alpha, which put a
%   spurious fast branch on the spiral curve.
%
%   Real lateral roots are separated by SPEED, not by which of p or phi is
%   larger.  The bank angle can outweigh the roll rate in the roll subsidence
%   eigenvector at moderate alpha, which sent it into the spiral bin and put
%   a spurious fast branch on the spiral curve.  The physical distinction is
%   that roll subsidence is the fast root and spiral the slow one.
w   = [1/100 1 1 1 1 1 1 1]';      % weight V down; it is numerically large
n   = numel(ev);
cls = zeros(n,1);
islat = false(n,1);

for i = 1:n
    v   = abs(w.*Vec(:,i));  v = v/sum(v);
    lon = v(1) + v(2) + v(5) + v(8);      % V, alpha, q, theta
    lat = v(3) + v(4) + v(6) + v(7);      % beta, p, r, phi
    osc = abs(imag(ev(i))) > 1e-6;
    islat(i) = lat > lon;

    if ~islat(i)
        if osc, cls(i) = 1; else, cls(i) = 7; end   % split by speed below
    else
        if osc, cls(i) = 3; else, cls(i) = 5; end
    end
end

% Longitudinal oscillations: fastest is short period, slowest is phugoid.
j = find(cls == 1 & imag(ev) > 1e-6);
if numel(j) > 1
    [~,o] = sort(abs(ev(j)), 'descend');
    cls(j(o(1)))     = 1;
    cls(j(o(2:end))) = 2;
    for kk = 1:numel(j)                   % keep conjugates consistent
        cls(conj_index(ev, j(kk))) = cls(j(kk));
    end
elseif numel(j) == 1 && abs(ev(j)) < 0.6
    cls(j) = 2;  cls(conj_index(ev,j)) = 2;
end

% Lateral oscillations: fastest is dutch roll, slowest is roll-spiral.
j = find(cls == 3 & imag(ev) > 1e-6);
if numel(j) > 1
    [~,o] = sort(abs(ev(j)), 'descend');
    cls(j(o(1)))     = 3;
    cls(j(o(2:end))) = 4;
    for kk = 1:numel(j)
        cls(conj_index(ev, j(kk))) = cls(j(kk));
    end
end

% Real lateral roots: fastest is roll subsidence, slowest is spiral.
j = find(cls == 5);
if numel(j) > 1
    [~,o] = sort(abs(real(ev(j))), 'descend');
    cls(j(o(1)))     = 5;
    cls(j(o(2:end))) = 6;
elseif numel(j) == 1 && abs(real(ev(j))) < 0.2
    cls(j) = 6;
end

% Real longitudinal roots: the more positive one is the divergence.
j = find(cls == 7);
if numel(j) > 1
    [~,o] = sort(real(ev(j)), 'descend');
    cls(j(o(1)))     = 7;
    cls(j(o(2:end))) = 8;
end
end

function k = conj_index(ev, i)
%CONJ_INDEX  Index of the complex conjugate partner of eigenvalue i.
[~,k] = min(abs(ev - conj(ev(i))));
if k == i
    d = abs(ev - conj(ev(i)));  d(i) = Inf;  [~,k] = min(d);
end
end

function ax = full_axes(parent)
ax = axes('Parent', parent, 'Position', [0.10 0.11 0.80 0.80]);
end

function J = numerical_jacobian(fun, x)
n = numel(x);  J = zeros(n);
for i = 1:n
    d  = 1e-6*max(1, abs(x(i)));
    xp = x; xp(i) = xp(i) + d;
    xm = x; xm(i) = xm(i) - d;
    J(:,i) = (fun(xp) - fun(xm))/(2*d);
end
end