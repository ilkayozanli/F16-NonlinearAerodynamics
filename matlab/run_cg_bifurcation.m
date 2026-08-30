%RUN_CG_BIFURCATION  Continuation in c.g. position: where the deep stall appears.
%
%   Companion to run_pitch_equilibria.  Same root-finding, different sweep.
%
%   run_pitch_equilibria fixes the c.g. and sweeps the stabilator.  This
%   script does the opposite: it pins the stabilator at full nose-down and
%   sweeps the c.g.  That single change turns a qualitative comparison of
%   three c.g. cases into one number.
%
%   Why full nose-down.  With delta_h held at +25 deg the pilot has already
%   applied all the nose-down authority the aircraft has.  Any statically
%   stable equilibrium that exists at high angle of attack under that
%   condition is therefore an equilibrium the pilot cannot fly out of using
%   the stabilator.  The c.g. at which such an equilibrium first appears is
%   the critical c.g. for the locked-in deep stall.
%
%   What you should see: nothing at high alpha for a forward c.g., then a
%   stable and an unstable equilibrium appearing together at a well-defined
%   c.g. and separating as the c.g. moves aft.  Two equilibria being born
%   together as a parameter crosses a threshold is a saddle-node (fold)
%   bifurcation.  This is the cleanest bifurcation in the whole model and
%   the easiest one to write up.
%
%   The script also bisects to the critical c.g. and prints it.

clear; close all; clc;

dh   = 25;                          % [deg] full nose-down stabilator
V    = 150;                         % [m/s]
dlef = 0;                           % leading edge flap retracted

xcg_sweep  = 0.28:0.002:0.40;       % c.g. as a fraction of cbar
alpha_grid = linspace(-20, 90, 441);
ALPHA_HIGH = 40;                    % [deg] threshold for "high alpha"

branch_x = []; branch_a = []; branch_s = [];

fprintf('C.g. continuation at delta_h = %+d deg (full nose-down), V = %.0f m/s\n', dh, V);
fprintf('--------------------------------------------------------------------\n');

for xcg = xcg_sweep
    Cm = zeros(size(alpha_grid));
    for i = 1:numel(alpha_grid)
        Cm(i) = cm_of(alpha_grid(i), dh, dlef, V, xcg);
    end

    idx = find(Cm(1:end-1).*Cm(2:end) < 0);

    for i = idx
        a0 = fzero(@(a) cm_of(a, dh, dlef, V, xcg), ...
                   [alpha_grid(i) alpha_grid(i+1)]);
        hs    = 0.05;
        slope = (cm_of(a0+hs,dh,dlef,V,xcg) - cm_of(a0-hs,dh,dlef,V,xcg))/(2*hs);

        branch_x(end+1) = xcg;    %#ok<SAGROW>
        branch_a(end+1) = a0;     %#ok<SAGROW>
        branch_s(end+1) = slope;  %#ok<SAGROW>
    end
end

%% ---- Bisect to the critical c.g. ---------------------------------------
% Bracket: the smallest swept c.g. with a high-alpha equilibrium, and the
% largest one without.
has = arrayfun(@(x) any(branch_x==x & branch_a>ALPHA_HIGH), xcg_sweep);
k   = find(has, 1);

if isempty(k) || k == 1
    xcg_crit = NaN;
    fprintf('No clean threshold inside the swept range.\n');
else
    lo = xcg_sweep(k-1);   hi = xcg_sweep(k);
    for n = 1:40
        mid = 0.5*(lo + hi);
        if deep_exists(mid, dh, dlef, V, ALPHA_HIGH)
            hi = mid;
        else
            lo = mid;
        end
    end
    xcg_crit = hi;
    fprintf(['Critical c.g. for the deep-stall equilibrium pair: x_cg = %.4f c-bar\n' ...
             'Forward of this the aircraft always has nose-down authority.\n' ...
             'Aft of it a statically stable high-alpha trim exists at full\n' ...
             'nose-down stabilator, which the pilot cannot fly out of.\n'], xcg_crit);
end

%% ---- Plot --------------------------------------------------------------
figure('Name','C.g. bifurcation','Position',[100 100 880 580]); hold on;

st = branch_s < 0;
plot(branch_x(st),  branch_a(st),  '.', 'MarkerSize', 12, ...
    'DisplayName','stable equilibrium');
plot(branch_x(~st), branch_a(~st), 'x', 'MarkerSize', 6, 'LineWidth', 0.8, ...
    'DisplayName','unstable equilibrium');

if ~isnan(xcg_crit)
    yl = ylim;
    plot([xcg_crit xcg_crit], [-20 90], 'k--', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('saddle-node at x_{cg} = %.3f', xcg_crit));
    ylim(yl);
end

xlabel('Centre of gravity  x_{cg}  (fraction of mean aerodynamic chord)');
ylabel('Equilibrium angle of attack \alpha (deg)');
title(sprintf(['Deep-stall onset: equilibria at full nose-down stabilator ' ...
               '(\\delta_h = %+d\\circ)'], dh));
legend('Location','northwest'); grid on; box on;
ylim([-20 90]);

%% ---- Caveat worth repeating in the report ------------------------------
fprintf(['\nCaveat for the write-up: the stable branch sits above 70 deg, where\n' ...
         'the alpha table has nodes only at 70, 80 and 90 deg. The existence of\n' ...
         'the branch and the critical c.g. are robust; the exact angle is not.\n' ...
         'The branch also leaves the swept range at the aft end because it runs\n' ...
         'past the 90 deg table limit, which is a data boundary, not physics.\n']);

% ------------------------------------------------------------------------
function tf = deep_exists(xcg, dh, dlef, V, a_hi)
%DEEP_EXISTS  True if any equilibrium sits above a_hi at this c.g.
g = linspace(a_hi, 90, 401);
Cm = zeros(size(g));
for i = 1:numel(g)
    Cm(i) = cm_of(g(i), dh, dlef, V, xcg);
end
tf = any(Cm(1:end-1).*Cm(2:end) < 0);
end

function Cm = cm_of(alpha, dh, dlef, V, xcg)
[~,~,~,~,Cm,~] = f16_aero(alpha, 0, dh, 0, 0, dlef, 0, 0, 0, V, xcg);
end
