%PLOT_COEFFICIENTS  Nonlinear F-16 aerodynamic coefficients vs angle of attack.
%
%   Replacement for the original f16_parameters.m.  Everything appears in a
%   single window with tabs along the top rather than six separate figure
%   windows.  Tabs, in order:
%
%       C_X, C_Z, C_m       longitudinal coefficients, three panels each
%       Lift and drag       C_L, C_D and L/D in wind axes
%       Lateral-directional CY, Cl, Cn against sideslip
%       Departure           Cn_beta_dynamic against alpha
%       Interpolation       linear against spline, methodology figure
%
%   What changed from the original script, beyond the layout:
%
%   1. Coefficients come from the full buildup in f16_aero, so what is
%      plotted is what the aircraft actually experiences, not a raw table
%      entry.  The c.g. transfer terms, the trim increment, the deep-stall
%      increment and the LEF blending are all included.
%
%   2. Interpolation is linear because that is what the source model
%      specifies.  The last tab quantifies what that choice costs: the
%      difference from a cubic spline peaks near alpha = 85 deg, reaching
%      about 0.16 on CZ at full nose-down stabilator but only about 0.009 on
%      Cm at zero deflection.  Report it that way rather than claiming the
%      spline is badly wrong everywhere, because it is not.
%
%   3. Lateral-directional coefficients are plotted.  Departure and spin are
%      lateral-directional phenomena, so CY, Cl and Cn are not optional here.
%
%   4. Cn_beta_dynamic is computed.  It is the classical departure criterion
%      and the most direct single answer to the stability half of the
%      research question: where it goes negative is the predicted departure
%      boundary.
%
%   Coefficient vector convention: C = [CX CY CZ Cl Cm Cn].
%
%   If you would rather have separate windows again, set SINGLE_WINDOW to
%   false below.

clear; close all; clc;

SINGLE_WINDOW = true;

c = f16_constants();

alpha       = linspace(-20, 90, 441);
dh_values   = [-25 -10 0 10 25];
beta_values = [-30 -15 0 15 30];
V   = 150;
xcg = 0.35;

if SINGLE_WINDOW
    fig = figure('Name','F-16 aerodynamic coefficients', ...
                 'Position',[80 60 900 760],'NumberTitle','off');
    tg  = uitabgroup(fig);
else
    tg = [];
end

%% ---- Longitudinal coefficients -----------------------------------------
long_idx   = [1 3 5];
long_names = {'C_X','C_Z','C_m'};
long_tabs  = {'C_X','C_Z','C_m'};

for ci = 1:3
    idx = long_idx(ci);
    parent = new_page(tg, long_tabs{ci}, long_names{ci});

    % panel 1: baseline
    ax = stack_axes(parent, 3, 1);
    y  = zeros(size(alpha));
    for i = 1:numel(alpha)
        C = coefs(alpha(i),0,0,V,xcg);  y(i) = C(idx);
    end
    plot(ax, alpha, y, 'LineWidth', 1.8); grid(ax,'on');
    xlabel(ax,'\alpha (deg)'); ylabel(ax,long_names{ci});
    title(ax, sprintf('%s vs \\alpha  (\\beta = 0, \\delta_h = 0, x_{cg} = %.2f)', ...
        long_names{ci}, xcg));

    % panel 2: stabilator effect
    ax = stack_axes(parent, 3, 2); hold(ax,'on');
    for d = dh_values
        for i = 1:numel(alpha)
            C = coefs(alpha(i),0,d,V,xcg);  y(i) = C(idx);
        end
        plot(ax, alpha, y, 'LineWidth', 1.2);
    end
    grid(ax,'on'); xlabel(ax,'\alpha (deg)'); ylabel(ax,long_names{ci});
    title(ax, ['Effect of stabilator on ' long_names{ci}]);
    legend(ax, arrayfun(@(d) sprintf('\\delta_h = %+d\\circ', d), dh_values, ...
        'UniformOutput', false), 'Location','best');

    % panel 3: sideslip effect
    ax = stack_axes(parent, 3, 3); hold(ax,'on');
    for b = beta_values
        for i = 1:numel(alpha)
            C = coefs(alpha(i),b,0,V,xcg);  y(i) = C(idx);
        end
        plot(ax, alpha, y, 'LineWidth', 1.2);
    end
    grid(ax,'on'); xlabel(ax,'\alpha (deg)'); ylabel(ax,long_names{ci});
    title(ax, ['Effect of sideslip on ' long_names{ci}]);
    legend(ax, arrayfun(@(b) sprintf('\\beta = %+d\\circ', b), beta_values, ...
        'UniformOutput', false), 'Location','best');
end

%% ---- Lateral-directional coefficients vs sideslip -----------------------
parent = new_page(tg, 'Lateral-directional', 'Lateral-directional');
alpha_cases = [0 15 30 45 60];
beta_fine   = linspace(-30, 30, 121);
lat_names   = {'C_Y','C_l','C_n'};
lat_idx     = [2 4 6];
yb = zeros(size(beta_fine));

for li = 1:3
    ax = stack_axes(parent, 3, li); hold(ax,'on');
    for a = alpha_cases
        for i = 1:numel(beta_fine)
            C = coefs(a, beta_fine(i), 0, V, xcg);  yb(i) = C(lat_idx(li));
        end
        plot(ax, beta_fine, yb, 'LineWidth', 1.2);
    end
    grid(ax,'on'); xlabel(ax,'\beta (deg)'); ylabel(ax,lat_names{li});
    title(ax, [lat_names{li} ' vs sideslip at several angles of attack']);
    legend(ax, arrayfun(@(a) sprintf('\\alpha = %d\\circ', a), alpha_cases, ...
        'UniformOutput', false), 'Location','best');
end

%% ---- Lift and drag: the wind-axis view ---------------------------------
% The dataset contains body-axis forces only.  Lift and drag are wind-axis
% quantities and must be obtained by rotating through the angle of attack:
%
%     C_N = -C_Z          normal force        C_A = -C_X   axial force
%     C_L =  C_N cos(alpha) - C_A sin(alpha)
%     C_D =  C_N sin(alpha) + C_A cos(alpha)
%
% This matters for the argument of the report.  A linear model assumes C_L
% proportional to alpha and C_D quadratic in C_L.  Neither holds here: the
% lift curve bends over and reverses, and the drag rise is not quadratic.
CL = zeros(size(alpha));  CD = zeros(size(alpha));
for i = 1:numel(alpha)
    C  = coefs(alpha(i), 0, 0, V, xcg);
    CN = -C(3);  CA = -C(1);
    ar = deg2rad(alpha(i));
    CL(i) = CN*cos(ar) - CA*sin(ar);
    CD(i) = CN*sin(ar) + CA*cos(ar);
end

parent = new_page(tg, 'Lift and drag', 'Lift and drag');

ax = stack_axes(parent, 3, 1);
plot(ax, alpha, CL, 'LineWidth', 1.8); grid(ax,'on');
hold(ax,'on');
lin = alpha >= 0 & alpha <= 10;
pc  = polyfit(alpha(lin), CL(lin), 1);
plot(ax, alpha, polyval(pc, alpha), '--', 'LineWidth', 1.2);
ylim(ax, [min(CL)-0.2, max(CL)+0.4]);
xlabel(ax,'\alpha (deg)'); ylabel(ax,'C_L');
legend(ax, {'C_L from the model', ...
    sprintf('linear extrapolation of the 0-10\\circ slope (%.3f /deg)', pc(1))}, ...
    'Location','southeast');
[CLmax, im] = max(CL);
title(ax, sprintf(['Lift curve: departs from linear beyond about 10\\circ, ' ...
                   'C_{L,max} = %.2f at \\alpha = %.0f\\circ'], CLmax, alpha(im)));

ax = stack_axes(parent, 3, 2);
plot(ax, alpha, CD, 'LineWidth', 1.8); grid(ax,'on');
xlabel(ax,'\alpha (deg)'); ylabel(ax,'C_D');
title(ax,'Drag coefficient: monotonic rise, not a quadratic in C_L');

ax = stack_axes(parent, 3, 3);
plot(ax, alpha, CL./max(CD,1e-6), 'LineWidth', 1.8); grid(ax,'on');
[LDmax, il] = max(CL./max(CD,1e-6));
xlabel(ax,'\alpha (deg)'); ylabel(ax,'L/D');
title(ax, sprintf('Lift-to-drag ratio, maximum %.1f at \\alpha = %.0f\\circ', ...
    LDmax, alpha(il)));

fprintf('C_L,max = %.3f at alpha = %.1f deg;  (L/D)_max = %.1f at alpha = %.1f deg\n', ...
    CLmax, alpha(im), LDmax, alpha(il));

%% ---- Departure criterion: Cn_beta_dynamic ------------------------------
% Cn_beta_dyn = Cn_beta*cos(alpha) - (Izz/Ixx)*Cl_beta*sin(alpha)
% Negative means directional divergence, i.e. nose-slice departure.
db  = 1;                                  % [deg] sideslip perturbation
Cnb = zeros(size(alpha));
Clb = zeros(size(alpha));
for i = 1:numel(alpha)
    Cp  = coefs(alpha(i), +db, 0, V, xcg);
    Cn_ = coefs(alpha(i), -db, 0, V, xcg);
    Cnb(i) = (Cp(6) - Cn_(6))/(2*db);     % [1/deg]
    Clb(i) = (Cp(4) - Cn_(4))/(2*db);
end
ar      = deg2rad(alpha);
Cnb_dyn = Cnb.*cos(ar) - (c.Izz/c.Ixx)*Clb.*sin(ar);

parent = new_page(tg, 'Departure', 'Departure criterion');
ax = stack_axes(parent, 1, 1);
plot(ax, alpha, Cnb_dyn, 'LineWidth', 1.8); hold(ax,'on');
plot(ax, [alpha(1) alpha(end)], [0 0], 'k--', 'LineWidth', 1);
grid(ax,'on'); xlabel(ax,'\alpha (deg)'); ylabel(ax,'C_{n\beta,dyn}  (1/deg)');
title(ax, 'Dynamic directional stability: negative means departure-prone');

k = find(Cnb_dyn(1:end-1) > 0 & Cnb_dyn(2:end) < 0, 1);
if ~isempty(k)
    fprintf('C_n_beta_dyn first goes negative near alpha = %.1f deg.\n', alpha(k));
end

%% ---- Interpolation method comparison -----------------------------------
% The primary reason for using linear interpolation is that it is what the
% source model specifies, so matching it keeps the results comparable with
% the published F-16 model.  The plots below are the secondary, empirical
% argument, and they are shown honestly: the difference between linear and
% spline is small for most slices and only becomes substantial for a few.
%
% The largest deviations in the whole table set occur near alpha = 85 deg,
% in the sparse 70-80-90 region.  On CZ at full nose-down stabilator the
% difference reaches about 0.16, which is a large fraction of the
% coefficient itself.  On the Cm slice at beta = 0, delta_h = 0 it is only
% about 0.009, which is why that slice looks identical either way.

run('F16aerodata.m');
a_tab  = f16data.alpha1;
a_fine = linspace(-20, 90, 1101);
j0 = find(f16data.beta == 0);

CZ_raw = squeeze(f16data.CZ(:, j0, f16data.de1 == 25));   % worst case
Cm_raw = squeeze(f16data.Cm(:, j0, f16data.de1 ==  0));   % typical case

parent = new_page(tg, 'Interpolation', 'Interpolation');

% panel 1: the case where the choice visibly matters
ax = stack_axes(parent, 2, 1); hold(ax,'on');
plot(ax, a_tab, CZ_raw, 'ko', 'MarkerSize', 6, 'DisplayName','table points');
plot(ax, a_fine, interp1(a_tab, CZ_raw, a_fine, 'linear'), ...
    'LineWidth', 1.6, 'DisplayName','linear (used here)');
plot(ax, a_fine, interp1(a_tab, CZ_raw, a_fine, 'spline'), '--', ...
    'LineWidth', 1.6, 'DisplayName','spline');
grid(ax,'on'); xlabel(ax,'\alpha (deg)'); ylabel(ax,'C_Z table value');
title(ax, 'C_Z at \beta = 0, \delta_h = +25\circ: the largest disagreement in the table set');
legend(ax, 'Location','best');

% panel 2: how big the difference is, and where
d_CZ = interp1(a_tab, CZ_raw, a_fine, 'spline') - interp1(a_tab, CZ_raw, a_fine, 'linear');
d_Cm = interp1(a_tab, Cm_raw, a_fine, 'spline') - interp1(a_tab, Cm_raw, a_fine, 'linear');

ax = stack_axes(parent, 2, 2); hold(ax,'on');
plot(ax, a_fine, d_CZ, 'LineWidth', 1.6, ...
    'DisplayName','C_Z, \delta_h = +25\circ');
plot(ax, a_fine, d_Cm, 'LineWidth', 1.6, ...
    'DisplayName','C_m, \delta_h = 0\circ');
plot(ax, [a_fine(1) a_fine(end)], [0 0], 'k--', 'LineWidth', 0.8, ...
    'HandleVisibility','off');
grid(ax,'on'); xlabel(ax,'\alpha (deg)'); ylabel(ax,'spline minus linear');
title(ax, 'Difference concentrates in the sparse 70-80-90 deg region');
legend(ax, 'Location','best');

fprintf(['Interpolation check: max |spline - linear| is %.4f for C_Z at ' ...
         'delta_h = +25 deg\n   and %.4f for C_m at delta_h = 0 deg, ' ...
         'both peaking near alpha = %.0f deg.\n'], ...
         max(abs(d_CZ)), max(abs(d_Cm)), a_fine(find(abs(d_CZ)==max(abs(d_CZ)),1)));

% ========================================================================
function parent = new_page(tg, tab_title, fig_name)
%NEW_PAGE  A tab in the shared window, or a fresh figure if tabs are off.
if isempty(tg)
    parent = figure('Name', fig_name, 'Position', [80 80 760 720]);
else
    parent = uitab(tg, 'Title', tab_title);
end
end

function ax = stack_axes(parent, n, i)
%STACK_AXES  Axes i of n stacked vertically, i = 1 at the top.
%   Positions are set explicitly rather than via subplot or tiledlayout so
%   that this works inside a uitab on any MATLAB version.
left   = 0.10;
width  = 0.84;
bottom = 0.07;
top    = 0.05;
slot   = (1 - bottom - top)/n;
if n == 1
    height = slot*0.86;      % a lone panel should fill the tab
else
    height = slot*0.66;      % stacked panels need room for titles and labels
end
pos_y  = bottom + (n - i)*slot + (slot - height)*0.55;
ax = axes('Parent', parent, 'Position', [left pos_y width height]);
end

function C = coefs(alpha, beta, dh, V, xcg)
%COEFS  Row vector [CX CY CZ Cl Cm Cn] at zero rates, clean configuration.
[CX,CY,CZ,Cl,Cm,Cn] = f16_aero(alpha, beta, dh, 0, 0, 0, 0, 0, 0, V, xcg);
C = [CX CY CZ Cl Cm Cn];
end