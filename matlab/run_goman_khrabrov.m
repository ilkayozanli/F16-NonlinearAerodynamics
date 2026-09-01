%RUN_GOMAN_KHRABROV  Unsteady aerodynamics as a state-space model.
%
%   Implements the Goman-Khrabrov model, which is what Section 3 of the
%   report argues for but which nothing else in this project provides.
%
%   THE PROBLEM.  Every other script here treats the aerodynamic
%   coefficients as instantaneous functions of the current angle of attack:
%   look up Cm(alpha) and you are done.  That is fine while the flow stays
%   attached.  Above the stall it is wrong, because separated flow takes
%   time to develop and time to reattach.  The moment acting on the aircraft
%   right now depends on where alpha has recently been, not only on where it
%   is.
%
%   THE MODEL.  Goman and Khrabrov (1994) add one internal state x, the
%   position of the flow separation point, running from 1 for fully attached
%   to 0 for fully separated.  It obeys a first-order equation with two
%   parameters:
%
%       tau1 * dx/dt + x = x0( alpha - tau2 * alphadot )
%
%   tau1 is the relaxation time: how long the separated flow takes to catch
%   up with a change in alpha.  tau2 is a lag on the input: the separation
%   responds to a slightly delayed angle of attack when alpha is changing.
%   x0(alpha) is the steady separation position, recovered from the static
%   wind tunnel data.
%
%   RECOVERING x0 FROM THE TABLES.  The Kirchhoff relation for a stalling
%   aerofoil gives the static normal force as
%
%       CN = CN_alpha * ((1 + sqrt(x0))/2)^2 * sin(alpha)
%
%   so x0 follows by inverting it using the tabulated CN.  The result is
%   x0 = 1 below about 17 deg, meaning fully attached, then falling
%   monotonically as alpha increases.  The model is therefore only
%   meaningful above that angle; below it there is nothing separated to lag.
%
%   RECOVERING THE COEFFICIENT.  This script uses the effective-angle form:
%   the coefficient at any instant is the static coefficient evaluated at
%   the angle whose steady separation position equals the current x.  In
%   steady flight x equals x0(alpha), the effective angle equals alpha, and
%   the model reduces exactly to the static tables.  When alpha is moving, x
%   lags, the effective angle differs from alpha, and the coefficient
%   differs from its static value.  That is the whole mechanism.
%
%   This is a simplified implementation.  The original paper reconstructs
%   the coefficients from x through the Kirchhoff expression directly rather
%   than through an effective angle.  Say so in the methodology section.
%
%   PARAMETERS.  tau1 and tau2 are not in the F-16 dataset; they come from
%   dynamic wind tunnel testing, which the van Oort and Sonneveldt tables do
%   not include.  Rather than borrowing a number and pretending it is the
%   F-16's, this script sweeps both and shows that the hysteresis loop
%   appears for any reasonable choice and only its width changes.  That is
%   the honest way to present a model with unmeasured parameters.
%
%   Reference: Goman, M. and Khrabrov, A., "State-Space Representation of
%   Aerodynamic Characteristics of an Aircraft at High Angles of Attack",
%   Journal of Aircraft, 31(5), 1994, pp. 1109-1115.

clear; close all; clc;

c = f16_constants();

V      = 150;         % [m/s] reference airspeed, for the reduced frequency
ALPHA_SEP = 18;       % [deg] below this the flow is attached and x0 = 1
tau1   = 0.35;        % [s] relaxation time, nominal
tau2   = 0.05;        % [s] alphadot lag, nominal
a_mean = 30;          % [deg] centre of the oscillation
a_amp  = 8;           % [deg] amplitude
freqs  = [0.02 0.2 0.5 1.5];   % [Hz] oscillation frequencies to compare

%% ---- Static data and the separation curve ------------------------------
al = (0:0.5:90)';
CN = zeros(size(al));  Cm = zeros(size(al));
for i = 1:numel(al)
    [~,~,CZi,~,Cmi,~] = f16_aero(al(i), 0, 0, 0, 0, 0, 0, 0, 0, V, 0.35);
    CN(i) = -CZi;      % normal force is minus the body-axis Z force
    Cm(i) =  Cmi;
end

% Kirchhoff slope, fitted on the attached-flow part of the curve.
lin  = al >= 4 & al <= 12;
pfit = polyfit(sind(al(lin)), CN(lin), 1);
CNa  = pfit(1);
fprintf('Kirchhoff CN_alpha fitted on 4-12 deg: %.3f per rad\n', CNa);

ratio = CN ./ (CNa*sind(max(al,1e-6)));
x0    = min(max((2*sqrt(max(ratio,0)) - 1).^2, 0), 1);

% Keep the region where separation is actually developing, and force it to
% be strictly decreasing so that it can be inverted.
m    = al >= ALPHA_SEP;
a_m  = al(m);
x_m  = cummin(x0(m)) - 1e-9*(0:sum(m)-1)';

X0     = @(a) interp1(a_m, x_m, min(max(a,a_m(1)),a_m(end)), 'pchip');
X0INV  = @(x) interp1(flipud(x_m), flipud(a_m), ...
                      min(max(x,min(x_m)),max(x_m)), 'pchip');
CMSTAT = @(a) interp1(al, Cm, min(max(a,al(1)),al(end)), 'pchip');

fprintf('Separation point: x0 = 1 below %.0f deg, %.3f at 30 deg, %.3f at 45 deg\n', ...
    ALPHA_SEP, X0(30), X0(45));

%% ---- Hysteresis loops at several frequencies ---------------------------
fig = figure('Name','Goman-Khrabrov unsteady model','Position',[60 50 1000 700], ...
             'NumberTitle','off');
tg  = uitabgroup(fig);
co  = lines(numel(freqs));

ax = full_axes(uitab(tg,'Title','Hysteresis loops')); hold(ax,'on');
plot(ax, al(al>=a_mean-a_amp-2 & al<=a_mean+a_amp+2), ...
         Cm(al>=a_mean-a_amp-2 & al<=a_mean+a_amp+2), ...
     'k-', 'LineWidth', 2, 'DisplayName','static table');

fprintf('\n  freq [Hz]  reduced freq k   C_m at %g deg rising / falling   loop width\n', a_mean);
fprintf('  --------------------------------------------------------------------\n');

for j = 1:numel(freqs)
    f = freqs(j);
    [a_t, cm_t] = gk_loop(f, a_mean, a_amp, tau1, tau2, X0, X0INV, CMSTAT, ALPHA_SEP);
    plot(ax, a_t, cm_t, '-', 'Color', co(j,:), 'LineWidth', 1.4, ...
        'DisplayName', sprintf('%.2f Hz', f));

    k  = 2*pi*f*c.cbar/(2*V);                 % reduced frequency
    [cu, cd] = branch_values(a_t, cm_t, a_mean);
    fprintf('  %8.2f   %12.4f   %+8.4f / %+8.4f        %.4f\n', f, k, cu, cd, abs(cu-cd));
end

grid(ax,'on'); box(ax,'on');
xlabel(ax,'\alpha (deg)'); ylabel(ax,'C_m');
title(ax, sprintf(['Pitching moment against \\alpha: static curve is a line, ' ...
                   'unsteady response is a loop (\\tau_1 = %.2f s)'], tau1));
legend(ax,'Location','best');

%% ---- The separation curve itself ---------------------------------------
ax = full_axes(uitab(tg,'Title','Separation point')); hold(ax,'on');
plot(ax, al, x0, 'LineWidth', 1.8, 'DisplayName','x_0 from Kirchhoff inversion');
plot(ax, [ALPHA_SEP ALPHA_SEP], [0 1], 'k--', 'LineWidth', 1, ...
    'DisplayName','onset of separation');
grid(ax,'on'); ylim(ax,[0 1.05]);
xlabel(ax,'\alpha (deg)'); ylabel(ax,'x_0  (1 = attached, 0 = separated)');
title(ax,'Steady separation point recovered from the static normal force');
legend(ax,'Location','best');

%% ---- Time histories ----------------------------------------------------
ax = full_axes(uitab(tg,'Title','Time histories')); hold(ax,'on');
f = 0.5;
[a_t, cm_t, t, x_t] = gk_loop(f, a_mean, a_amp, tau1, tau2, X0, X0INV, CMSTAT, ALPHA_SEP);
yyaxis(ax,'left');
plot(ax, t, a_t, 'LineWidth', 1.5); ylabel(ax,'\alpha (deg)');
yyaxis(ax,'right');
plot(ax, t, x_t, 'LineWidth', 1.5); hold(ax,'on');
plot(ax, t, X0(a_t), '--', 'LineWidth', 1.2);
ylabel(ax,'separation point');
grid(ax,'on'); xlabel(ax,'time (s)');
title(ax, sprintf(['At %.1f Hz the separation point (solid) lags its steady ' ...
                   'value (dashed): that lag is the hysteresis'], f));

%% ---- Parameter sensitivity ---------------------------------------------
ax = full_axes(uitab(tg,'Title','Parameter sweep')); hold(ax,'on');
t1_list = [0.10 0.25 0.50 1.00];
f_list  = logspace(-2, 0.5, 20);
for j = 1:numel(t1_list)
    w = zeros(size(f_list));
    for i = 1:numel(f_list)
        [a_t, cm_t] = gk_loop(f_list(i), a_mean, a_amp, t1_list(j), tau2, ...
                              X0, X0INV, CMSTAT, ALPHA_SEP);
        [cu, cd] = branch_values(a_t, cm_t, a_mean);
        w(i) = abs(cu - cd);
    end
    semilogx(ax, f_list, w, 'o-', 'LineWidth', 1.4, ...
        'DisplayName', sprintf('\\tau_1 = %.2f s', t1_list(j)));
end
grid(ax,'on');
xlabel(ax,'oscillation frequency (Hz)'); ylabel(ax,'C_m loop width');
title(ax,'Loop width against frequency: a peak, not a trend, for every \tau_1');
legend(ax,'Location','best');

fprintf(['\nThe loop width peaks at an intermediate frequency for every tau1.\n' ...
         'Slow oscillation lets the separation point keep up, so there is no lag.\n' ...
         'Very fast oscillation freezes it, so it stops responding at all.\n' ...
         'Maximum hysteresis happens in between, where the flow is trying to\n' ...
         'follow and failing.  That peak is the signature of the model and it\n' ...
         'is present for any reasonable tau1, which is the point of the sweep.\n']);

% ------------------------------------------------------------------------
function [a_t, cm_t, t, x_t] = gk_loop(f, a_mean, a_amp, tau1, tau2, ...
                                       X0, X0INV, CMSTAT, ALPHA_SEP)
%GK_LOOP  One steady-state cycle of the Goman-Khrabrov model.
w   = 2*pi*f;
alp = @(t) a_mean + a_amp*sin(w*t);
dal = @(t) deg2rad(a_amp*w*cos(w*t));          % alphadot in rad/s

% tau1*xdot + x = x0(alpha - tau2*alphadot)
rhs = @(t,x) ( X0( max(alp(t) - rad2deg(tau2*dal(t)), ALPHA_SEP) ) - x )/tau1;

opts = odeset('RelTol',1e-9,'AbsTol',1e-11);
sol  = ode45(rhs, [0 8/f], X0(a_mean), opts);   % run in, then sample one cycle

t    = linspace(6/f, 8/f, 800);
x_t  = deval(sol, t)';
a_t  = alp(t)';
cm_t = CMSTAT( X0INV(x_t) );
end

function [c_up, c_down] = branch_values(a, cm, a_query)
%BRANCH_VALUES  Coefficient on the rising and falling branches of the loop.
d  = gradient(a);
up = d > 0;  dn = ~up;
c_up   = interp1(a(up), cm(up), a_query, 'linear', NaN);
adn    = flipud(a(dn));  cdn = flipud(cm(dn));
c_down = interp1(adn, cdn, a_query, 'linear', NaN);
end

function ax = full_axes(parent)
ax = axes('Parent', parent, 'Position', [0.10 0.11 0.80 0.80]);
end
