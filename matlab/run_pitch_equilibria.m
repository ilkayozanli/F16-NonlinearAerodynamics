%RUN_PITCH_EQUILIBRIA  Pitch-axis equilibrium map of the nonlinear F-16.

%   This is pitch-moment equilibrium at a specified flight condition, not
%   full flight trim.  Use f16_trim for that.  The two answer different
%   questions: this one maps where the aircraft can sit in pitch, which is
%   what the deep-stall argument needs; f16_trim gives the controls required
%   to actually fly at a given speed.
%
%   RUNTIME.  Set QUALITY below.  'draft' takes well under a minute and is
%   what you want while iterating.  'final' is about four times slower and
%   gives smooth curves for the report.  Both find the same equilibria; only
%   the resolution of the plotted branches differs.

clear; close all; clc;

QUALITY = 'final';         % 'draft' or 'final'

switch QUALITY
    case 'draft'
        dh_sweep   = -25:1:25;
        alpha_grid = linspace(-20, 90, 441);
    case 'final'
        dh_sweep   = -25:0.25:25;
        alpha_grid = linspace(-20, 90, 1101);
    otherwise
        error('QUALITY must be ''draft'' or ''final''.');
end

V    = 150;                       % [m/s] flight condition for the rate terms
dlef = 0;                         % leading edge flap retracted
xcg_cases = [0.30 0.35 0.38];     % c.g. positions [fraction of cbar]

colors = lines(numel(xcg_cases));
figure('Name','Pitch equilibrium map','Position',[100 100 900 600]); hold on;

fprintf('Pitch equilibria, V = %.0f m/s, beta = 0, LEF = %.0f deg, quality = %s\n', ...
    V, dlef, QUALITY);
fprintf('----------------------------------------------------------------------\n');

results = struct('xcg',{},'dh',{},'alpha',{},'dCm_dalpha',{},'stable',{});
t0 = tic;

for k = 1:numel(xcg_cases)
    xcg = xcg_cases(k);
    branch_a = []; branch_d = []; branch_s = [];

    for dh = dh_sweep
        % Vectorised evaluation of Cm over the alpha grid.
        Cm = zeros(size(alpha_grid));
        for i = 1:numel(alpha_grid)
            Cm(i) = cm_of(alpha_grid(i), dh, dlef, V, xcg);
        end

        idx = find(Cm(1:end-1).*Cm(2:end) < 0);

        for i = idx
            a0 = fzero(@(a) cm_of(a, dh, dlef, V, xcg), ...
                       [alpha_grid(i) alpha_grid(i+1)]);

            hstep = 0.05;
            slope = (cm_of(a0+hstep,dh,dlef,V,xcg) - ...
                     cm_of(a0-hstep,dh,dlef,V,xcg))/(2*hstep);

            branch_d(end+1) = dh;      %#ok<SAGROW>
            branch_a(end+1) = a0;      %#ok<SAGROW>
            branch_s(end+1) = slope;   %#ok<SAGROW>

            results(end+1) = struct('xcg',xcg,'dh',dh,'alpha',a0, ...
                'dCm_dalpha',slope,'stable',slope<0);  %#ok<SAGROW>
        end
    end

    st = branch_s < 0;
    plot(branch_d(st),  branch_a(st),  '.', 'Color', colors(k,:), ...
        'MarkerSize', 12, 'DisplayName', sprintf('x_{cg} = %.2f, stable', xcg));
    plot(branch_d(~st), branch_a(~st), 'x', 'Color', colors(k,:), ...
        'MarkerSize', 5, 'LineWidth', 0.5, ...
        'DisplayName', sprintf('x_{cg} = %.2f, unstable', xcg));

    n_deep = sum(branch_a > 40 & st);
    fprintf('x_cg = %.2f c-bar : %4d equilibria, %3d of them stable above alpha = 40 deg   [%.1f s]\n', ...
        xcg, numel(branch_a), n_deep, toc(t0));
end

xlabel('Stabilator deflection \delta_h (deg)   [positive = trailing edge down]');
ylabel('Equilibrium angle of attack \alpha (deg)');
title('F-16 pitch-axis equilibria: dots are statically stable, crosses unstable');
legend('Location','northwest'); grid on; box on;
ylim([-20 90]);

fprintf('Total sweep time: %.1f s\n', toc(t0));

% ---- Report a representative deep-stall point ---------------------------
deep = results([results.xcg]==0.35 & [results.alpha]>40 & [results.stable]);
if ~isempty(deep)
    [~,j] = max([deep.dh]);
    fprintf(['\nDeep-stall equilibrium at x_cg = 0.35: alpha = %.1f deg with the\n' ...
             'stabilator at %+.1f deg, i.e. full nose-down control applied and the\n' ...
             'aircraft still sits in a statically stable trim. dCm/dalpha = %.4f per deg.\n'], ...
             deep(j).alpha, deep(j).dh, deep(j).dCm_dalpha);
end

% ---- Verification against the flight-trim solver ------------------------
[thtl, dh_t, a_t, res] = f16_trim(152.4, 0, 0, 0.35);
fprintf(['\nCheck case, level flight at 152.4 m/s (500 ft/s), sea level, x_cg = 0.35:\n' ...
         '   alpha = %.2f deg, delta_h = %+.2f deg, throttle = %.3f, residual = %.1e\n'], ...
         a_t, dh_t, thtl, res);

% ------------------------------------------------------------------------
function Cm = cm_of(alpha, dh, dlef, V, xcg)
[~,~,~,~,Cm,~] = f16_aero(alpha, 0, dh, 0, 0, dlef, 0, 0, 0, V, xcg);
end