%% Trim search: Cm(alpha) = 0 for each stabilator deflection

clear all;
close all;

run('F16aerodata.m')

beta = 0;
de_values = f16data.de1;

alpha_data = f16data.alpha1;

% Find beta index corresponding to beta = 0
[~, beta_index] = min(abs(f16data.beta - beta));

% Fine alpha grid
alpha_fine = linspace(min(alpha_data), max(alpha_data), 2000);

fprintf('Trim search results (beta = 0):\n');
fprintf('--------------------------------\n');

figure;
hold on;

trim_results = struct('de', {}, 'alpha_trim', {}, 'slope', {});

for d = 1:length(de_values)

    de = de_values(d);

    % Find stabilator-deflection index
    [~, de_index] = min(abs(f16data.de1 - de));

    % ---------------------------------------------------------
    % Extract the actual 1-D Cm(alpha) data
    % beta = 0, delta_h = current deflection
    % ---------------------------------------------------------
    Cm_data = squeeze(f16data.Cm(:, beta_index, de_index));
     % SHOW ORIGINAL DATA POINTS
    plot(alpha_data, Cm_data, 'ko', ...
        'MarkerSize', 4, ...
        'HandleVisibility', 'off');

    % 1-D PCHIP interpolation
    Cm_fine = interp1(alpha_data, Cm_data, alpha_fine, 'pchip');

    % Plot interpolated curve
    plot(alpha_fine, Cm_fine, ...
        'DisplayName', sprintf('\\delta_h = %d°', de));

    % ---------------------------------------------------------
% Find trim points: Cm = 0
% ---------------------------------------------------------
trims_this_de = [];

for i = 1:length(alpha_fine)-1

    % Case 1: point is already extremely close to zero
    if abs(Cm_fine(i)) < 1e-6

        trims_this_de(end+1) = alpha_fine(i);

    % Case 2: Cm changes sign between two points
    elseif Cm_fine(i) * Cm_fine(i+1) < 0

        % 1-D PCHIP function
        f = @(a) interp1(alpha_data, Cm_data, a, 'pchip');

        % Find exact zero
        a_root = fzero(f, ...
            [alpha_fine(i), alpha_fine(i+1)]);

        trims_this_de(end+1) = a_root;
    end
end
    % ---------------------------------------------------------
    % Calculate dCm/dalpha at every trim point
    % ---------------------------------------------------------
    for k = 1:length(trims_this_de)

        a0 = trims_this_de(k);
        fprintf('FOUND TRIM: delta_h = %+4d deg, alpha = %.3f deg\n', de, a0);

        h = 0.01;   % degrees

        Cm_minus = interp1(alpha_data, Cm_data, ...
            a0-h, 'pchip');

        Cm_plus = interp1(alpha_data, Cm_data, ...
            a0+h, 'pchip');

        slope = (Cm_plus - Cm_minus)/(2*h);

        % Static longitudinal stability
        if slope < 0
            stability_str = 'STABLE';
        else
            stability_str = 'UNSTABLE';
        end

        fprintf(['delta_h = %+4d deg  ->  ' ...
                 'alpha_trim = %6.2f deg   ' ...
                 'dCm/dalpha = %8.4f  [%s]\n'], ...
                 de, a0, slope, stability_str);

        trim_results(end+1) = struct( ...
            'de', de, ...
            'alpha_trim', a0, ...
            'slope', slope);
    end
end

% Zero Cm line
yline(0, 'k-', 'LineWidth', 1);

xlabel('\alpha (deg)');
ylabel('C_m');
title('Pitching Moment vs \alpha, with Trim Points');

legend show;
grid on;

% -------------------------------------------------------------
% Mark trim points on the plot
% -------------------------------------------------------------

for k = 1:length(trim_results)

    if trim_results(k).slope < 0
        marker = 'go';       % stable
    else
        marker = 'rx';       % unstable
    end

    plot(trim_results(k).alpha_trim, 0, marker, ...
        'MarkerSize', 6, ...
        'LineWidth', 2, ...
        'HandleVisibility', 'off');

end