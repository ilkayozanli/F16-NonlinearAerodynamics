clear all;
close all;
%%%--------F-16 High Angle of Attack Analysis---
%% F-16 Specifications
L= 15.027 ; %[m]
b = 9.449; % [m]
h = 5.09; %[m]
empty_weight = 9207 ; %[kg]
Max_Togw = 21772 ; %[kg]
n = 9  ; %[g] design load factor
Service_ceiling= 15000; %[m]

run('F16aerodata.m')

alpha_range = f16data.alpha1;   % smooth curve
de_values   = f16data.de1;      % full elevator/stabilator range
beta_values = [-30 -15 0 15 30];

% ------------------------------------------------------------------
% Coefficients to process: name, data field, y-axis label, plot title tag
% ------------------------------------------------------------------
coeffs = {
    'CX', 'C_X', 'C_X'
    'CZ', 'C_Z', 'C_Z'
    'Cm', 'C_m', 'C_m'
};

for c = 1:size(coeffs,1)
    field_name = coeffs{c,1};
    y_label    = coeffs{c,2};
    title_tag  = coeffs{c,3};

    data = f16data.(field_name);   % CX, CZ, or Cm array

    figure('Name', field_name);

    % --- Subplot 1: baseline curve, beta = 0, de = 0 ---
    beta = 0; de = 0;
    vals = zeros(size(alpha_range));
    for i = 1:length(alpha_range)
        vals(i) = interpn(f16data.alpha1, f16data.beta, f16data.de1, ...
            data, alpha_range(i), beta, de, 'spline');
    end
    subplot(3,1,1);
    plot(alpha_range, vals, 'LineWidth', 2);
    xlabel('Angle of Attack \alpha (deg)');
    ylabel(y_label);
    title(['F-16 ' title_tag ' vs Angle of Attack for \beta=0 deflection=0']);
    grid on;

    % --- Subplot 2: effect of elevator/stabilator deflection ---
    beta = 0;
    subplot(3,1,2); hold on;
    for d = 1:length(de_values)
        de = de_values(d);
        vals = zeros(size(alpha_range));
        for i = 1:length(alpha_range)
            vals(i) = interpn(f16data.alpha1, f16data.beta, f16data.de1, ...
                data, alpha_range(i), beta, de, 'spline');
        end
        plot(alpha_range, vals, '--', 'LineWidth', 1);
    end
    xlabel('\alpha (deg)');
    ylabel(y_label);
    title(['Effect of Elevator on ' title_tag ' (Full Range)']);
    legend('\delta_e = -25°','-10°','0°','10°','25°');
    grid on;

    % --- Subplot 3: effect of sideslip ---
    de = 0;
    subplot(3,1,3); hold on;
    for b = 1:length(beta_values)
        beta = beta_values(b);
        vals = zeros(size(alpha_range));
        for i = 1:length(alpha_range)
            vals(i) = interpn(f16data.alpha1, f16data.beta, f16data.de1, ...
                data, alpha_range(i), beta, de, 'spline');
        end
        plot(alpha_range, vals, '-.', 'LineWidth', 1);
    end
    xlabel('\alpha (deg)');
    ylabel(y_label);
    title(['Effect of Sideslip on ' title_tag]);
    legend('\beta = -30°','\beta = -15°','\beta = 0°','\beta = 15°','\beta = 30°');
    grid on;
end