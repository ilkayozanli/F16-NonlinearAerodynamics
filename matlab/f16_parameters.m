clear all;
close all;
%%%--------F-16 High Angle of Attack Analysis---
%% F-16 Specifications
L= 15.027 ; %[m]
b = 9.449; % [m]
h = 5.09; %[m]
empty_weight = 9.207 ; %[kg]
Max_Togw = 21.772 ; %[kg]
n = 9  ; %[g] design load factor
Service_ceiling= 15000; %[m]

run('F16aerodata.m')

% Define input ranges
alpha_range = f16data.alpha1;   % smooth curve
beta = 0;
de = 0;

% Preallocate
CX_vals = zeros(size(alpha_range));


% Loop through alpha values
for i = 1:length(alpha_range)
    alpha = alpha_range(i);

    CX_vals(i) = interpn( ...
        f16data.alpha1, ...
        f16data.beta, ...
        f16data.de1, ...
        f16data.CX, ...
        alpha, beta, de, ...
        'spline');
end

% Plot
subplot(3,1,1);
plot(alpha_range, CX_vals, 'LineWidth', 2);
xlabel('Angle of Attack \alpha (deg) ');
ylabel('C_X');
title('F-16 C_X vs Angle of Attack for \beta=0 deflection=0');
grid on;

%% Deflection effects
de_values = f16data.de1;   % use full range
beta = 0;

subplot(3,1,2); hold on;

for d = 1:length(de_values)
    de = de_values(d);
    CX_vals = zeros(size(alpha_range));

    for i = 1:length(alpha_range)
        alpha = alpha_range(i);

        CX_vals(i) = interpn( ...
            f16data.alpha1, ...
            f16data.beta, ...
            f16data.de1, ...
            f16data.CX, ...
            alpha, beta, de, ...
            'spline');
    end

    plot(alpha_range, CX_vals,'--' ,'LineWidth', 1);
end

xlabel('\alpha (deg)');
ylabel('C_X');
title('Effect of Elevator on C_X (Full Range)');
legend('\delta_e = -25°','-10°','0°','10°','25°');
grid on;
%% Account for beta
beta_values = [-30 -15 0 15 30];
de = 0;

subplot(3,1,3); hold on;

for b = 1:length(beta_values)
    beta = beta_values(b);
    CX_vals = zeros(size(alpha_range));

    for i = 1:length(alpha_range)
        alpha = alpha_range(i);

        CX_vals(i) = interpn( ...
            f16data.alpha1, ...
            f16data.beta, ...
            f16data.de1, ...
            f16data.CX, ...
            alpha, beta, de, ...
            'spline');
    end

    plot(alpha_range, CX_vals, '-.','LineWidth', 1);
end

xlabel('\alpha (deg)');
ylabel('C_X');
title('Effect of Sideslip on C_X');
legend('\beta = -30°','\beta = -15°','\beta = 0°','\beta = 15°','\beta = 30°');
grid on;


