
clear all;
close all;

%%%--------F-16 High Angle of Attack Analysis---
%% F-16 Specifications
L= 15.027 ; %[m]
b= 9.449;% [m]
h = 5.09; %[m]
empty_weight = 9.207 ; %[kg]
Max_Togw = 21.772 ; %[kg]
n = 9  ; %[g] design load factor
Service_ceiling= 15000; %[m]


% Define input ranges
alpha_range = -20:1:90;   % smooth curve
beta = 0;
de = 0;

% Preallocate
CX_vals = zeros(size(alpha_range));
run('F16aerodata.m')

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
figure;
plot(alpha_range, CX_vals, 'LineWidth', 2);
xlabel('Angle of Attack \alpha (deg)');
ylabel('C_X');
title('F-16 C_X vs Angle of Attack');
grid on;



