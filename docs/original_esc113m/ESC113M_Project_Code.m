% Cooling of Hot Liquid in a Stirred Tank: Solved using RK4 & compared to Analytical

%  Parameters 
U = 150;              % Overall heat transfer coefficient (W/m^2·K)
A = 2;                % Heat transfer area (m^2)
rho = 1000;           % Density of liquid (kg/m^3)
V = 0.5;              % Volume of liquid (m^3)
Cp = 4180;            % Specific heat capacity (J/kg·K)
Tc = 300;             % Coolant temperature (K)
T0 = 370;             % Initial liquid temperature (K)

% Derived constant
k = (U*A)/(rho*V*Cp);  % 1/s

% Time settings
t0 = 0; tf = 1800; dt = 10;
t = t0:dt:tf;
N = length(t);
T = zeros(size(t));
T(1) = T0;

% ODE definition
f = @(T) -k*(T - Tc);

%  RK4 Integration 
for i = 1:N-1
    k1 = f(T(i));
    k2 = f(T(i) + 0.5*dt*k1);
    k3 = f(T(i) + 0.5*dt*k2);
    k4 = f(T(i) + dt*k3);
    
    T(i+1) = T(i) + (dt/6)*(k1 + 2*k2 + 2*k3 + k4);
end

%  Analytical Solution 
T_exact = Tc + (T0 - Tc) * exp(-k * t);

%  Print & Compare Selected Values
% Select 10 equally spaced points
n_points = 10;
idx = round(linspace(1, N, n_points));
t_selected = t(idx);
T_rk4_selected = T(idx);
T_exact_selected = T_exact(idx);
abs_error_selected = abs(T_rk4_selected - T_exact_selected);

% Print table
fprintf('\nComparison of RK4 and Analytical Solutions at Selected Times:\n');
fprintf('Time (s)    T (RK4)    T (Analytical)    Absolute Error\n');
fprintf('--------    --------   --------------    --------------\n');
for i = 1:n_points
    fprintf('%8.0f    %9.3f      %12.3f      %10.5f\n', ...
        t_selected(i), T_rk4_selected(i), T_exact_selected(i), abs_error_selected(i));
end

% Max error
max_error = max(abs(T - T_exact));
fprintf('\nMaximum absolute error over all points: %.16f K\n', max_error);

%  Plot: RK4 vs Analytical
figure;
plot(t, T, 'b-', 'LineWidth', 2); hold on;
plot(t, T_exact, 'r--', 'LineWidth', 2);
plot(t_selected, T_rk4_selected, 'ko', 'MarkerFaceColor', 'k');  % markers

% Annotate selected points
for i = 1:n_points
    text(t_selected(i), T_rk4_selected(i)+0.5, ...
        sprintf('%.1f', T_rk4_selected(i)), ...
        'FontSize', 8, 'HorizontalAlignment', 'center');
end

legend('RK4 Numerical', 'Analytical', 'Selected Points', 'Location', 'best');
xlabel('Time (s)');
ylabel('Temperature (K)');
title('Cooling of a Hot Liquid in a Stirred Tank: RK4 vs Analytical');
grid on;
