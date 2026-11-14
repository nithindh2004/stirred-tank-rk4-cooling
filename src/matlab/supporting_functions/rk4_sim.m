function [t, Y] = rk4_sim(rpm, mdot, p, dt)
%RK4_SIM  Classical fourth-order Runge-Kutta trajectory of (T_L, T_J).
n = round(p.tf/dt); h = p.tf/n; t = (0:n)'*h;
UA = U_coef(rpm, mdot, p, p.Uref)*p.A;
Y = zeros(n+1, 2); Y(1,:) = [p.T0, p.Tin];
f = @(y) [-UA*(y(1)-y(2))/(p.ML*p.CpL); (UA*(y(1)-y(2)) + mdot*p.CpC*(p.Tin-y(2)))/(p.MJ*p.CpC)];
for k = 1:n
    y = Y(k,:)';
    k1 = f(y); k2 = f(y + 0.5*h*k1); k3 = f(y + 0.5*h*k2); k4 = f(y + h*k3);
    Y(k+1,:) = (y + (h/6)*(k1 + 2*k2 + 2*k3 + k4))';
end
end
