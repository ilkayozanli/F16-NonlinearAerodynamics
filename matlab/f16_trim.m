function [thtl, dh, alpha, resid] = f16_trim(V, h, gamma_deg, xcg, guess)
%F16_TRIM  Steady straight-flight trim of the nonlinear F-16 model.
%
%   [thtl, dh, alpha, resid] = F16_TRIM(V, h, gamma_deg, xcg, guess)
%
%   Solves the three longitudinal equilibrium equations simultaneously:
%
%       Vdot     = 0      (thrust balances drag along the flight path)
%       alphadot = 0      (lift balances weight normal to it)
%       qdot     = 0      (pitching moment about the c.g. is zero)
%
%   for throttle, stabilator and angle of attack, with sideslip and all
%   body rates zero.  This is what "trim" means.  Solving only qdot = 0,
%   which is the same as finding the zeros of Cm(alpha), gives points that
%   satisfy moment balance but at which the airplane is not in equilibrium
%   at all: the reported alpha will generally not hold the aircraft up.
%
%   Inputs:  V         true airspeed                       [m/s]
%            h         altitude                            [m]      (default 0)
%            gamma_deg flight path angle                   [deg]    (default 0)
%            xcg       c.g. as a fraction of cbar          [-]      (default 0.35)
%            guess     [thtl dh alpha] initial guess                (optional)
%
%   Outputs: throttle [0..1], stabilator [deg], alpha [deg], and the largest
%   absolute residual, which should be at machine precision on convergence.
%
%   Check case:  V = 500 ft/s = 152.4 m/s, sea level, level flight,
%   xcg = 0.35 gives approximately alpha = 2.0 deg, dh = -0.2 deg,
%   throttle = 0.13.  If your run does not reproduce this, the aerodynamic
%   buildup or the constants are wrong; fix that before going further.

if nargin < 2 || isempty(h),         h = 0;            end
if nargin < 3 || isempty(gamma_deg), gamma_deg = 0;    end
if nargin < 4 || isempty(xcg),       xcg = 0.35;       end
if nargin < 5 || isempty(guess),     guess = [0.2 -1 3]; end

c = f16_constants();
gam = deg2rad(gamma_deg);

fun = @(u) trim_residual(u, V, h, gam, xcg, c);

if exist('fsolve','file') == 2
    opts = optimoptions('fsolve','Display','off', ...
                        'FunctionTolerance',1e-12,'StepTolerance',1e-12);
    u = fsolve(fun, guess, opts);
else
    % No Optimization Toolbox: minimise the sum of squares instead.
    cost = @(u) sum(fun(u).^2);
    u = fminsearch(cost, guess, optimset('Display','off', ...
                   'TolX',1e-10,'TolFun',1e-14,'MaxFunEvals',5e4,'MaxIter',5e4));
end

thtl  = u(1);  dh = u(2);  alpha = u(3);
resid = max(abs(fun(u)));

end

% ------------------------------------------------------------------------
function res = trim_residual(u, V, h, gam, xcg, c)

thtl  = min(max(u(1),0),1);
dh    = u(2);
alpha = u(3);

[rho, a] = f16_engine('atmos', h);
qbar = 0.5*rho*V^2;
mach = V/a;
T    = f16_engine('thrust', f16_engine('tgear', thtl), h, mach);

[CX,~,CZ,~,Cm,~] = f16_aero(alpha, 0, dh, 0, 0, 0, 0, 0, 0, V, xcg);

theta = deg2rad(alpha) + gam;                 % beta = 0, wings level

Vdot     = (qbar*c.S*CX + T)/c.mass - c.g*sin(theta);
alphadot = (qbar*c.S*CZ)/(c.mass*V) + c.g*cos(theta)/V;
qdot     = qbar*c.S*c.cbar*Cm/c.Iyy;

res = [Vdot; alphadot; qdot];

end
