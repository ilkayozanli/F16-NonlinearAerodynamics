function [dh, alpha, gamma_deg, resid] = f16_trim_glide(V, thtl, h, xcg, guess)
%F16_TRIM_GLIDE  Steady flight at fixed throttle with the flight path free.
%
%   [dh, alpha, gamma_deg, resid] = F16_TRIM_GLIDE(V, thtl, h, xcg, guess)
%
%   Solves the same three longitudinal equilibrium equations as f16_trim,
%   but swaps which quantities are known:
%
%       f16_trim        given V and gamma, solve for throttle, dh, alpha
%       f16_trim_glide  given V and throttle, solve for dh, alpha, gamma
%
%   WHY THIS EXISTS.  Level-flight trim cannot reach high angle of attack.
%   Holding altitude at high alpha means flying slower than the aircraft can
%   sustain, so the reachable range at sea level stops near alpha = 7 deg.
%   That is useless for a high-alpha study.
%
%   Letting the flight path angle float fixes it.  The aircraft descends,
%   gravity supplies the energy the engine is not, and the equilibrium
%   family continues to alpha = 40 deg and beyond.  These are genuine
%   equilibria of the full eight-state system: all rates are zero, sideslip
%   and bank are zero, and the aircraft descends along a straight path at
%   constant speed.  This is the same family Kolb (2017) continues for the
%   F-18, where the throttle is held and the elevator varied.
%
%   Typical results at idle throttle, sea level, x_cg = 0.35:
%       V = 150 m/s  ->  alpha =  2.1 deg,  gamma =  -5.8 deg
%       V =  80 m/s  ->  alpha = 10.5 deg,  gamma =  -8.6 deg
%       V =  50 m/s  ->  alpha = 40.2 deg,  gamma = -36.8 deg

if nargin < 2 || isempty(thtl),  thtl = 0;    end
if nargin < 3 || isempty(h),     h    = 0;    end
if nargin < 4 || isempty(xcg),   xcg  = 0.35; end
if nargin < 5 || isempty(guess), guess = [0 5 0]; end   % [dh alpha gamma_deg]

c   = f16_constants();
fun = @(v) glide_residual(v, V, thtl, h, xcg, c);

if exist('fsolve','file') == 2
    opts = optimoptions('fsolve','Display','off', ...
                        'FunctionTolerance',1e-12,'StepTolerance',1e-12);
    v = fsolve(fun, guess, opts);
else
    v = fminsearch(@(vv) sum(fun(vv).^2), guess, ...
        optimset('Display','off','TolX',1e-10,'TolFun',1e-14, ...
                 'MaxFunEvals',5e4,'MaxIter',5e4));
end

dh = v(1);  alpha = v(2);  gamma_deg = v(3);
resid = max(abs(fun(v)));

end

% ------------------------------------------------------------------------
function res = glide_residual(v, V, thtl, h, xcg, c)

dh = v(1);  alpha = v(2);  gam = deg2rad(v(3));

[rho, a_s] = f16_engine('atmos', h);
qbar = 0.5*rho*V^2;
T    = f16_engine('thrust', f16_engine('tgear', min(max(thtl,0),1)), h, V/a_s);

[CX,~,CZ,~,Cm,~] = f16_aero(alpha, 0, dh, 0, 0, 0, 0, 0, 0, V, xcg);

theta = deg2rad(alpha) + gam;

Vdot     = (qbar*c.S*CX + T)/c.mass - c.g*sin(theta);
alphadot = (qbar*c.S*CZ)/(c.mass*V) + c.g*cos(theta)/V;
qdot     = qbar*c.S*c.cbar*Cm/c.Iyy;

res = [Vdot; alphadot; qdot];

end
