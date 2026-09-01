function xdot = f16_dynamics(x, u, h, xcg, smooth)
%F16_DYNAMICS  Six-degree-of-freedom equations of motion for the NASA F-16.
%
%   xdot = F16_DYNAMICS(x, u, h, xcg)
%
%   State vector (8 elements), matching the formulation used by Kolb (2017)
%   for the F-18 HARV so the two studies are directly comparable:
%
%     x(1) = V       true airspeed                    [m/s]
%     x(2) = alpha   angle of attack                  [rad]
%     x(3) = beta    sideslip angle                   [rad]
%     x(4) = p       roll rate,  body axis            [rad/s]
%     x(5) = q       pitch rate, body axis            [rad/s]
%     x(6) = r       yaw rate,   body axis            [rad/s]
%     x(7) = phi     bank angle                       [rad]
%     x(8) = theta   pitch attitude                   [rad]
%
%   Control vector (4 elements):
%     u(1) = throttle    [0..1]
%     u(2) = stabilator  [deg], positive trailing edge down
%     u(3) = aileron     [deg]
%     u(4) = rudder      [deg]
%
%   Altitude is a fixed parameter rather than a state, and heading is
%   omitted because nothing in the model depends on it.  Engine spool
%   dynamics are omitted so that thrust follows the throttle directly; add a
%   ninth state with a first-order lag if you need throttle transients.
%
%   Validation, level flight at 152.4 m/s, sea level, x_cg = 0.35:
%   xdot is zero to machine precision at the trim point and the eigenvalues
%   of the Jacobian are approximately
%
%       -3.60                roll subsidence
%       -0.40 +- 2.85i       dutch roll
%       -0.019 +- 0.123i     phugoid
%       -0.016               spiral
%       +0.571, -2.671       pitch, split into two real roots
%
%   The positive real root is the relaxed static stability of the bare
%   airframe: no oscillatory short period at all, just a divergence with a
%   time to double of about 1.2 s.  This is the dynamic counterpart of the
%   unstable low-alpha equilibrium found in run_pitch_equilibria, and it is
%   why the real aircraft cannot be flown without its flight control system.

if nargin < 3 || isempty(h),      h   = 0;     end
if nargin < 4 || isempty(xcg),    xcg = 0.35;  end
if nargin < 5 || isempty(smooth), smooth = false; end

% smooth = true refines the alpha axis so that the Jacobian is continuous.
% Required for eigenvalue and continuation work; see the header of f16_aero.

c = f16_constants();

V     = max(x(1), 20);      % guard against the solver driving V through zero
alpha = x(2);  beta = x(3);
p     = x(4);  q    = x(5);  r = x(6);
phi   = x(7);  theta = x(8);

thtl = min(max(u(1),0),1);
dh   = u(2);  da = u(3);  dr = u(4);

% ---- Forces ------------------------------------------------------------
[rho, a_s] = f16_engine('atmos', h);
qbar = 0.5*rho*V^2;
T    = f16_engine('thrust', f16_engine('tgear', thtl), h, V/a_s);

[CX,CY,CZ,Cl,Cm,Cn] = f16_aero(rad2deg(alpha), rad2deg(beta), dh, da, dr, 0, ...
                               p, q, r, V, xcg, smooth);

Fx = qbar*c.S*CX + T;
Fy = qbar*c.S*CY;
Fz = qbar*c.S*CZ;

% ---- Body-axis velocity components and their derivatives ---------------
ca = cos(alpha);  sa = sin(alpha);
cb = cos(beta);   sb = sin(beta);
uu = V*ca*cb;     vv = V*sb;      ww = V*sa*cb;

ud = r*vv - q*ww - c.g*sin(theta)                 + Fx/c.mass;
vd = p*ww - r*uu + c.g*cos(theta)*sin(phi)        + Fy/c.mass;
wd = q*uu - p*vv + c.g*cos(theta)*cos(phi)        + Fz/c.mass;

Vdot     = (uu*ud + vv*vd + ww*wd)/V;
alphadot = (uu*wd - ww*ud)/(uu^2 + ww^2);
betadot  = (V*vd - vv*Vdot)/(V^2*cb);

% ---- Moments, including the Ixz cross-product and engine gyroscopics ----
L =  qbar*c.S*c.b*Cl;
M =  qbar*c.S*c.cbar*Cm - c.Heng*r;
N =  qbar*c.S*c.b*Cn    + c.Heng*q;

Gam = c.Ixx*c.Izz - c.Ixz^2;

pdot = ( c.Izz*L + c.Ixz*N ...
       - ( c.Ixz*(c.Iyy - c.Ixx - c.Izz)*p ...
         + (c.Ixz^2 + c.Izz*(c.Izz - c.Iyy))*r )*q ) / Gam;

qdot = ( M - (c.Ixx - c.Izz)*p*r - c.Ixz*(p^2 - r^2) ) / c.Iyy;

rdot = ( c.Ixz*L + c.Ixx*N ...
       + ( c.Ixz*(c.Iyy - c.Ixx - c.Izz)*r ...
         + (c.Ixz^2 + c.Ixx*(c.Ixx - c.Iyy))*p )*q ) / Gam;

% ---- Attitude kinematics ------------------------------------------------
phidot   = p + tan(theta)*(q*sin(phi) + r*cos(phi));
thetadot = q*cos(phi) - r*sin(phi);

xdot = [Vdot; alphadot; betadot; pdot; qdot; rdot; phidot; thetadot];

end