function [CX,CY,CZ,Cl,Cm,Cn] = f16_aero(alpha,beta,dh,da,dr,dlef,p,q,r,V,xcg,smooth)
%F16_AERO  Complete nonlinear aerodynamic coefficients of the NASA F-16 model.
%
%   [CX,CY,CZ,Cl,Cm,Cn] = F16_AERO(alpha,beta,dh,da,dr,dlef,p,q,r,V,xcg,smooth)
%
%   Inputs
%     alpha, beta   angle of attack and sideslip           [deg]
%     dh            stabilator deflection, +TE down        [deg], |dh| <= 25
%     da            aileron deflection                     [deg], |da| <= 21.5
%     dr            rudder deflection                      [deg], |dr| <= 30
%     dlef          leading edge flap deflection           [deg], 0 .. 25
%     p, q, r       body angular rates                     [rad/s]
%     V             true airspeed                          [m/s]
%     xcg           c.g. position as a fraction of cbar    [-]
%     smooth        false (default) or true, see below
%
%   Outputs: body-axis force and moment coefficients, moments taken about xcg.
%
%   THE SMOOTH OPTION, AND WHY IT EXISTS.
%
%   With smooth = false the alpha dependence is interpolated linearly on the
%   raw 5-degree table grid.  That is what the source model specifies and it
%   gives the correct coefficient VALUES, so it is right for plotting
%   coefficients, for root-finding on Cm, and for time simulation.
%
%   But linear interpolation makes the SLOPE piecewise constant, jumping at
%   every table node.  Anything that differentiates the model inherits those
%   jumps.  Linearising about a trim point at alpha = 9.9 deg and again at
%   10.1 deg gives dCm/dalpha of 0.00122 and 0.00060, a factor of two apart,
%   purely from crossing the node at 10 deg.  The eigenvalues jump with it.
%   Nothing physical happens at 10 degrees; the interpolation is jumping.
%
%   With smooth = true the alpha axis is first refined to 0.25 deg using
%   pchip, and linear interpolation is then applied on the refined grid.
%   The values are essentially unchanged, but the slope becomes effectively
%   continuous, which is what eigenvalue analysis and continuation need.
%
%   Use smooth = true for: linearisation, eigenvalues, continuation.
%   Use smooth = false for: coefficient plots, Cm root-finding, simulation.
%
%   Coefficient buildup follows Sonneveldt, "Nonlinear F-16 Model
%   Description", TU Delft 2006, Sec. 2, the companion document to the
%   F16aerodata.m tables.  The same equations appear in Russell, "Non-linear
%   F-16 Simulation using Simulink and Matlab", Univ. of Minnesota 2003.
%   Both restate NASA TP-1538.
%
%   Three structural points the NASA formulation requires and that are easy
%   to miss:
%
%   1. The "_lef" tables are the TOTAL coefficient with the leading edge flap
%      fully deployed, not an increment.  The increment is formed as
%      C_lef(alpha,beta) - C(alpha,beta,dh=0) and scaled by (1 - dlef/25),
%      so dlef = 0 gives the clean configuration and dlef = 25 the deployed one.
%
%   2. The tabulated moments are referenced to xcg_ref = 0.35*cbar.  Moving
%      the c.g. requires CZ*(xcg_ref-xcg) in Cm and -CY*(xcg_ref-xcg)*(cbar/b)
%      in Cn.  Omitting these makes every stability result valid only at the
%      reference c.g., and c.g. position is what controls whether a
%      deep-stall trim exists.
%
%   3. dCm_ds is the deep-stall pitching moment increment.  Zero below
%      alpha = 35 deg, large above it.  Without it the deep-stall equilibrium
%      branch is misplaced.
%
%   PERFORMANCE.  Interpolants are built once as griddedInterpolant objects
%   and cached, separately for each smoothing mode.  Call CLEAR FUNCTIONS if
%   you edit F16aerodata.m, otherwise the cached tables stay stale.

persistent Graw Gsm
if nargin < 12 || isempty(smooth), smooth = true; end

if smooth
    if isempty(Gsm),  Gsm  = build_interpolants(0.25); end
    G = Gsm;
else
    if isempty(Graw), Graw = build_interpolants([]);   end
    G = Graw;
end

if nargin < 11, xcg  = 0.35; end
if nargin < 10, V    = 150;  end
if nargin <  9, r    = 0;    end
if nargin <  8, q    = 0;    end
if nargin <  7, p    = 0;    end
if nargin <  6, dlef = 0;    end
if nargin <  5, dr   = 0;    end
if nargin <  4, da   = 0;    end

xcgr = 0.35;                    % reference c.g. of the tables
cbar = 3.450336;  b = 9.144;    % must match f16_constants

% Clamp to table limits.  griddedInterpolant would extrapolate otherwise, and
% the high-alpha branches of a sweep routinely ask for alpha beyond 90 deg.
a  = min(max(alpha, G.a1lo), G.a1hi);   % coarse alpha grid,  -20 .. 90
aL = min(max(alpha, G.a2lo), G.a2hi);   % LEF alpha grid,     -20 .. 45
bt = min(max(beta,  G.blo),  G.bhi);
e1 = min(max(dh,    G.e1lo), G.e1hi);   % -25 .. 25, five points
e2 = min(max(dh,    G.e2lo), G.e2hi);   % -25, 0, 25  (Cn and Cl)
e3 = min(max(dh,    G.e3lo), G.e3hi);   % deep stall grid
kl = 1 - dlef/25;                       % LEF blending factor

% ---- Axial force --------------------------------------------------------
CX0    = G.CX(a,bt,e1);
CX0_e0 = G.CX(a,bt,0);
CX = CX0 + (G.CX_lef(aL,bt) - CX0_e0)*kl ...
   + (cbar/(2*V))*q*( G.CXq(a) + G.dCXq_lef(aL)*kl );

% ---- Normal force -------------------------------------------------------
CZ0    = G.CZ(a,bt,e1);
CZ0_e0 = G.CZ(a,bt,0);
CZ = CZ0 + (G.CZ_lef(aL,bt) - CZ0_e0)*kl ...
   + (cbar/(2*V))*q*( G.CZq(a) + G.dCZq_lef(aL)*kl );

% ---- Pitching moment ----------------------------------------------------
Cm0    = G.Cm(a,bt,e1);
Cm0_e0 = G.Cm(a,bt,0);
Cm = Cm0 ...
   + CZ*(xcgr - xcg) ...                                  % c.g. transfer
   + (G.Cm_lef(aL,bt) - Cm0_e0)*kl ...
   + (cbar/(2*V))*q*( G.Cmq(a) + G.dCmq_lef(aL)*kl ) ...
   + G.dCm(a) ...                                         % trim increment
   + G.dCm_ds(a,e3);                                      % deep stall

% ---- Side force ---------------------------------------------------------
CY0        = G.CY(a,bt);
CY_lef     = G.CY_lef(aL,bt);
dCY_a20    = G.CY_da20(a,bt) - CY0;
dCY_a20lef = G.CY_da20lef(aL,bt) - CY_lef - dCY_a20;
CY = CY0 + (CY_lef - CY0)*kl ...
   + (dCY_a20 + dCY_a20lef*kl)*(da/20) ...
   + (G.CY_dr30(a,bt) - CY0)*(dr/30) ...
   + (b/(2*V))*( r*( G.CYr(a) + G.dCYr_lef(aL)*kl ) ...
               + p*( G.CYp(a) + G.dCYp_lef(aL)*kl ) );

% ---- Yawing moment ------------------------------------------------------
Cn0        = G.Cn(a,bt,e2);
Cn0_e0     = G.Cn(a,bt,0);
Cn_lef     = G.Cn_lef(aL,bt);
dCn_a20    = G.Cn_da20(a,bt) - Cn0_e0;
dCn_a20lef = G.Cn_da20lef(aL,bt) - Cn_lef - dCn_a20;
Cn = Cn0 + (Cn_lef - Cn0_e0)*kl ...
   - CY*(xcgr - xcg)*(cbar/b) ...                         % c.g. transfer
   + (dCn_a20 + dCn_a20lef*kl)*(da/20) ...
   + (G.Cn_dr30(a,bt) - Cn0_e0)*(dr/30) ...
   + (b/(2*V))*( r*( G.Cnr(a) + G.dCnr_lef(aL)*kl ) ...
               + p*( G.Cnp(a) + G.dCnp_lef(aL)*kl ) ) ...
   + G.dCnbeta(a)*bt;

% ---- Rolling moment -----------------------------------------------------
Cl0        = G.Cl(a,bt,e2);
Cl0_e0     = G.Cl(a,bt,0);
Cl_lef     = G.Cl_lef(aL,bt);
dCl_a20    = G.Cl_da20(a,bt) - Cl0_e0;
dCl_a20lef = G.Cl_da20lef(aL,bt) - Cl_lef - dCl_a20;
Cl = Cl0 + (Cl_lef - Cl0_e0)*kl ...
   + (dCl_a20 + dCl_a20lef*kl)*(da/20) ...
   + (G.Cl_dr30(a,bt) - Cl0_e0)*(dr/30) ...
   + (b/(2*V))*( r*( G.Clr(a) + G.dClr_lef(aL)*kl ) ...
               + p*( G.Clp(a) + G.dClp_lef(aL)*kl ) ) ...
   + G.dClbeta(a)*bt;

end

% ------------------------------------------------------------------------
function G = build_interpolants(dalpha)
%BUILD_INTERPOLANTS  Read F16aerodata.m once and wrap every table.
%   dalpha empty  -> raw table grid, linear interpolation (source model)
%   dalpha scalar -> alpha axis refined to that spacing with pchip first

run('F16aerodata.m');          % defines f16data in this workspace
D = f16data;

a1 = D.alpha1(:);  a2 = D.alpha2(:);  be = D.beta(:);
e1 = D.de1(:);     e2 = D.de2(:);     e3 = D.de3(:);

if ~isempty(dalpha)
    a1n = (a1(1):dalpha:a1(end))';
    a2n = (a2(1):dalpha:a2(end))';

    for nm = {'CX','CZ','Cm','Cn','Cl'}
        D.(nm{1}) = refine3(D.(nm{1}), a1, a1n);
    end
    for nm = {'CY','CY_da20','CY_dr30','Cn_da20','Cn_dr30','Cl_da20','Cl_dr30'}
        D.(nm{1}) = refine2(D.(nm{1}), a1, a1n);
    end
    for nm = {'CX_lef','CZ_lef','Cm_lef','CY_lef','Cn_lef','Cl_lef', ...
              'CY_da20lef','Cn_da20lef','Cl_da20lef'}
        D.(nm{1}) = refine2(D.(nm{1}), a2, a2n);
    end
    D.dCm_ds = refine2(D.dCm_ds, a1, a1n);
    for nm = {'CXq','CZq','Cmq','CYr','CYp','Cnr','Cnp','Clr','Clp', ...
              'dCm','dCnbeta','dClbeta'}
        D.(nm{1}) = interp1(a1, D.(nm{1})(:), a1n, 'pchip');
    end
    for nm = {'dCXq_lef','dCZq_lef','dCmq_lef','dCYr_lef','dCYp_lef', ...
              'dCnr_lef','dCnp_lef','dClr_lef','dClp_lef'}
        D.(nm{1}) = interp1(a2, D.(nm{1})(:), a2n, 'pchip');
    end

    a1 = a1n;  a2 = a2n;
end

G.a1lo = a1(1);  G.a1hi = a1(end);
G.a2lo = a2(1);  G.a2hi = a2(end);
G.blo  = be(1);  G.bhi  = be(end);
G.e1lo = e1(1);  G.e1hi = e1(end);
G.e2lo = e2(1);  G.e2hi = e2(end);
G.e3lo = e3(1);  G.e3hi = e3(end);

I3 = @(T,x,y,z) griddedInterpolant({x,y,z}, T, 'linear', 'nearest');
I2 = @(T,x,y)   griddedInterpolant({x,y},   T, 'linear', 'nearest');
I1 = @(T,x)     griddedInterpolant(x, T(:), 'linear', 'nearest');

% three-dimensional tables
G.CX = I3(D.CX, a1, be, e1);
G.CZ = I3(D.CZ, a1, be, e1);
G.Cm = I3(D.Cm, a1, be, e1);
G.Cn = I3(D.Cn, a1, be, e2);
G.Cl = I3(D.Cl, a1, be, e2);

% two-dimensional tables on the coarse alpha grid
for nm = {'CY','CY_da20','CY_dr30','Cn_da20','Cn_dr30','Cl_da20','Cl_dr30'}
    G.(nm{1}) = I2(D.(nm{1}), a1, be);
end

% two-dimensional tables on the LEF alpha grid
for nm = {'CX_lef','CZ_lef','Cm_lef','CY_lef','Cn_lef','Cl_lef', ...
          'CY_da20lef','Cn_da20lef','Cl_da20lef'}
    G.(nm{1}) = I2(D.(nm{1}), a2, be);
end

% deep stall table, alpha against the seven-point stabilator grid
G.dCm_ds = I2(D.dCm_ds, a1, e3);

% one-dimensional tables on the coarse alpha grid
for nm = {'CXq','CZq','Cmq','CYr','CYp','Cnr','Cnp','Clr','Clp', ...
          'dCm','dCnbeta','dClbeta'}
    G.(nm{1}) = I1(D.(nm{1}), a1);
end

% one-dimensional increment tables on the LEF alpha grid
for nm = {'dCXq_lef','dCZq_lef','dCmq_lef','dCYr_lef','dCYp_lef', ...
          'dCnr_lef','dCnp_lef','dClr_lef','dClp_lef'}
    G.(nm{1}) = I1(D.(nm{1}), a2);
end

end

% ------------------------------------------------------------------------
function Tn = refine3(T, aold, anew)
%REFINE3  Pchip-refine the first dimension of a 3-D table.
[~, nb, nd] = size(T);
Tn = zeros(numel(anew), nb, nd);
for j = 1:nb
    for k = 1:nd
        Tn(:,j,k) = interp1(aold, T(:,j,k), anew, 'pchip');
    end
end
end

function Tn = refine2(T, aold, anew)
%REFINE2  Pchip-refine the first dimension of a 2-D table.
nb = size(T,2);
Tn = zeros(numel(anew), nb);
for j = 1:nb
    Tn(:,j) = interp1(aold, T(:,j), anew, 'pchip');
end
end