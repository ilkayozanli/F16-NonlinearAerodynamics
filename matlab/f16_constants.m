function c = f16_constants()
%F16_CONSTANTS  Reference geometry, mass and inertia for the NASA F-16 model.
%
%   All values are those of the simulation model documented in NASA TP-1538
%   (Nguyen et al., 1979) and reproduced in Stevens & Lewis, "Aircraft
%   Control and Simulation", Appendix A.  The originals are in US customary
%   units; the imperial value is given in the comment on each line so the
%   conversion can be checked against the source.
%
%   IMPORTANT: these are the *reference* values of the aerodynamic model,
%   not the physical dimensions of a real F-16.  The aerodynamic
%   coefficients in F16aerodata.m are non-dimensionalised with exactly
%   these S, b and cbar, so no other values may be used with them.

% ---- Reference geometry -------------------------------------------------
c.S    = 27.87;        % [m^2]  wing reference area      (300 ft^2)
c.b    = 9.144;        % [m]    reference span           (30 ft)
c.cbar = 3.450336;     % [m]    mean aerodynamic chord   (11.32 ft)

% ---- Mass and inertia (nominal loading, 20 500 lb) -----------------------
c.mass = 9295.44;      % [kg]         (636.94 slug)
c.Ixx  = 12874.8;      % [kg m^2]     (9 496 slug ft^2)
c.Iyy  = 75673.6;      % [kg m^2]     (55 814 slug ft^2)
c.Izz  = 85552.1;      % [kg m^2]     (63 100 slug ft^2)
c.Ixz  =  1331.4;      % [kg m^2]     (982 slug ft^2)
c.Heng =   216.9;      % [kg m^2/s]   engine angular momentum about x (160 slug ft^2/s)

% ---- Centre of gravity --------------------------------------------------
c.xcg_ref = 0.35;      % [-] reference c.g. of the aero tables, fraction of cbar
c.xcg     = 0.35;      % [-] actual c.g. used in the analysis; VARY THIS

% ---- Environment --------------------------------------------------------
c.g = 9.80665;         % [m/s^2]

% ---- Control surface limits (NASA TP-1538) ------------------------------
c.dh_lim   = 25.0;     % [deg] stabilator,  +- 25   (positive = trailing edge down)
c.da_lim   = 21.5;     % [deg] ailerons,    +- 21.5
c.dr_lim   = 30.0;     % [deg] rudder,      +- 30
c.dlef_lim = 25.0;     % [deg] leading edge flap, 0 to 25

% ---- Physical airframe data (NOT used in the aerodynamic model) ----------
% Kept only for the descriptive part of the report.  Note the corrected
% units: the original script had these as kg when they are tonnes.
c.length_physical  = 15.03;   % [m]
c.span_physical    =  9.45;   % [m]  span over wingtip launchers
c.height_physical  =  5.09;   % [m]
c.empty_mass       = 9207;    % [kg]  ~ 9.207 t
c.max_togw_mass    = 21772;   % [kg]  ~ 21.772 t
c.design_load      = 9;       % [g]
c.service_ceiling  = 15000;   % [m]

end
