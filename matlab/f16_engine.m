function varargout = f16_engine(mode, varargin)
%F16_ENGINE  Thrust model and standard atmosphere for the NASA F-16 model.
%
%   T   = f16_engine('thrust', power, h, mach)   thrust [N]
%   P   = f16_engine('tgear',  throttle)         commanded power level [%]
%   [rho, a] = f16_engine('atmos', h)            density [kg/m^3], speed of sound [m/s]
%
%   Engine model from NASA TP-1538 (Nguyen et al., 1979), as tabulated in
%   Stevens & Lewis, "Aircraft Control and Simulation", Appendix A.  Three
%   thrust tables (idle, military, maximum afterburner) against Mach number
%   and altitude, blended by the commanded power level:
%       power <  50 : between idle and military
%       power >= 50 : between military and maximum
%
%   The throttle-to-power gearing has a break at 0.77 throttle, which is the
%   afterburner detent on the real aircraft.
%
%   >>> VERIFY THESE TABLES AGAINST YOUR COPY OF THE SOURCE BEFORE USING <<<
%   The aerodynamic tables come straight from F16aerodata.m and are exact.
%   The thrust tables below are transcribed by hand from the reference and
%   are the one part of this model that has not been checked against a
%   primary document.  Trim throttle and any thrust-limited manoeuvre result
%   depends on them.  The stability, trim alpha/dh and bifurcation results
%   are almost insensitive to them, so this does not block the main analysis.
%
%   Original tables are in pounds force; the conversion to newtons is applied
%   below so the rest of the model can stay in SI.

LB2N = 4.4482216152605;
FT2M = 0.3048;

switch lower(mode)

case 'tgear'
    throttle = min(max(varargin{1}, 0), 1);
    if throttle <= 0.77
        varargout{1} = 64.94*throttle;
    else
        varargout{1} = 217.38*throttle - 117.38;
    end

case 'thrust'
    power = varargin{1};  h = varargin{2};  mach = varargin{3};

    alt_grid  = [0 10000 20000 30000 40000 50000]*FT2M;    % [m]
    mach_grid = [0 0.2 0.4 0.6 0.8 1.0];

    % rows = Mach, columns = altitude, values in lbf
    T_idle = [ 1060   670   880  1140  1500  1860
                635   425   690  1070  1300  1660
                 60    25   345   755  1130  1525
              -1020  -710  -300   350   910  1360
              -2700 -1900 -1300  -247   600  1100
              -3600 -1400  -595  -342  -200   700];

    T_mil  = [12680  9150  6200  3950  2450  1400
              12680  9150  6313  4040  2470  1400
              12610  9312  6610  4290  2600  1560
              12640  9839  7090  4660  2840  1660
              12390 10176  7750  5320  3250  1930
              11680  9848  8050  6100  3800  2310];

    T_max  = [20000 15000 10800  7000  4000  2500
              21420 15700 11225  7323  4435  2600
              22700 16860 12250  8154  5000  2835
              24240 18910 13760  9285  5700  3215
              26070 21075 15975 11115  6860  3950
              28886 23319 18300 13484  8642  5057];

    m = min(max(mach, mach_grid(1)), mach_grid(end));
    z = min(max(h,    alt_grid(1)),  alt_grid(end));

    Ti = interpn(mach_grid, alt_grid, T_idle, m, z, 'linear')*LB2N;
    Tm = interpn(mach_grid, alt_grid, T_mil,  m, z, 'linear')*LB2N;
    Tx = interpn(mach_grid, alt_grid, T_max,  m, z, 'linear')*LB2N;

    if power < 50
        varargout{1} = Ti + (Tm - Ti)*power/50;
    else
        varargout{1} = Tm + (Tx - Tm)*(power - 50)/50;
    end

case 'atmos'
    % ISA, troposphere and lower stratosphere only.
    h = varargin{1};
    if h < 11000
        T = 288.15 - 0.0065*h;
        p = 101325*(T/288.15)^5.2559;
    else
        T = 216.65;
        p = 22632*exp(-9.80665*(h - 11000)/(287.05*216.65));
    end
    varargout{1} = p/(287.05*T);            % density
    varargout{2} = sqrt(1.4*287.05*T);      % speed of sound

otherwise
    error('f16_engine: unknown mode "%s"', mode);
end

end
