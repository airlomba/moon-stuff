% Function to calculate the position of the Moon
%
% Sytnax:
%   [Az, Alt, Dist] = find_moon (LAT, LON, YY, MM, DD, hh, mm, ss)
%
% Inputs:
%   LAT - Latitude in decimal degree (North is positive)
%   LON - Longitude in decimal degrees (East is positive)
%   YY  - Year (four digits)
%   MM  - Month number (e.g., 7 for July)
%   DD  - Day number (e.g., 19 for 19th)
%   hh  - hours (0 to 23)
%   mm  - minutes (0 to 59)
%   ss  - seconds (0 to 59)
%
% Outputs:
%   Az  - Azimuth in decimal degrees (North: 0, East: 90, South: 180, West: 270)
%   Alt - Altitude above horizon, in decimal degrees (negative = below horizon)
%   Dist- Distance from Earth coordinates to Moon, in kilometers
%

function [Az, Alt, Dist] = find_moon (LAT, LON, YY, MM, DD, hh, mm, ss)
  UT = (hh + mm / 60 + ss/3600) / 24;

  %Julian Date:
  if MM < 3
    YY = YY -1;
    MM = MM +12;
  end
  A=fix(YY/100);
  B=fix(A/4);
  C=2-A+B;
  E = fix(365.25*(YY+4716));
  F=fix(30.6001*(MM+1));
  JD=C+DD+E+F-1524.5;
  JD_t= JD + UT;

  d = 367*YY - fix((7*(YY + fix((MM+9)/12)))/4) ...
        + fix((275*MM)/9) + DD - 730530 + UT;

  %Sidereal Time
  T_t = (JD_t - 2451545.0)/36525;
  GMST = wrapTo360(280.46061837 + 360.98564736629*(JD_t - 2451545.0) + ...
              0.000387933*T_t*T_t - (T_t*T_t*T_t)/38710000);
  LST = wrapTo360(GMST + LON);

  oblecl =  wrapTo360(23.4393 - 3.563E-7 * d);   % obliquity of the ecliptic

  % Sun orbital elements:
  w_s = wrapTo360(282.9404 + 4.70935E-5 * d);    %(longitude of perihelion)
  a_s = 1.000000;                               %(mean distance, a.u.)
  e_s = 0.016709 - 1.151E-9 * d;    %(eccentricity)
  M_s = wrapTo360(356.0470 + 0.9856002585 * d);    %(mean anomaly)
  L_s = wrapTo360(w_s + M_s);   % Sun mean longitude
  E_s = M_s + (180/pi) * e_s * sin(M_s*(pi/180)) * ...
          (1 + e_s * cos(M_s*(pi/180))); % Eccentric anomaly

  x_s = cos(E_s*(pi/180)) - e_s;   %Sun's rectangular coordinates
  y_s = sin(E_s*(pi/180)) * sqrt(1 - e_s*e_s);   %Sun's rectangular coordinates
  r_s = sqrt(x_s*x_s + y_s*y_s);     % Sun's distance
  v_s = wrapTo360(atan2d(y_s, x_s)); % Sun's true anomaly
  lon_s = wrapTo360(v_s + w_s);    % Sun's longitude

  x_s = r_s * cos(lon_s*(pi/180));   % Sun's ecliptic rectangular coordinates
  y_s = r_s * sin(lon_s*(pi/180));   % Sun's ecliptic rectangular coordinates
  z_s = 0.0;                         % Sun's ecliptic rectangular coordinates

  xequat_s = x_s;    % Sun's rotated ecliptic rectangular coordinates
  yequat_s = y_s * cos(oblecl*(pi/180)) - z_s * sin(oblecl*(pi/180));
  zequat_s = y_s * sin(oblecl*(pi/180)) + z_s * cos(oblecl*(pi/180));

  r_s    =  sqrt( xequat_s*xequat_s + yequat_s*yequat_s + zequat_s*zequat_s );
  RA_s   =  atan2d( yequat_s, xequat_s ); % [deg]
  Decl_s =  atan2d( zequat_s, sqrt( xequat_s*xequat_s + yequat_s*yequat_s) );

  %Hour Angle. Altitude and Azimuth
  HA_s = wrapTo360(LST - RA_s);

  x_temp = cos(HA_s*(pi/180)) * cos(Decl_s*(pi/180));
  y_temp = sin(HA_s*(pi/180)) * cos(Decl_s*(pi/180));
  z_temp = sin(Decl_s*(pi/180));

  xhor_s = x_temp * sin(LAT*(pi/180)) - z_temp * cos(LAT*(pi/180));
  yhor_s = y_temp;
  zhor_s = x_temp * cos(LAT*(pi/180)) + z_temp * sin(LAT*(pi/180));

  SUN_azimuth  = wrapTo360(atan2d( yhor_s, xhor_s ) + 180);
  SUN_altitude = atan2d( zhor_s, sqrt(xhor_s*xhor_s+yhor_s*yhor_s));

  % Moon orbital elements:
  N_m = wrapTo360(125.1228 - 0.0529538083  * d);    %(Long asc. node)
  i_m =   5.1454;                        %(Inclination)
  w_m = wrapTo360(318.0634 + 0.1643573223  * d);    %(Arg. of perigee)
  a_m =  60.2666;                        %(Mean distance) Earth equatorial radii
  e_m = 0.054900;                        %(Eccentricity)
  M_m = wrapTo360(115.3654 + 13.0649929509 * d);    %(Mean anomaly)

  E0_m = M_m + (180/pi) * e_m * sin(M_m*(pi/180)) * ...
          (1 + e_m * cos(M_m*(pi/180)));
  E1_m = E0_m - (E0_m - (180/pi) * e_m * sin(E0_m*(pi/180)) ...
          - M_m) / (1 - e_m * cos(E0_m*(pi/180)));
  Max_Error = abs(E0_m - E1_m);
  while Max_Error > 0.005
    E0_m = E1_m;
    E1_m = E0_m - (E0_m - (180/pi) * e_m * sin(E0_m*(pi/180)) - M_m) ...
          / (1 - e_m * cos(E0_m*(pi/180)));
    Max_Error = abs(E0_m - E1_m);
  end
  E_m = E1_m;    % Eccentricyti of the Moon

  x_m = a_m * (cos(E_m*(pi/180)) - e_m);
  y_m = a_m * sqrt(1 - e_m*e_m) * sin(E_m*(pi/180));

  r_m = sqrt( x_m*x_m + y_m*y_m ); % Moon distance
  v_m = wrapTo360(atan2d(y_m,x_m));  % True anomaly

  xeclip_m = r_m * ( cos(N_m*(pi/180)) * cos((v_m+w_m)*(pi/180)) ...
          - sin(N_m*(pi/180)) * sin((v_m+w_m)*(pi/180)) * cos(i_m*(pi/180)) );
  yeclip_m = r_m * ( sin(N_m*(pi/180)) * cos((v_m+w_m)*(pi/180)) ...
          + cos(N_m*(pi/180)) * sin((v_m+w_m)*(pi/180)) * cos(i_m*(pi/180)) );
  zeclip_m = r_m * sin((v_m+w_m)*(pi/180)) * sin(i_m*(pi/180));

  long_m =  wrapTo360(atan2d( yeclip_m, xeclip_m ));
  lat_m  =  atan2d( zeclip_m, sqrt( xeclip_m*xeclip_m + yeclip_m*yeclip_m ) );
  r_m    =  sqrt( xeclip_m*xeclip_m + yeclip_m*yeclip_m + zeclip_m*zeclip_m );

  % Moon position higher accuracy
  L_m  =  wrapTo360(N_m + w_m + M_m);   % Moon's Mean Longitude
  D   =  wrapTo360(L_m - L_s);   % Moon's Mean elongation
  F   =  wrapTo360(L_m - N_m);     % Moon's argument of latitude

  % Perturbations in longitude
  Pert_lon =  -1.274 * sin((M_m - 2*D)*(pi/180)) ... % evection
              +0.658 * sin(2*D*(pi/180)) ...  % variation
              -0.186 * sin(M_s*(pi/180)) ... %yearly equation
              -0.059 * sin((2*M_m - 2*D)*(pi/180)) ...
              -0.057 * sin((M_m - 2*D + M_s)*(pi/180)) ...
              +0.053 * sin((M_m + 2*D)*(pi/180)) ...
              +0.046 * sin((2*D - M_s)*(pi/180)) ...
              +0.041 * sin((M_m - M_s)*(pi/180)) ...
              -0.035 * sin(D*(pi/180))  ...   % parallactic equation
              -0.031 * sin((M_m + M_s)*(pi/180))  ...
              -0.015 * sin((2*F - 2*D)*(pi/180))  ...
              +0.011 * sin((M_m - 4*D)*(pi/180));
  % Perturbations in Latitude
  Pert_lat =  -0.173 * sin((F - 2*D)*(pi/180)) ...
              -0.055 * sin((M_m - F - 2*D)*(pi/180)) ...
              -0.046 * sin((M_m + F - 2*D)*(pi/180)) ...
              +0.033 * sin((F + 2*D)*(pi/180)) ...
              +0.017 * sin((2*M_m + F)*(pi/180));
  % Perturbations in distance            
  Pert_r = -0.58 * cos((M_m - 2*D)*(pi/180)) ...
              -0.46 * cos((2*D)*(pi/180));
              
  long_m = long_m + Pert_lon;      % [deg]
  lat_m = lat_m + Pert_lat;  % [deg]
  r_m = r_m + Pert_r;    % [Earth radii]

  xeclip_m = r_m * cos(long_m*(pi/180)) * cos(lat_m*(pi/180));
  yeclip_m = r_m * sin(long_m*(pi/180)) * cos(lat_m*(pi/180));
  zeclip_m = r_m * sin(lat_m*(pi/180));

  xequat_m = xeclip_m;
  yequat_m = yeclip_m * cos(oblecl*(pi/180)) - zeclip_m * sin(oblecl*(pi/180));
  zequat_m = yeclip_m * sin(oblecl*(pi/180)) + zeclip_m * cos(oblecl*(pi/180));

  RA_m   = wrapTo360(atan2d( yequat_m, xequat_m ));
  Decl_m = atan2d(zequat_m, sqrt(xequat_m*xequat_m + yequat_m*yequat_m));

  % Moon topocentric position
  m_par = asind( 1/r_m );   % Moon parallax
  gclat = LAT - 0.1924 * sin(2*(LAT*pi/180));
  rho   = 0.99833 + 0.00167 * cos(2*(LAT*pi/180));
  HA_m = wrapTo360(LST - RA_m);
  g = wrapTo360(atand( tan(gclat*(pi/180)) / cos(HA_m*(pi/180))));  % aux angle
  topo_RA_m   = RA_m  - m_par * rho * cos(gclat*(pi/180)) * ...
            sin(HA_m*(pi/180)) / cos(Decl_m*(pi/180));
  topo_Decl_m = Decl_m - m_par * rho * sin(gclat*(pi/180)) * ...
            sin((g - Decl_m)*(pi/180)) / sin(g*(pi/180));

  %Hour Angle. Altitude and Azimuth
  HA_m = wrapTo360(LST - topo_RA_m);  % [degrees]

  x_temp = cos(HA_m*(pi/180)) * cos(topo_Decl_m*(pi/180));
  y_temp = sin(HA_m*(pi/180)) * cos(topo_Decl_m*(pi/180));
  z_temp = sin(topo_Decl_m*(pi/180));

  xhor_m = x_temp * sin(LAT*(pi/180)) - z_temp * cos(LAT*(pi/180));
  yhor_m = y_temp;
  zhor_m = x_temp * cos(LAT*(pi/180)) + z_temp * sin(LAT*(pi/180));

  Az  = atan2d( yhor_m, xhor_m ) + 180;
  Alt = atan2d( zhor_m, sqrt(xhor_m*xhor_m+yhor_m*yhor_m) );
  Dist = round(r_m * 6371);
end
