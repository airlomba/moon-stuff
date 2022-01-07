% Function to find next moon rise from a preset date and time
% Emmanuel Lomba (CT7AFR)
%
% Syntax:
%   [Y, M, D, H, M, S, Az, Al] = find_next_moonrise (LAT, LON, YY, MM, DD, hh, mm, ss)
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
%   Year  - Year of next moonrise
%   Month - Month of next moonrise
%   Day   - Day of next moonrise
%   Hours - Hours of next moonrise
%   Minutes - Minutes of next moonrise
%   Seconds - Seconds of next moonrise
%   Azimuth - Azimuth of the Moon at moonrise in decimal degrees (N=0, E=90)
%   Altitude - Altitude above horizon in decimanl degrees
%
% Dependencies:
%   Function find_moon()
%

function [Year, Month, Day, Hour, Minute, Second, Azimuth, Altitude] = find_next_moonrise (LAT, LON, YY, MM, DD, hh, mm, ss)
  agora = [YY, MM, DD, hh, mm, ss];
  [Az, Alt, d] = find_moon(LAT, LON, ...
      agora(1), agora(2), agora(3), agora(4), agora(5), agora(6));
  
  if Alt < 0.0
    i = 1;
    while Alt < 0.0
      depois_secs = datenum(agora) + (i/1440);   % increment ONE minute
      depois_vect = datevec(depois_secs);
      YY = depois_vect(1);
      MM = depois_vect(2);
      DD = depois_vect(3);
      hh = depois_vect(4);
      mm = depois_vect(5);
      ss = depois_vect(6);
      [Az, Alt, d] = find_moon(LAT, LON, YY, MM, DD, hh, mm, ss);
      i = i + 1;
    end
    Year = YY; Month = MM; Day = DD;
    Hour = hh; Minute = mm; Second = ss;
    Azimuth = Az;
  else
    i = 1;
    while Alt > 0.0
      depois_secs = datenum(agora) + (i/1440);   % increment ONE minute
      depois_vect = datevec(depois_secs);
      YY = depois_vect(1);
      MM = depois_vect(2);
      DD = depois_vect(3);
      hh = depois_vect(4);
      mm = depois_vect(5);
      ss = depois_vect(6);
      [Az, Alt, d] = find_moon(LAT, LON, YY, MM, DD, hh, mm, ss);
      i = i + 1;
    end
    while Alt < 0.0
      depois_secs = datenum(agora) + (i/1440);   % increment ONE minute
      depois_vect = datevec(depois_secs);
      YY = depois_vect(1);
      MM = depois_vect(2);
      DD = depois_vect(3);
      hh = depois_vect(4);
      mm = depois_vect(5);
      ss = depois_vect(6);
      [Az, Alt, d] = find_moon(LAT, LON, YY, MM, DD, hh, mm, ss);
      i = i + 1;
    end
    Year = YY; Month = MM; Day = DD;
    Hour = hh; Minute = mm; Second = ss;
    Azimuth = Az;
    Altitude = Alt;
  end                      
end