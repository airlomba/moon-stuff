% Julian Date calculation with time
% for GNU Octave (or Matlab)
%
% Either let the script use actual (GMT)time or edit YY, MM, DD, hh, mm and ss
% with your prefered date and time.

time_struct = gmtime(time());         
%time_struct = localtime(time());   %uncomment this line for actual local time
YY = year(date);
MM = month(date);
DD = day(date);
hh = time_struct.hour;
mm = time_struct.min;
ss = time_struct.sec;

UT = (hh + mm/60 + ss/3600)/24;
  
if MM < 3
  YY = YY -1;
  MM = MM +12;
end

A=fix(YY/100);
B=fix(A/4);
C=2-A+B;
E = fix(365.25*(YY+4716));
F=fix(30.6001*(MM+1));
JD=C+DD+E+F-1524.5 + UT;

printf("Date and Time: %s\n", strftime("%d/%m/%Y - %T", time_struct));
printf("Julian Date:   %f\n\n", JD);


