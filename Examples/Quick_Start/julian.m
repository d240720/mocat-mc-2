function jd = julian(yr, mo, dy, hr, mn, sc)
% Convert calendar date to Julian Date
% julian(year, month, day, hour, minute, second)
t = datetime(yr, mo, dy, hr, mn, sc);
jd = juliandate(t);
end
