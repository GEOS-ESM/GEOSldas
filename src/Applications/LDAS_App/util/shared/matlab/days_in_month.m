
function [n_days] = days_in_month( year, month )

% reichle, 22 Jun 2005

days_in_month_leap    = [31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

days_in_month_nonleap = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31];

if (is_leap_year(year)) 
  n_days = days_in_month_leap(month);
else
  n_days = days_in_month_nonleap(month);
end 

% **************************** EOF ********************************