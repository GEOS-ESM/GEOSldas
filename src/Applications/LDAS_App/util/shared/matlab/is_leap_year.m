function [ leap ] = is_leap_year(year) 

% determine whether a given year is a leap yearb
%
% input:  year, must be SCALAR !
%
% output: leap = 0 if year is not leap year
%         leap = 1 if year is leap year
%
% reichle, 1 Mar 2001
%
% ---------------------------------------------------------------------

if (length(year) ~= 1)
  disp('error, input to is_leap_year() must be scalar, exiting...')
  return
end

if (mod(year,4) ~= 0) 
  leap = 0;
elseif (mod(year,400) == 0)
  leap = 1;
elseif (mod(year,100) == 0) 
  leap = 0;
else
  leap = 1;
end 

% ========= EOF =========================================================
  