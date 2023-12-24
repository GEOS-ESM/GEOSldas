
function pentad = pentad_of_year(day_of_year, year)

if (is_leap_year(year) & day_of_year>=59) 
  
  pentad = floor((day_of_year-2)/5)+1;
  
else
  
  pentad = floor((day_of_year-1)/5)+1;
  
end

% ======================= EOF ==================================
