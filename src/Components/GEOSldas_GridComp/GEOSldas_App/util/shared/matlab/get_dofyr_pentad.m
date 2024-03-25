
function [ date_time ] = get_dofyr_pentad( date_time )

% compute dofyr and pentad for date_time

date_time.dofyr = date_time.day;

% add up days in months prior to current month

for i=1:(date_time.month-1)
  
  date_time.dofyr = date_time.dofyr + days_in_month(date_time.year,i);
  
end 

date_time.pentad = pentad_of_year(date_time.dofyr, date_time.year);


% ======================= EOF ==================================