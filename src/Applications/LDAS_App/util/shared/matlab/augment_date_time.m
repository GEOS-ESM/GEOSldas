
function [date_time] = augment_date_time( dtstep, date_time_old )

% reichle, 22 Jun 2005

% dtstep in seconds

date_time = date_time_old;

dtstep_left = dtstep;

if isnan(dtstep)
  
  error('ERROR: dtstep=NaN in augment_date_time.m');
  
elseif dtstep==0   % trivial case

  date_time = get_dofyr_pentad(date_time);
  
  return
  
elseif dtstep>0
  
  while (dtstep_left>0)
    
    % increase by one day at a time
    
    dtstep_tmp  = min( dtstep_left, 86400 );
    
    dtstep_left = round( dtstep_left - dtstep_tmp );
    
    % compute secs_in_day from hh:mm:ss
  
    secs_in_day = date_time.hour*3600 + date_time.min*60 + date_time.sec;
    
    % augment 
    
    secs_in_day = secs_in_day + dtstep_tmp;
    
    % compute new hh:mm:ss from secs_in_day
    
    date_time.hour = (floor(mod(secs_in_day,86400)/3600));
    date_time.min  = (floor(mod(secs_in_day,86400)/60) - date_time.hour*60);
    date_time.sec  = (floor(mod(secs_in_day,86400))    ...
		      - date_time.hour*3600            ...
		      - date_time.min*60);
    
    % augment year/month/day and dofyr as necessary
    
    if ( secs_in_day >= 86400 ) 
      
      % get number of days in month  
      
      last_day = days_in_month( date_time.year, date_time.month);
      
      if (date_time.day==last_day) 
	
	if (date_time.month==12) 
	  
	  date_time.year  = date_time.year + 1;
	  date_time.month     = 1;
	  date_time.day       = 1;
	  
	else
	  
	  date_time.month = date_time.month + 1;
	  date_time.day   = 1;
	  
	end
	
      else
	
	date_time.day   = date_time.day   + 1;
	
      end
      
    end
    
  end
  
else
  
  while (dtstep_left<0)
    
    % decrease by one day at a time
    
    dtstep_tmp  = max( dtstep_left, -86400 );
    
    dtstep_left = round( dtstep_left - dtstep_tmp );
    
    % compute secs_in_day from hh:mm:ss
  
    secs_in_day = date_time.hour*3600 + date_time.min*60 + date_time.sec;
    
    % augment 
    
    secs_in_day = secs_in_day + dtstep_tmp;
    
    % compute new hh:mm:ss from secs_in_day
    
    secs_in_day_tmp = secs_in_day + 86400;
    
    date_time.hour = (floor(mod(secs_in_day_tmp,86400)/3600));
    date_time.min  = (floor(mod(secs_in_day_tmp,86400)/60)-date_time.hour*60);
    date_time.sec  = (floor(mod(secs_in_day_tmp,86400))    ...
		      - date_time.hour*3600            ...
		      - date_time.min*60);
    
    % augment year/month/day and dofyr as necessary
    
    if ( secs_in_day < 0 ) 
      
      if (date_time.day==1) 
	
	if (date_time.month==1) 
	  
	  date_time.year  = date_time.year - 1;
	  date_time.month     = 12;
	  date_time.day       = 31;
	  
	else
	  
	  date_time.month = date_time.month - 1;
	  
	  % get number of days in previous month  
	  
	  date_time.day = days_in_month( date_time.year, date_time.month);
	  
	end
	
      else
	
	date_time.day   = date_time.day   - 1;
	
      end
      
    end
    
  end
  
end

% get dofyr and pentad

date_time = get_dofyr_pentad(date_time);

  
% ****** EOF *******************************************************

