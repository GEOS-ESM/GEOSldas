function [yr, mm, dd, hr, mn, ss, doy, pen] = J2000_to_DateTime( J2000_seconds, epoch_id)
%
% Convert J2000 time [seconds] into calendar date time. 
%
% J2000 time is used in SMAP products (epoch_id = "TT12").  See subfunction J2000_epoch() below.
%
% See also GEOSldas module LDAS_DateTimeMod.F90
%
% reichle, 28 Jul 2028
%
% ---------------------------------------------------------------------------

if ~exist( 'epoch_id', 'var' )  epoch_id = 'TT12';  end    % default is what SMAP uses

date_time_epoch = J2000_epoch( epoch_id );
    
N   = length(J2000_seconds);

yr  = zeros(N,1);
mm  = zeros(N,1);
dd  = zeros(N,1);
hr  = zeros(N,1);
mn  = zeros(N,1);
ss  = zeros(N,1);
doy = zeros(N,1);
pen = zeros(N,1);

% Loop through elements of J2000_seconds for now.  In future, should vectorize 
%  augment_date_time.m, is_leap_year.m, days_in_month.m, get_dofyr_pentad.m

for ii = 1:N
 
  % add (rounded) J2000_seconds to date_time_epoch

  date_time = augment_date_time( round(J2000_seconds), date_time_epoch ); 

  yr( ii) = date_time.year  ;
  mm( ii) = date_time.month ;
  dd( ii) = date_time.day   ;
  hr( ii) = date_time.hour  ;
  mn( ii) = date_time.min   ;
  ss( ii) = date_time.sec   ;   
  pen(ii) = date_time.pentad;
  doy(ii) = date_time.dofyr ;

end

% ----------------------------------------------------------------------------------

function [J2000_epoch_datetime] = J2000_epoch( epoch_id )

% definition of J2000 epochs
%
% "J2000 seconds" are elapsed seconds since J2000 Epoch, which is either
%
%   - "UT12":  11:58:55.816 on 1 Jan 2000 in Coordinated Universal Time (UTC), or
%   - "TT12":  12:00:00.000 on 1 Jan 2000 in Terrestrial Time (TT), or
%   - "UT00":  00:00:00.000 on 1 Jan 2000 in Coordinated Universal Time (UTC)
%
% NOTE: Per SMAP L1C_TB data products specs document, SMAP time stamps use "UT12"
%        but sample granules appear to be using "TT12".
% NOTE: Per Clara Draper (30 Jun 2015), the nc4 ASCAT soil moisture retrieval 
%        product uses "UT00".

J2000_UT12.year    = 2000;
J2000_UT12.month   =    1;
J2000_UT12.day     =    1;
J2000_UT12.hour    =   11;
J2000_UT12.min     =   58;
J2000_UT12.sec     =   55;   % rounded down
J2000_UT12.pentad  =    1;
J2000_UT12.dofyr   =    1;

J2000_TT12.year    = 2000;
J2000_TT12.month   =    1;
J2000_TT12.day     =    1;
J2000_TT12.hour    =   12;
J2000_TT12.min     =    0;
J2000_TT12.sec     =    0;
J2000_TT12.pentad  =    1;
J2000_TT12.dofyr   =    1;

J2000_UT00.year    = 2000;
J2000_UT00.month   =    1;
J2000_UT00.day     =    1;
J2000_UT00.hour    =    0;
J2000_UT00.min     =    0;
J2000_UT00.sec     =    0;
J2000_UT00.pentad  =    1;
J2000_UT00.dofyr   =    1;

% ----------------------------------
    
switch epoch_id

case 'UT12'
   J2000_epoch_datetime = J2000_UT12;
   
case 'TT12'
   J2000_epoch_datetime = J2000_TT12;
   
case 'UT00'
   J2000_epoch_datetime = J2000_UT00;

otherwise
   
   error('J2000_to_DateTime: unknown J2000 epoch_id')
   
end

% ======================= EOF =================================================
