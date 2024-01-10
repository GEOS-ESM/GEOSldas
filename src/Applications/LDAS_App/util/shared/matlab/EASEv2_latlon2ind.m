%
%  SMAPEASE2FORWARD  The principal function is to perform forward transformation
%                    from (lat,lon)'s to (row,col)'s for a set of nested EASE
%                    grids defined at 1, 3, 9, and 36km grid resolutions.  These
%                    grids are all based on the EASE-Grid 2.0 specification (WGS84
%                    ellipsoid).
%
%  SYNTAX            [row,col] = smapease2forward(lat,lon,gridid)
%
%                    where gridid is a 3-character string enclosed in single
%                    quotes, in the form of {M|N|S}{01,03,09,36}.  This subroutine
%                    accepts vector inputs and produce vector outputs.
%
%  HISTORY           This subroutine was adapted from the offical EASE-Grid-2.0
%                    conversion utilities (written in IDL) developed by the
%                    NSIDC.
%
%                    Note that in NSIDC's original implementation, (row,col) are
%                    zero-based.  In other words, the first cell is (0,0) and the
%                    last cell is (N-1,M-1), where N and M are the row and column
%                    dimensions of the array.  In this MATLAB implementation, the
%                    same convention is used.  In other words, the end point of
%                    the first cell is located at (r,c) = (-0.5,-0.5) whereas the
%                    end point of the last cell is located at (r,c) = (14615.5,
%                    34703.5).  Thus,
%
%                   [lat,lon] = smapease2inverse(-0.5,-0.5,'M01') returns:
%                    lat = 85.044566407398861
%                    lon = 1.799999999999994e+02
%
%                   [lat,lon] = smapease2inverse(14615.5,34703.5,'M01') returns:
%                    lat = -85.044566407398861
%                    lon = -1.799999999999994e+02
%
%                    The polar grids, on the other hand, are more complete in
%                    terms of latitude coverage:
%
%                   [lat,lon] = smapease2inverse(8999,8999,'N01')
%                    lat = 89.993669248945238 
%                    lon = -135
%                   [lat,lon] = smapease2inverse(9000,9000,'N01')
%                    lat = 89.993669248945238
%                    lon = 45
%
%                   [lat,lon] = smapease2inverse(8999,8999,'S01')
%                    lat = -89.993669248945238 
%                    lon = -45
%                   [lat,lon] = smapease2inverse(9000,9000,'S01')
%                    lat = -89.993669248945238
%                    lon = 135
%
%  UPDATE            North/south polar projections were added. (03/2012)
%
%  REFERENCE         Brodzik, M. J., B. Billingsley, T. Haran, B. Raup, and M. H.
%                    Savoie (2012): EASE-Grid 2.0: Incremental but Significant
%                    Improvements for Earth-Gridded Data Sets. ISPRS International
%                    Journal of Geo-Information, vol. 1, no. 1, pp. 32-45,
%                    http://www.mdpi.com/2220-9964/1/1/32/
%  
%  Steven Chan, 11/2011
%  Email: steven.k.chan@jpl.nasa.gov

function [row,col] = EASEv2_latlon2ind(lat,lon,gridid,return_rounded)

% By design, [row, col] are real numbers, with the fractional portion indicating
%   the position of the specified [lat, lon] coordinates between adjacent grid
%   cell centers.  E.g., col=0.5 indicates that the input longitude is on the
%   boundary between grid cells associated with col=0 and col=1.
% If the [optional] input argument 'return_rounded' is present and ~=0, then 
%   [row, col] are rounded to the nearest integer.

% Constants returned by EASE2_GRID_INFO.PRO
projection = gridid(1);
switch gridid
  case 'M36'
    map_scale_m = 36032.220840584;
    cols = 964;
    rows = 406;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'M09'
    map_scale_m = 9008.055210146;
    cols = 3856;
    rows = 1624;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'M03'
    map_scale_m = 3002.6850700487;
    cols = 11568;
    rows = 4872;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'M01'
    map_scale_m = 1000.89502334956;
    cols = 34704;
    rows = 14616;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'N36'
    map_scale_m = 36000.0;
    cols = 500;
    rows = 500;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'N09'
    map_scale_m = 9000.0;
    cols = 2000;
    rows = 2000;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'N03'
    map_scale_m = 3000.0;
    cols = 6000;
    rows = 6000;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'N01'
    map_scale_m = 1000.0;
    cols = 18000;
    rows = 18000;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'S36'
    map_scale_m = 36000.0;
    cols = 500;
    rows = 500;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'S09'
    map_scale_m = 9000.0;
    cols = 2000;
    rows = 2000;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'S03'
    map_scale_m = 3000.0;
    cols = 6000;
    rows = 6000;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  case 'S01'
    map_scale_m = 1000.0;
    cols = 18000;
    rows = 18000;
    r0 = (cols-1)/2;
    s0 = (rows-1)/2;
  otherwise
    disp(['ERROR: Incompatible grid specification.']);
end

% Constants returned by EASE2_MAP_INFO.PRO
epsilon = 1.0e-6;
map_equatorial_radius_m = 6378137.0;
map_eccentricity = 0.081819190843;
e2 = map_eccentricity^2;
switch projection
  case 'M'
    map_reference_latitude = 0.0;
    map_reference_longitude = 0.0;
    map_second_reference_latitude = 30.0;
    sin_phi1 = sin(map_second_reference_latitude*pi/180);
    cos_phi1 = cos(map_second_reference_latitude*pi/180);
    kz = cos_phi1/sqrt(1.0-e2*sin_phi1*sin_phi1);
  case 'N'
    map_reference_latitude = 90.0;
    map_reference_longitude = 0.0;
  case 'S'
    map_reference_latitude = -90.0;
    map_reference_longitude = 0.0;
end

% Selected calculations inside WGS84_CONVERT.PRO and WGS84_CONVERT_XY.PRO
dlon = lon - map_reference_longitude;
msk1 = dlon < -180.0; dlon(msk1) = dlon(msk1) + 360.0;
msk2 = dlon > +180.0; dlon(msk2) = dlon(msk2) - 360.0;
phi = lat*pi/180.0;
lam = dlon*pi/180.0;
sin_phi = sin(phi);
q = (1.0-e2)*((sin_phi./(1.0-e2*sin_phi.*sin_phi))-(1.0/(2.0*map_eccentricity))*log((1.0-map_eccentricity*sin_phi)./(1.0+map_eccentricity*sin_phi)));
qp = 1.0-((1.0-e2)/(2.0*map_eccentricity)*log((1.0-map_eccentricity)/(1.0+map_eccentricity)));
switch projection
  case 'M'
    x = map_equatorial_radius_m*kz*lam;
    y = (map_equatorial_radius_m*q)/(2.0*kz);
  case 'N'
    tmp = qp - q;
    tmp(abs(tmp) < epsilon) = 0.0;
    rho = map_equatorial_radius_m*sqrt(tmp);
    x =  rho.*sin(lam);
    y = -rho.*cos(lam);
  case 'S'
    tmp = qp + q;
    tmp(abs(tmp) < epsilon) = 0.0;
    rho = map_equatorial_radius_m*sqrt(tmp);
    x =  rho.*sin(lam);
    y =  rho.*cos(lam);
end
row = s0-(y/map_scale_m);
col = r0+(x/map_scale_m);
switch projection
  case 'N'
    idx = lat < 0.0;
    row(idx) = NaN;
    col(idx) = NaN;
  case 'S'
    idx = lat > 0.0;
    row(idx) = NaN;
    col(idx) = NaN;
end

if exist('return_rounded','var')
  if return_rounded
    col=round(col);
    row=round(row);
  end
end

% ========= EOF =========================================================
