%=========================================================================

function [dist_x_deg, dist_y_deg] = dist_km2deg( dist_km, lat )

    MAPL_PI     = 3.14159265358979323846;
    MAPL_RADIUS = 6371.0E3;

    % distance between latitudes is equal (always full meridional radius circle)
    % assumuming the Earth is a perfect ball.
    % NOTE: MAPL_radius (Earth radius) is in [m] and dist_km is in [km]

    dist_y_deg = dist_km .* (180./MAPL_PI) ./ (MAPL_RADIUS./1000.);

    % distance between longitudes decreases towards the poles
    % (radius of parallel circles decreases)
    % NOTE: cos() needs argument in [rad], lat is in [deg] (-90:90)

    dist_x_deg = dist_y_deg ./ cos( MAPL_PI./180. .* lat );

    if (any(dist_x_deg<0. | dist_y_deg<0.))  
         disp( 'encountered negative distance' );
    end
    
end

%=========================================================================
