PROGRAM tile_bin2nc4

  implicit none
  INCLUDE 'netcdf.inc'

  integer       :: i,k,  n, NTILES
  integer       :: NCFOutID, Vid, STATUS, CellID, TimID, nVars
  character(256) :: Usage="tile_bin2nc4.x BINFILE DESCRIPTOR TILECOORD"
  character(512) :: BINFILE, TILECOORD, DESCRIPTOR, arg(3)
  character(128) :: MYNAME, BUF
  integer, dimension(8)               :: date_time_values
  character (22)                      :: time_stamp
  real, allocatable, dimension (:)    :: lons, lats, var
  integer, allocatable, dimension (:) :: tileid, i_index, j_index 
  integer       :: myunit1, myunit2
  real :: undef
  ! processing command line agruments

  I = command_argument_count()
  
  if( I /=3 ) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     stop
  end if

  do n=1,I
     call get_command_argument(n,arg(n))
  enddo

  call get_environment_variable ("MYNAME"        ,MYNAME        )
  read(arg(1),'(a)') BINFILE
  read(arg(2),'(a)') DESCRIPTOR
  read(arg(3),'(a)') TILECOORD

!  print *,MYNAME
!  print *,trim(BINFILE)
!  print *,trim(DESCRIPTOR)
!  print *,trim(TILECOORD)

 ! reading TILECOORD

  open (newunit=myunit1, file = trim(TILECOORD), form = 'unformatted', action ='read')
  read (myunit1) NTILES
  allocate (lons   (1:NTILES))
  allocate (lats   (1:NTILES))
  allocate (tileid (1:NTILES))
  allocate (var    (1:NTILES))
  allocate (i_index(1:NTILES))
  allocate (j_index(1:NTILES))

  read (myunit1) tileid
  read (myunit1) tileid
  read (myunit1) tileid
  read (myunit1) lons
  read (myunit1) lats
  read (myunit1) var
  read (myunit1) var
  read (myunit1) var
  read (myunit1) var
  read (myunit1) i_index
  read (myunit1) j_index  

  close (myunit1,status = 'keep')

  ! read binary and write NC4

  open (newunit=myunit1, file = trim(DESCRIPTOR), form ='formatted', action = 'read')
  nVars = 0

  undef = 0.100000E+16
  k = 0
  do
     read(myunit1, '(a)', iostat=status) buf
     if (status /= 0) exit
     k = k + 1
     if(buf(1:index(buf,' ') -1) == 'vars') then
        i =  index(buf,' ') 
        read (buf(i:),*, IOSTAT = n) nVars
     endif
     if(buf(1:index(buf,' ') -1) == 'undef') then
        i =  index(buf,' ') 
        read (buf(i:),*, IOSTAT = n) undef
     endif
     if(nVars /= 0) exit

  end do

  status = NF_CREATE (trim(BINFILE)//'.nc4', NF_NETCDF4, NCFOutID)
  status = NF_DEF_DIM(NCFOutID, 'tile' , NTILES, CellID)
  status = NF_DEF_DIM(NCFOutID, 'time' , NF_UNLIMITED, TimID)

  status = NF_DEF_VAR(NCFOutID, 'lon'    , NF_FLOAT, 1 ,CellID, vid)
  status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',    &
         LEN_TRIM('longitude'), 'longitude') 
  status = NF_DEF_VAR(NCFOutID, 'lat'    , NF_FLOAT, 1 ,CellID, vid)
  status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',    &
         LEN_TRIM('latitude'), 'latitude') 
  status = NF_DEF_VAR(NCFOutID, 'IG'    , NF_INT, 1 ,CellID, vid)
  status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',    &
         LEN_TRIM('I_INDEX'), 'I_INDEX') 
  status = NF_DEF_VAR(NCFOutID, 'JG'    , NF_INT, 1 ,CellID, vid)
  status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',    &
         LEN_TRIM('J_INDEX'), 'J_INDEX') 
  do n = 1, nVars

     read(myunit1, '(a)', iostat=status) buf
     status = NF_DEF_VAR(NCFOutID,buf(1:index(buf,' ') -1) , NF_FLOAT, 2 ,(/CellID, TimID/), vid)
     status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',               &
          LEN_TRIM(getAttribute(buf(1:index(buf,' ') -1), LNAME = 1)),  &
          getAttribute(buf(1:index(buf,' ') -1), LNAME = 1)) 
     status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units',                   &
          LEN_TRIM(getAttribute(buf(1:index(buf,' ') -1), UNT = 1)),    &
          getAttribute(buf(1:index(buf,' ') -1), UNT = 1))     
     status = nf_put_att_real(NCFOutID, vid, '_FillValue',NF_FLOAT, 1, undef) 
  end do


  call date_and_time(VALUES=date_time_values)
          
  write (time_stamp,'(i4.4,a1,i2.2,a1,i2.2,1x,a2,1x,i2.2,a1,i2.2,a1,i2.2)')      &
       date_time_values(1),'-',date_time_values(2),'-',date_time_values(3),'at', &
       date_time_values(5),':',date_time_values(6),':',date_time_values(7)
  
  status = NF_PUT_ATT_TEXT(NCFOutID, NF_GLOBAL, 'CreatedBy', LEN_TRIM(MYNAME),  &
       trim(MYNAME))
  status = NF_PUT_ATT_TEXT(NCFOutID, NF_GLOBAL, 'Date'   , LEN_TRIM(time_stamp),trim(time_stamp))

  status = NF_ENDDEF(NCFOutID)  

  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,'lon'   ) ,(/1/),(/NTILES/),lons    )
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,'lat'   ) ,(/1/),(/NTILES/),lats    )
  status = NF_PUT_VARA_INT (NCFOutID,VarID(NCFOutID,'IG'    ) ,(/1/),(/NTILES/),i_index )
  status = NF_PUT_VARA_INT (NCFOutID,VarID(NCFOutID,'JG'    ) ,(/1/),(/NTILES/),j_index )

  ! reading and writing

  open (newunit=myunit2, file = trim(BINFILE)//'.bin', form = 'unformatted', action = 'read')

  rewind (myunit1)
  do i = 1, k
     read(myunit1, '(a)', iostat=status) buf
  end do

  do n = 1, nVars
     read (myunit1, '(a)', iostat=status) buf
     read (myunit2) var
     status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,buf(1:index(buf,' ') -1)) ,(/1,1/),(/NTILES,1/),var )
  end do  
  
  STATUS   = NF_CLOSE (NCFOutID)
  close (myunit1)
  close (myunit2)

  contains

   ! ----------------------------------------------------------------------

   integer function VarID (NCFID, VNAME) 
     
     integer, intent (in)      :: NCFID
     character(*), intent (in) :: VNAME
     integer                   :: status

     STATUS = NF_INQ_VARID (NCFID, trim(VNAME) ,VarID)
     IF (STATUS .NE. NF_NOERR) &
          CALL HANDLE_ERR(STATUS, trim(VNAME))  
     
   end function VarID

  ! -----------------------------------------------------------------------

   SUBROUTINE HANDLE_ERR(STATUS, Line)

     INTEGER,      INTENT (IN) :: STATUS
     CHARACTER(*), INTENT (IN) :: Line

     IF (STATUS .NE. NF_NOERR) THEN
        PRINT *, trim(Line),': ',NF_STRERROR(STATUS)
        STOP 'Stopped'
     ENDIF

   END SUBROUTINE HANDLE_ERR
 
 ! ***********************************************************************

  FUNCTION getAttribute (SHORT_NAME, LNAME, UNT) result (str_atr)

    character(*), intent(in)           :: SHORT_NAME
    integer, intent (in), optional     :: LNAME, UNT
    character(128)                     :: str_atr, LONG_NAME, UNITS

    SELECT case (trim(SHORT_NAME))

    ! For SM_L4 
    ! reichle, 20 May 2020: verified SHORT_NAME and corrected UNITS to match SMAP L4_SM Product Specs;  LONG_NAME (mostly) from GEOS_CatchGridComp.F90 
    ! reichle, 14 Feb 2022: added "WATERTABLED" (now: "PEATCLSM_WATERLEVEL") and "FSWCHANGE" (now: "PEATCLSM_FSWCHANGE")
    ! reichle, 21 Feb 2022: added "mwrtm_vegopacity"   
       
    case ('sm_surface');                       LONG_NAME = 'water_surface_layer';                                              UNITS = 'm3 m-3'
    case ('sm_rootzone');                      LONG_NAME = 'water_root_zone';                                                  UNITS = 'm3 m-3'
    case ('sm_profile');                       LONG_NAME = 'water_ave_prof';                                                   UNITS = 'm3 m-3'
    case ('sm_surface_wetness');               LONG_NAME = 'surface_soil_wetness';                                             UNITS = '1'
    case ('sm_rootzone_wetness');              LONG_NAME = 'root_zone_soil_wetness';                                           UNITS = '1'
    case ('sm_profile_wetness');               LONG_NAME = 'ave_prof_soil_wetness';                                            UNITS = '1'
    case ('surface_temp');                     LONG_NAME = 'ave_catchment_temp_incl_snw';                                      UNITS = 'K' 
    case ('soil_temp_layer1');                 LONG_NAME = 'soil_temperatures_layer_1';                                        UNITS = 'K' 
    case ('soil_temp_layer2');                 LONG_NAME = 'soil_temperatures_layer_2';                                        UNITS = 'K' 
    case ('soil_temp_layer3');                 LONG_NAME = 'soil_temperatures_layer_3';                                        UNITS = 'K' 
    case ('soil_temp_layer4');                 LONG_NAME = 'soil_temperatures_layer_4';                                        UNITS = 'K' 
    case ('soil_temp_layer5');                 LONG_NAME = 'soil_temperatures_layer_5';                                        UNITS = 'K' 
    case ('soil_temp_layer6');                 LONG_NAME = 'soil_temperatures_layer_6';                                        UNITS = 'K' 
    case ('snow_mass');                        LONG_NAME = 'snow_mass';                                                        UNITS = 'kg m-2'
    case ('snow_depth');                       LONG_NAME = 'snow_depth_in_snow_covered_area';                                  UNITS = 'm'
    case ('land_evapotranspiration_flux');     LONG_NAME = 'Evaporation_land';                                                 UNITS = 'kg m-2 s-1'
    case ('overland_runoff_flux');             LONG_NAME = 'runoff_flux';                                                      UNITS = 'kg m-2 s-1'
    case ('baseflow_flux');                    LONG_NAME = 'baseflow_flux';                                                    UNITS = 'kg m-2 s-1'
    case ('snow_melt_flux');                   LONG_NAME = 'Snowmelt_flux_land';                                               UNITS = 'kg m-2 s-1'
    case ('soil_water_infiltration_flux');     LONG_NAME = 'rainwater_infiltration_flux';                                      UNITS = 'kg m-2 s-1'
    case ('land_fraction_saturated');          LONG_NAME = 'fractional_area_of_saturated_zone';                                UNITS = '1'
    case ('land_fraction_unsaturated');        LONG_NAME = 'fractional_area_of_unsaturated_zone';                              UNITS = '1'
    case ('land_fraction_wilting');            LONG_NAME = 'fractional_area_of_wilting_zone';                                  UNITS = '1'
    case ('land_fraction_snow_covered');       LONG_NAME = 'fractional_area_of_land_snowcover';                                UNITS = '1'
    case ('heat_flux_sensible');               LONG_NAME = 'Sensible_heat_flux_land';                                          UNITS = 'W m-2'
    case ('heat_flux_latent');                 LONG_NAME = 'Latent_heat_flux_land';                                            UNITS = 'W m-2'
    case ('heat_flux_ground');                 LONG_NAME = 'Ground_heating_land';                                              UNITS = 'W m-2'
    case ('net_downward_shortwave_flux');      LONG_NAME = 'Net_shortwave_land';                                               UNITS = 'W m-2'
    case ('net_downward_longwave_flux');       LONG_NAME = 'Net_longwave_land';                                                UNITS = 'W m-2'
    case ('radiation_shortwave_downward_flux');LONG_NAME = 'Incident_shortwave_land';                                          UNITS = 'W m-2'
    case ('radiation_longwave_absorbed_flux'); LONG_NAME = 'surface_absorbed_longwave_flux';                                   UNITS = 'W m-2'
    case ('precipitation_total_surface_flux'); LONG_NAME = 'RainfSnowf';                                                       UNITS = 'kg m-2 s-1'
    case ('snowfall_surface_flux');            LONG_NAME = 'snowfall';                                                         UNITS = 'kg m-2 s-1'
    case ('surface_pressure');                 LONG_NAME = 'surface_pressure';                                                 UNITS = 'Pa'
    case ('height_lowatmmodlay');              LONG_NAME = 'reference_height_for_Tair_Qair_Wind';                              UNITS = 'm'
    case ('temp_lowatmmodlay');                LONG_NAME = 'air_temperature_at_RefH';                                          UNITS = 'K' 
    case ('specific_humidity_lowatmmodlay');   LONG_NAME = 'specific_humidity_at_RefH';                                        UNITS = 'kg kg-1'
    case ('windspeed_lowatmmodlay');           LONG_NAME = 'wind_speed_at_RefH';                                               UNITS = 'm s-1'
    case ('vegetation_greenness_fraction');    LONG_NAME = 'greeness_fraction';                                                UNITS = '1' 
    case ('leaf_area_index');                  LONG_NAME = 'leaf_area_index';                                                  UNITS = 'm2 m-2' 
    case ('depth_to_water_table_from_surface_in_peat');  LONG_NAME = 'depth_to_water_table_from_surface_in_peat';              UNITS = 'm'
    case ('free_surface_water_on_peat_flux');  LONG_NAME = 'change_in_free_surface_water_reservoir_on_peat';                   UNITS = 'kg m-2 s-1'
    case ('mwrtm_vegopacity');                 LONG_NAME = 'Lband_microwave_vegopacity_normalized_with_cos_inc_angle';         UNITS = '1'

    ! additional defintions for SMAP Nature Run - reichle, 20 May 2020

    case ('snow_temp_layer1');                 LONG_NAME = 'temperature_top_snow_layer';                                       UNITS = 'K' 
    case ('tb_h');                             LONG_NAME = 'brightness_temperature_land_1410MHz_40deg_Hpol';                   UNITS = 'K' 
    case ('tb_v');                             LONG_NAME = 'brightness_temperature_land_1410MHz_40deg_Vpol';                   UNITS = 'K' 
    case ('TB_LAND_1410MHZ_40DEG_HPOL');       LONG_NAME = 'brightness_temperature_land_1410MHz_40deg_Hpol';                   UNITS = 'K' 
    case ('TB_LAND_1410MHZ_40DEG_VPOL');       LONG_NAME = 'brightness_temperature_land_1410MHz_40deg_Vpol';                   UNITS = 'K' 
 
    ! Done for SM_L4

    case ('Tair');       LONG_NAME = 'air_temperature_at_RefH';                                          UNITS = 'K' 
    case ('TA');         LONG_NAME = 'air_temperature_at_RefH';                                          UNITS = 'K' 
    case ('Qair');       LONG_NAME = 'specific_humidity_at_RefH';                                        UNITS = 'kg kg-1'
    case ('QA');         LONG_NAME = 'specific_humidity_at_RefH';                                        UNITS = 'kg kg-1'
    case ('LWdown');     LONG_NAME = 'surface_absorbed_longwave_flux';                                   UNITS = 'W m-2'
    case ('LWDNSRF');    LONG_NAME = 'surface_absorbed_longwave_flux';                                   UNITS = 'W m-2'
    case ('SWdown');     LONG_NAME = 'downward_shortwave_radiation';                                     UNITS = 'W m-2'
    case ('Wind');       LONG_NAME = 'wind_speed_at_RefH';                                               UNITS = 'm s-1'
    case ('UU');         LONG_NAME = 'wind_speed_at_RefH';                                               UNITS = 'm s-1'
    case ('Psurf');      LONG_NAME = 'surface_pressure';                                                 UNITS = 'Pa'
    case ('PS');         LONG_NAME = 'surface_pressure';                                                 UNITS = 'Pa'
    case ('Rainf_C');    LONG_NAME = 'convective_rainfall';                                              UNITS = 'kg m-2 s-1'
    case ('Rainf');      LONG_NAME = 'liquid_water_precipitation';                                       UNITS = 'kg m-2 s-1'
    case ('Snowf');      LONG_NAME = 'total_snowfall';                                                   UNITS = 'kg m-2 s-1'
    case ('RainfSnowf'); LONG_NAME = 'RainfSnowf';                                                       UNITS = 'kg m-2 s-1'
    case ('SWnet');      LONG_NAME = 'downward_net_shortwave_radiation';                                 UNITS = 'W m-2'
    case ('RefH');       LONG_NAME = 'reference_height_for_Tair_Qair_Wind';                              UNITS = 'm'
    case ('DZ');         LONG_NAME = 'reference_height_for_Tair_Qair_Wind';                              UNITS = 'm'
    case ('CATDEF');     LONG_NAME = 'catchment_deficit';                                                UNITS = 'kg m-2'
    case ('RZEXC');      LONG_NAME = 'root_zone_excess';                                                 UNITS = 'kg m-2'
    case ('SRFEXC');     LONG_NAME = 'surface_excess';                                                   UNITS = 'kg m-2'
    case ('WESNN1');     LONG_NAME = 'snow_mass_layer_1';                                                UNITS = 'kg m-2'
    case ('WESNN2');     LONG_NAME = 'snow_mass_layer_2';                                                UNITS = 'kg m-2'
    case ('WESNN3');     LONG_NAME = 'snow_mass_layer_3';                                                UNITS = 'kg m-2'
    case ('HTSNNN1');    LONG_NAME = 'heat_content_snow_layer_1';                                        UNITS = 'J m-2'
    case ('HTSNNN2');    LONG_NAME = 'heat_content_snow_layer_2';                                        UNITS = 'J m-2'
    case ('HTSNNN3');    LONG_NAME = 'heat_content_snow_layer_3';                                        UNITS = 'J m-2'
    case ('SNDZN1');     LONG_NAME = 'snow_depth_layer_1';                                               UNITS = 'm'
    case ('SNDZN2');     LONG_NAME = 'snow_depth_layer_2';                                               UNITS = 'm'
    case ('SNDZN3');     LONG_NAME = 'snow_depth_layer_3';                                               UNITS = 'm'
    case ('FICE1');      LONG_NAME = 'snow_frozen_fraction_layer_1';                                     UNITS = '1'
    case ('FICE2');      LONG_NAME = 'snow_frozen_fraction_layer_2';                                     UNITS = '1'
    case ('FICE3');      LONG_NAME = 'snow_frozen_fraction_layer_3';                                     UNITS = '1'
    case ('ALBVR');      LONG_NAME = 'surface_reflectivity_for_visible_beam';                            UNITS = '1'
    case ('ALBVF');      LONG_NAME = 'surface_reflectivity_for_visible_diffuse';                         UNITS = '1'
    case ('ALBNR');      LONG_NAME = 'surface_reflectivity_for_near_infared_beam';                       UNITS = '1'
    case ('ALBNF');      LONG_NAME = 'surface_reflectivity_for_near_infrared_diffuse';                   UNITS = '1'
    case ('HLWUP');      LONG_NAME = 'surface_emitted_longwave_flux';                                    UNITS = 'W m-2'
    case ('GWETPROF');   LONG_NAME = 'ave_prof_soil_wetness';                                            UNITS = '1'
    case ('GWETROOT');   LONG_NAME = 'root_zone_soil_wetness';                                           UNITS = '1'
    case ('GWETTOP');    LONG_NAME = 'surface_soil_wetness';                                             UNITS = '1'
    case ('PRMC');       LONG_NAME = 'water_ave_prof';                                                   UNITS = 'm3 m-3'
    case ('RZMC');       LONG_NAME = 'water_root_zone';                                                  UNITS = 'm3 m-3'
    case ('SFMC');       LONG_NAME = 'water_surface_layer';                                              UNITS = 'm3 m-3'
    case ('TPSNOW');     LONG_NAME = 'temperature_top_snow_layer';                                       UNITS = 'K' 
    case ('TUNST');      LONG_NAME = 'temperature_unsaturated_zone';                                     UNITS = 'K' 
    case ('TSAT');       LONG_NAME = 'temperature_saturated_zone';                                       UNITS = 'K' 
    case ('TWLT');       LONG_NAME = 'temperature_wilted_zone';                                          UNITS = 'K' 
    case ('TSURF');      LONG_NAME = 'ave_catchment_temp_incl_snw';                                      UNITS = 'K' 
    case ('TPSURF');     LONG_NAME = 'ave_catchment_temp_incl_snw';                                      UNITS = 'K' 
    case ('GRN');        LONG_NAME = 'greeness_fraction';                                                UNITS = '1' 
    case ('LAI');        LONG_NAME = 'leaf_area_index';                                                  UNITS = '1' 
    case ('TP1');        LONG_NAME = 'soil_temperatures_layer_1';                                        UNITS = 'K'             ! units now K, rreichle & borescan, 6 Nov 2020
    case ('TP2');        LONG_NAME = 'soil_temperatures_layer_2';                                        UNITS = 'K'             ! units now K, rreichle & borescan, 6 Nov 2020
    case ('TP3');        LONG_NAME = 'soil_temperatures_layer_3';                                        UNITS = 'K'             ! units now K, rreichle & borescan, 6 Nov 2020
    case ('TP4');        LONG_NAME = 'soil_temperatures_layer_4';                                        UNITS = 'K'             ! units now K, rreichle & borescan, 6 Nov 2020
    case ('TP5');        LONG_NAME = 'soil_temperatures_layer_5';                                        UNITS = 'K'             ! units now K, rreichle & borescan, 6 Nov 2020
    case ('TP6');        LONG_NAME = 'soil_temperatures_layer_6';                                        UNITS = 'K'             ! units now K, rreichle & borescan, 6 Nov 2020
    case ('PRECTOTLAND');LONG_NAME = 'Total_precipitation_land';                                         UNITS = 'kg m-2 s-1'
    case ('PRECSNOLAND');LONG_NAME = 'snowfall_land';                                                    UNITS = 'kg m-2 s-1'
    case ('SNOWMASS')   ;LONG_NAME = 'snow_mass';                                                        UNITS = 'kg m-2'
    case ('SNOMAS')     ;LONG_NAME = 'snow_mass';                                                        UNITS = 'kg m-2'
    case ('SNO');        LONG_NAME = 'snowfall';                                                         UNITS = 'kg m-2 s-1'
    case ('SNODP');      LONG_NAME = 'snow_depth_in_snow_covered_area';                                  UNITS = 'm'
    case ('EVPSOIL');    LONG_NAME = 'baresoil_evap_energy_flux';                                        UNITS = 'W m-2'
    case ('EVPTRNS');    LONG_NAME = 'transpiration_energy_flux';                                        UNITS = 'W m-2'
    case ('EVPINTR');    LONG_NAME = 'interception_loss_energy_flux';                                    UNITS = 'W m-2'
    case ('EVPSBLN');    LONG_NAME = 'snow_ice_evaporation_energy_flux';                                 UNITS = 'W m-2'
    case ('RUNOFF');     LONG_NAME = 'runoff_flux';                                                      UNITS = 'kg m-2 s-1'
    case ('BASEFLOW');   LONG_NAME = 'baseflow_flux';                                                    UNITS = 'kg m-2 s-1'
    case ('SMLAND');     LONG_NAME = 'Snowmelt_flux_land';                                               UNITS = 'kg m-2 s-1'
    case ('QINFIL');     LONG_NAME = 'rainwater_infiltration_flux';                                      UNITS = 'kg m-2 s-1'
    case ('FRUNST');     LONG_NAME = 'fractional_area_of_unsaturated_zone';                              UNITS = '1'
    case ('FRSAT');      LONG_NAME = 'fractional_area_of_saturated_zone';                                UNITS = '1'
    case ('FRSNO');      LONG_NAME = 'fractional_area_of_land_snowcover';                                UNITS = '1'
    case ('FRWLT');      LONG_NAME = 'fractional_area_of_wilting_zone';                                  UNITS = '1'
    case ('PARDFLAND');  LONG_NAME = 'surface_downwelling_par_diffuse_flux';                             UNITS = 'W m-2'
    case ('PARDRLAND');  LONG_NAME = 'surface_downwelling_par_beam_flux';                                UNITS = 'W m-2'
    case ('SHLAND');     LONG_NAME = 'Sensible_heat_flux_land';                                          UNITS = 'W m-2'
    case ('LHLAND');     LONG_NAME = 'Latent_heat_flux_land';                                            UNITS = 'W m-2'
    case ('EVLAND');     LONG_NAME = 'Evaporation_land';                                                 UNITS = 'kg m-2 s-1'
    case ('LWLAND');     LONG_NAME = 'Net_longwave_land';                                                UNITS = 'W m-2'
    case ('SWLAND');     LONG_NAME = 'Net_shortwave_land';                                               UNITS = 'W m-2'
    case ('SWDOWNLAND'); LONG_NAME = 'Incident_shortwave_land';                                          UNITS = 'W m-2'
    case ('GHLAND');     LONG_NAME = 'Ground_heating_land';                                              UNITS = 'W m-2'
    case ('TWLAND');     LONG_NAME = 'Avail_water_storage_land';                                         UNITS = 'kg m-2'
    case ('TSLAND');     LONG_NAME = 'Total_snow_storage_land';                                          UNITS = 'kg m-2'
    case ('TELAND');     LONG_NAME = 'Total_energy_storage_land';                                        UNITS = 'J m-2'
    case ('WCHANGE');    LONG_NAME = 'rate_of_change_of_total_land_water';                               UNITS = 'kg m-2 s-1'
    case ('ECHANGE');    LONG_NAME = 'rate_of_change_of_total_land_energy';                              UNITS = 'W m-2'
    case ('SPLAND');     LONG_NAME = 'rate_of_spurious_land_energy_source';                              UNITS = 'W m-2'
    case ('SPWATR');     LONG_NAME = 'rate_of_spurious_land_water_source';                               UNITS = 'kg m-2 s-1'
    case ('SPSNOW');     LONG_NAME = 'rate_of_spurious_snow_energy';                                     UNITS = 'W m-2'
    case ('PEATCLSM_WATERLEVEL');LONG_NAME = 'depth_to_water_table_from_surface_in_peat';                UNITS = 'm'
    case ('PEATCLSM_FSWCHANGE'); LONG_NAME = 'change_in_free_surface_water_reservoir_on_peat';           UNITS = 'kg m-2 s-1'
    case ('CNLAI');      LONG_NAME = 'CN_exposed_leaf-area_index';                                       UNITS = '1'
    case ('CNTLAI');     LONG_NAME = 'CN_total_leaf-area_index';                                         UNITS = '1'
    case ('CNSAI');      LONG_NAME = 'CN_exposed_stem-area_index';                                       UNITS = '1'
    case ('CNTOTC');     LONG_NAME = 'CN_total_carbon';                                                  UNITS = 'kg m-2'
    case ('CNVEGC');     LONG_NAME = 'CN_total_vegetation_carbon';                                       UNITS = 'kg m-2'
    case ('CNROOT');     LONG_NAME = 'CN_total_root_carbon';                                             UNITS = 'kg m-2'
    case ('CNNPP');      LONG_NAME = 'CN_net_primary_production';                                        UNITS = 'kg m-2 s-1'
    case ('CNGPP');      LONG_NAME = 'CN_gross_primary_production';                                      UNITS = 'kg m-2 s-1'
    case ('CNSR');       LONG_NAME = 'CN_total_soil_respiration';                                        UNITS = 'kg m-2 s-1'
    case ('CNNEE');      LONG_NAME = 'CN_net_ecosystem_exchange';                                        UNITS = 'kg m-2 s-1'
    case ('CNXSMR');     LONG_NAME = 'abstract_C_pool_to_meet_excess_MR_demand';                         UNITS = 'kg m-2'
    case ('CNADD');      LONG_NAME = 'CN_added_to_maintain_positive_C';                                  UNITS = 'kg m-2 s-1'
    case ('PARABS');     LONG_NAME = 'absorbed_PAR';                                                     UNITS = 'W m-2'
    case ('PARINC');     LONG_NAME = 'incident_PAR';                                                     UNITS = 'W m-2'
    case ('SCSAT');      LONG_NAME = 'saturated_stomatal_conductance';                                   UNITS = 'm s-1'
    case ('SCUNS');      LONG_NAME = 'unstressed_stomatal_conductance';                                  UNITS = 'm s-1'
    case ('BTRAN');      LONG_NAME = 'transpiration coefficient';                                        UNITS = '1'
    case ('SIF');        LONG_NAME = 'solar induced fluorescence';                                       UNITS = 'umol m-2 sm s-1'
    case ('CLOSS');      LONG_NAME = 'CN_carbon_loss_to_fire';                                           UNITS = 'kg m-2 s-1'
    case ('BURN');       LONG_NAME = 'CN_fractional_area_burn_rate';                                     UNITS = 's-1'
    case ('FSEL');       LONG_NAME = 'fire season length';                                               UNITS = 'days'
    case ('EVPSNO');     LONG_NAME = 'snowpack_evaporation_energy_flux';                                 UNITS = 'W m-2'
    case ('GHTSKIN');    LONG_NAME = 'Ground_heating_skin_temp';                                         UNITS = 'W m-2'
    case ('WAT10CM');    LONG_NAME = 'soil moisture in Upper 10cm';                                      UNITS = 'kg m-2'
    case ('WATSOI');     LONG_NAME = 'totoal soil moisture';                                             UNITS = 'kg m-2'
    case ('ICESOI');     LONG_NAME = 'soil frozen water content';                                        UNITS = 'kg m-2'
    case ('RMELTDU001'); LONG_NAME = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_1';           UNITS = 'kg m-2 s-1'
    case ('RMELTDU002'); LONG_NAME = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_2';           UNITS = 'kg m-2 s-1'
    case ('RMELTDU003'); LONG_NAME = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_3';           UNITS = 'kg m-2 s-1'
    case ('RMELTDU004'); LONG_NAME = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_4';           UNITS = 'kg m-2 s-1'
    case ('RMELTDU005'); LONG_NAME = 'flushed_out_dust_mass_flux_from_the_bottom_layer_bin_5';           UNITS = 'kg m-2 s-1'
    case ('RMELTBC001'); LONG_NAME = 'flushed_out_black_carbon_mass_flux_from_the_bottom_layer_bin_1';   UNITS = 'kg m-2 s-1'
    case ('RMELTBC002'); LONG_NAME = 'flushed_out_black_carbon_mass_flux_from_the_bottom_layer_bin_2';   UNITS = 'kg m-2 s-1'
    case ('RMELTOC001'); LONG_NAME = 'flushed_out_organic_carbon_mass_flux_from_the_bottom_layer_bin_1'; UNITS = 'kg m-2 s-1'
    case ('RMELTOC002'); LONG_NAME = 'flushed_out_organic_carbon_mass_flux_from_the_bottom_layer_bin_2'; UNITS = 'kg m-2 s-1'
 
    ! land assimilation increments for Catchment prognostic variables in coupled land-atmosphere DAS (#sqz 2020-01)

    case ('TCFSAT_INCR');  LONG_NAME = 'increment_canopy_temperature_saturated_zone';                    UNITS = 'K' 
    case ('TCFTRN_INCR');  LONG_NAME = 'increment_canopy_temperature_transition_zone';                   UNITS = 'K'
    case ('TCFWLT_INCR');  LONG_NAME = 'increment_canopy_temperature_wilting_zone';                      UNITS = 'K'
    case ('QCFSAT_INCR');  LONG_NAME = 'increment_canopy_specific_humidity_saturated_zone';              UNITS = 'kg kg-1'
    case ('QCFTRN_INCR');  LONG_NAME = 'increment_canopy_specific_humidity_transition_zone';             UNITS = 'kg kg-1'
    case ('QCFWLT_INCR');  LONG_NAME = 'increment_canopy_specific_humidity_wilting_zone';                UNITS = 'kg kg-1'
    case ('CAPAC_INCR');   LONG_NAME = 'increment_interception_reservoir_capac';                         UNITS = 'kg m-2'
    case ('CATDEF_INCR');  LONG_NAME = 'increment_catchment_deficit';                                    UNITS = 'kg m-2'
    case ('RZEXC_INCR');   LONG_NAME = 'increment_root_zone_excess';                                     UNITS = 'kg m-2'
    case ('SRFEXC_INCR');  LONG_NAME = 'increment_surface_excess';                                       UNITS = 'kg m-2'
    case ('GHTCNT1_INCR'); LONG_NAME = 'increment_soil_heat_content_layer_1';                            UNITS = 'J m-2'
    case ('GHTCNT2_INCR'); LONG_NAME = 'increment_soil_heat_content_layer_2';                            UNITS = 'J m-2'
    case ('GHTCNT3_INCR'); LONG_NAME = 'increment_soil_heat_content_layer_3';                            UNITS = 'J m-2' 
    case ('GHTCNT4_INCR'); LONG_NAME = 'increment_soil_heat_content_layer_4';                            UNITS = 'J m-2'
    case ('GHTCNT5_INCR'); LONG_NAME = 'increment_soil_heat_content_layer_5';                            UNITS = 'J m-2'
    case ('GHTCNT6_INCR'); LONG_NAME = 'increment_soil_heat_content_layer_6';                            UNITS = 'J m-2'
    case ('WESNN1_INCR');  LONG_NAME = 'increment_snow_mass_layer_1';                                    UNITS = 'kg m-2'
    case ('WESNN2_INCR');  LONG_NAME = 'increment_snow_mass_layer_2';                                    UNITS = 'kg m-2'
    case ('WESNN3_INCR');  LONG_NAME = 'increment_snow_mass_layer_3';                                    UNITS = 'kg m-2'
    case ('HTSNNN1_INCR'); LONG_NAME = 'increment_heat_content_snow_layer_1';                            UNITS = 'J m-2'
    case ('HTSNNN2_INCR'); LONG_NAME = 'increment_heat_content_snow_layer_2';                            UNITS = 'J m-2'
    case ('HTSNNN3_INCR'); LONG_NAME = 'increment_heat_content_snow_layer_3';                            UNITS = 'J m-2'
    case ('SNDZN1_INCR');  LONG_NAME = 'increment_snow_depth_layer_1';                                   UNITS = 'm'
    case ('SNDZN2_INCR');  LONG_NAME = 'increment_snow_depth_layer_2';                                   UNITS = 'm'
    case ('SNDZN3_INCR');  LONG_NAME = 'increment_snow_depth_layer_3';                                   UNITS = 'm'

    ! land assimilation forecast and analysis for Catchment model diagnostics

    case ('SFMC_FCST');          LONG_NAME = 'soil_moisture_surface_forecast';                           UNITS = 'm3 m-3'
    case ('RZMC_FCST');          LONG_NAME = 'soil_moisture_rootzone_forecast';                          UNITS = 'm3 m-3'
    case ('PRMC_FCST');          LONG_NAME = 'soil_moisture_profile_forecast';                           UNITS = 'm3 m-3'
    case ('TSURF_FCST');         LONG_NAME = 'ave_catchment_temp_incl_snw_forecast';                     UNITS = 'K'
    case ('TSOIL1_FCST');        LONG_NAME = 'soil_temperatures_layer_1_forecast';                       UNITS = 'K'

    case ('SFMC_FCST_ENSSTD');   LONG_NAME = 'soil_moisture_surface_forecast_ensstd';                    UNITS = 'm3 m-3'
    case ('RZMC_FCST_ENSSTD');   LONG_NAME = 'soil_moisture_rootzone_forecast_ensstd';                   UNITS = 'm3 m-3'
    case ('PRMC_FCST_ENSSTD');   LONG_NAME = 'soil_moisture_profile_forecast_ensstd';                    UNITS = 'm3 m-3'
    case ('TSURF_FCST_ENSSTD');  LONG_NAME = 'ave_catchment_temp_incl_snw_forecast_ensstd';              UNITS = 'K'
    case ('TSOIL1_FCST_ENSSTD'); LONG_NAME = 'soil_temperatures_layer_1_forecast_ensstd';                UNITS = 'K'

    case ('SFMC_ANA');           LONG_NAME = 'soil_moisture_surface_analysis';                           UNITS = 'm3 m-3'
    case ('RZMC_ANA');           LONG_NAME = 'soil_moisture_rootzone_analysis';                          UNITS = 'm3 m-3'
    case ('PRMC_ANA');           LONG_NAME = 'soil_moisture_profile_analysis';                           UNITS = 'm3 m-3'
    case ('TSURF_ANA');          LONG_NAME = 'ave_catchment_temp_incl_snw_analysis';                     UNITS = 'K'
    case ('TSOIL1_ANA');         LONG_NAME = 'soil_temperatures_layer_1_analysis';                       UNITS = 'K'

    case ('SFMC_ANA_ENSSTD');    LONG_NAME = 'soil_moisture_surface_analysis_ensstd';                    UNITS = 'm3 m-3'
    case ('RZMC_ANA_ENSSTD');    LONG_NAME = 'soil_moisture_rootzone_analysis_ensstd';                   UNITS = 'm3 m-3'
    case ('PRMC_ANA_ENSSTD');    LONG_NAME = 'soil_moisture_profile_analysis_ensstd';                    UNITS = 'm3 m-3'
    case ('TSURF_ANA_ENSSTD');   LONG_NAME = 'ave_catchment_temp_incl_snw_analysis_ensstd';              UNITS = 'K'
    case ('TSOIL1_ANA_ENSSTD');  LONG_NAME = 'soil_temperatures_layer_1_analysis_ensstd';                UNITS = 'K'

    ! other land assimilation fields
       
    case ('MWRTM_VEGOPACITY'); LONG_NAME = 'Lband_microwave_vegopacity_normalized_with_cos_inc_angle';   UNITS = '1'       
       
    ! default LONG_NAME and UNITS for nc4 files created by tile_bin2nc4.F90 (used for any SHORT_NAME not listed above):

    case default;              LONG_NAME = 'not defined in tile_bin2nc4.F90';                            UNITS = 'not defined in tile_bin2nc4.F90';

    end select

    if (present(LNAME)) str_atr = trim (LONG_NAME)
    if (present(UNT))   str_atr = trim (UNITS    )

  END FUNCTION getAttribute

END PROGRAM tile_bin2nc4
