function [ veg_lookup, soil_lookup ] = get_mwRTM_lookup(option);

% Gabrielle De Lannoy, GSFC, 22Jul11
%
% Same lookup table as in LDASsa
%
% Updated: 
% GDL, 6  Sept 2011:  16 IGBP vegetation classes
% GDL, 27 Sept 2011:  h=1.2 for conif forest (paper J. Grant), 
%                     rather than h=1.6 (CMEM)
%==================================================================
% VEGETATION 
%==================================================================
%
% N_vegcls = 8;
% 
%    ! 1 Broadleaf evergreen trees   0.1    0.1  0.12  0.14  0.10  0.30  1.0    0.0
%    ! 2 Broadleaf deciduous trees   0.1    0.1  0.12  0.14  0.10  0.2   1.0    2.0
%    ! 3 Needleleaf trees            0.1    0.1  0.12  0.12  0.08  0.1   1.75   0.0
%    ! 4 Grassland                   0.1    0.1  0.05  0.11  0.09  0.15  1.0    0.0
%    ! 5 Broadleaf shrubs            0.1    0.1  0.12  0.12  0.10  0.2   1.0    1.0
%    ! 6 Dwarf trees                 0.1    0.1  0.05  0.11  0.09  0.15  1.0    1.0
%    ! 7 Bare soil                   0.1    0.1  0.00  0.00  0.00  0.    0.0   -1.0
%    ! 8 Desert soil                 0.1    0.1  0.00  0.00  0.00  0.    0.0   -1.0
%
% veg_lookup.rgh_hmin = [ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 ]*0.5;
% veg_lookup.rgh_hmax = [ 1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0,  1.0 ]*0.5;
% veg_lookup.omega   = [ 0.12,  0.12,  0.12,  0.05,  0.12,  0.05,  0.00,  0.00 ];
% veg_lookup.lewt = [ 0.30,  0.20,  0.10,  0.15,  0.20,  0.15,  0.00,  0.00 ];
% veg_lookup.bh   = [ 0.10,  0.10,  0.08,  0.09,  0.10,  0.09,  0.00,  0.00 ];
% veg_lookup.bv   = [ 0.14,  0.14,  0.12,  0.11,  0.12,  0.11,  0.00,  0.00 ];
% veg_lookup.rgh_Nrh  = [ 1.00,  1.00,  1.75,  1.00,  1.00,  1.00,  0.00,  0.00 ];
% veg_lookup.rgh_Nrv  = [ 0.00,  2.00,  0.00,  0.00,  1.00,  1.00, -1.00, -1.00 ];
% veg_lookup.tag  = {'rgh_hmin [-]', 'h_ [-]', '\omega [-]', 'lewt [-]',...
%                    'b_h [-]', 'b_v [-]', 'Nr_h [-]', 'Nr_h [-]'};

%1	Evergreen Needleleaf Forest
%2	Evergreen Broadleaf Forest
%3	Deciduous Needleleaf Forest
%4	Deciduous Broadleaf Forest
%5	Mixed Forest
%6	Closed Shrublands
%7	Open Shrublands
%8	Woody Savannas
%9	Savannas
%10	Grasslands
%11	Permanent Wetlands
%12	Croplands
%13	Urban and Built-Up
%14	Cropland & Natural Vegetation
%15	Snow and Ice
%16	Barren or Sparsely Vegetated


N_vegcls = 16;
veg_lookup.tag  = {'rgh_hmin [-]', 'rgh_hmax [-]',  '\omega [-]', 'lewt [-]',...
                   'b_h [-]',  'b_v [-]', 'Nr_h [-]',   'Nr_h [-]'};

if (strcmp(option,'CMEM') || strcmp(option,'Lit2'))
    
    %ECMWF CMEM-code
    %veg_lookup.rgh_hmin = [ 1.6, 1.3, 1.6, 1, 1.3, 0.7, 0.7, 0.7, 0.5, 0.1, 0.1, 0.5, 0, 0.7, 0, 0.1 ];
    %veg_lookup.rgh_hmax = [ 1.6, 1.3, 1.6, 1, 1.3, 0.7, 0.7, 0.7, 0.5, 0.1, 0.1, 0.5, 0, 0.7, 0, 0.1 ];
    veg_lookup.rgh_hmin  = [ 1.2, 1.3, 1.2, 1, 1.3, 0.7, 0.7, 0.7, 0.5, 0.1, 0.1, 0.5, 0, 0.7, 0, 0.1 ];
    veg_lookup.rgh_hmax  = [ 1.2, 1.3, 1.2, 1, 1.3, 0.7, 0.7, 0.7, 0.5, 0.1, 0.1, 0.5, 0, 0.7, 0, 0.1 ];
    
    veg_lookup.omega     = [ 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05 ];
    veg_lookup.lewt      = [ 1,   1,   1,   1,   1,    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,   0, 0.5,   0, 0 ];
    veg_lookup.bh        = [ 0.33,0.33,0.33,0.33,0.33, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.15,  0, 0.15,  0, 0 ];
    veg_lookup.bv        = [ 0.33,0.33,0.33,0.33,0.33, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.15,  0, 0.15,  0, 0 ];
    veg_lookup.rgh_Nrh   = [ 1,   1.75,1,   1,   1,    1,   1,   1,   1,   1,   1,    0,   1,    0,   1, 0 ];
    veg_lookup.rgh_Nrv   = [ 0,     0, 0,   2,   1,    0,   0,   0,   0,   0,   0,   -1,   1,   -1,   1,-1 ];
    veg_lookup.st_scale  = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ];
    
elseif (strcmp(option,'CMEM_SMOS') || strcmp(option,'Lit3'))
    
    %ECMWF SMOS monitoring-setup
    veg_lookup.rgh_hmin = [ 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66 ];
    veg_lookup.rgh_hmax = [ 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66, 1.66 ];
    
    veg_lookup.omega    = [ 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05 ];
    veg_lookup.lewt     = [ 1,   1,   1,   1,   1,    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0, 0.5, 0, 0 ];
    veg_lookup.bh       = [ 0.33,0.33,0.33,0.33,0.33, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.15,0, 0.15,0, 0 ];
    veg_lookup.bv       = [ 0.33,0.33,0.33,0.33,0.33, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.15,0, 0.15,0, 0 ];
    veg_lookup.rgh_Nrh  = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ];
    veg_lookup.rgh_Nrv  = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ];
    veg_lookup.st_scale = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ];

elseif (strcmp(option,'SMAP') || strcmp(option,'Lit1'))
    
    %Peggy O'Neill's SMAP ATBD, Table 2               
    veg_lookup.rgh_hmin = [ 0.16, 0.16, 0.16, 0.16, 0.16, 0.11, 0.11, 0.125, 0.156, 0.156, 0.156, 0.108, 0, 0.13, 0, 0.15 ];
    veg_lookup.rgh_hmax = [ 0.16, 0.16, 0.16, 0.16, 0.16, 0.11, 0.11, 0.125, 0.156, 0.156, 0.156, 0.108, 0, 0.13, 0, 0.15 ];
    veg_lookup.omega    = [ 0.12, 0.12, 0.12, 0.12, 0.08, 0.05, 0.05, 0.12,  0.08,  0.05,  0.05,  0.05,  0, 0.065,0, 0 ];
    veg_lookup.lewt     = [ 0.3,  0.3,  0.2,  0.2,  0.2,  0.2,  0.2,  0.15,  0.15,  0.15,  0.15,  0.15,  0, 0.15, 0, 0 ];
    veg_lookup.bh       = [ 0.1,  0.1,  0.12, 0.12, 0.12, 0.11, 0.11, 0.11,  0.11,   0.1,   0.1,  0.11,  0, 0.11, 0, 0 ];
    veg_lookup.bv       = [ 0.1,  0.1,  0.12, 0.12, 0.12, 0.11, 0.11, 0.11,  0.11,   0.1,   0.1,  0.11,  0, 0.11, 0, 0 ];
    veg_lookup.rgh_Nrh  = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ];
    veg_lookup.rgh_Nrv  = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ];
    veg_lookup.st_scale = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ];

elseif (strcmp(option,'Lit4'))

    %ECMWF CMEM-code; same as Lit2, but with new lewt values to go with new LAI values
    %veg_lookup.rgh_hmin = [ 1.6, 1.3, 1.6, 1, 1.3, 0.7, 0.7, 0.7, 0.5, 0.1, 0.1, 0.5, 0, 0.7, 0, 0.1 ];
    %veg_lookup.rgh_hmax = [ 1.6, 1.3, 1.6, 1, 1.3, 0.7, 0.7, 0.7, 0.5, 0.1, 0.1, 0.5, 0, 0.7, 0, 0.1 ];
    veg_lookup.rgh_hmin = [ 1.2, 1.3, 1.2, 1, 1.3, 0.7, 0.7, 0.7, 0.5, 0.1, 0.1, 0.5, 0, 0.7, 0, 0.1 ];
    veg_lookup.rgh_hmax = [ 1.2, 1.3, 1.2, 1, 1.3, 0.7, 0.7, 0.7, 0.5, 0.1, 0.1, 0.5, 0, 0.7, 0, 0.1 ];

    veg_lookup.omega    = [ 0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05 ];
    veg_lookup.lewt     = [ 1.9, 1.5, 1.5, 1.5, 1.7,  0.9, 0.8, 0.9, 0.8, 0.9, 0.8, 0.8,   0, 0.9,   0, 0 ];
   %veg_lookup.lewt = [ 1,   1,   1,   1,   1,    0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,   0, 0.5,   0, 0 ];
    veg_lookup.bh       = [ 0.33,0.33,0.33,0.33,0.33, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.15,  0, 0.15,  0, 0 ];
    veg_lookup.bv       = [ 0.33,0.33,0.33,0.33,0.33, 0.3, 0.3, 0.3, 0.2, 0.2, 0.2, 0.15,  0, 0.15,  0, 0 ];
    veg_lookup.rgh_Nrh  = [ 1,   1.75,1,   1,   1,    1,   1,   1,   1,   1,   1,    0,   1,    0,   1, 0 ];
    veg_lookup.rgh_Nrv  = [ 0,     0, 0,   2,   1,    0,   0,   0,   0,   0,   0,   -1,   1,   -1,   1,-1 ];
    veg_lookup.st_scale = [ 0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0 ];


elseif (strcmp(option,'Cal3c'))

   load('/hydro/gdelanno/proc_analysis/initial_lookup_3c_rgh_hmin_rgh_hmax_rgh_Nrh_rgh_Nrv_lewt_omega_bh_bv_st_scale.mat');

%was Cal2c before
elseif (strcmp(option,'CalD2'))

   load('/hydro/gdelanno/proc_analysis/initial_lookup_D2_rgh_hmin_rgh_hmax_rgh_Nrh_rgh_Nrv_lewt_omega_bh_bv_st_scale.mat');
   
else    
    
    error(['ERROR: Option ',option,' does not exist'])        
    
end

%==================================================================
% SOIL 
%==================================================================              
%1	Sand
%2	Loamy Sand
%3	Sandy Loam
%4	Loam (F)       ==> Silt Loam
%5	Silt Loam (F)  ==> Silt
%6	Silt (F)       ==> Loam
%7	Sandy Clay Loam
%8	Clay Loam (F)       ==> Silty Clay Loam 
%9	Silty Clay Loam (F) ==> Clay Loam
%10	Sandy Clay Loam
%11	Silty Clay Loam
%12	Clay


N_soilcls = 12;

soil_lookup.sf    = [.92, .82, .58, .17, .10, .43, .58, .10, .32, .52, .06, .22];
soil_lookup.cf    = [.03, .06, .10, .13, .05, .18, .27, .34, .34, .42, .47, .58];
soil_lookup.fc    = [0.132, 0.156, 0.196, 0.27 , 0.361, 0.25 , 0.253, 0.334, 0.301, 0.288, 0.363, 0.353];
soil_lookup.wp    = [0.033, 0.051, 0.086, 0.169, 0.045, 0.148, 0.156, 0.249, 0.211, 0.199, 0.286, 0.276];
soil_lookup.poros = [0.373, 0.386, 0.419, 0.476, 0.471, 0.437, 0.412, 0.478, 0.447, 0.415, 0.478, 0.45];
soil_lookup.b     = [3.3,   3.8,   4.34,  5.25,  3.63,  5.96,  7.32,  8.41,  8.34,  9.7,  10.78, 12.93]; 
soil_lookup.PsiS  = [-0.05,	-0.07, -0.16, -0.65, -0.84, -0.24, -0.12, -0.63, -0.28, -0.12, -0.58, -0.27];
soil_lookup.Ksat  = [2.45E-05, 1.75E-05, 8.35E-06, 2.36E-06, 1.10E-06, 4.66E-06, 6.31E-06, 1.44E-06, 2.72E-06, 4.25E-06, 1.02E-06, 1.33E-06];
soil_lookup.tag   = {'Sand [-]', 'Clay [-]', 'Field Capacity [m3/m3]', ...
                     'Wilting Point [m3/m3]', 'Porosity [m3/m3]', 'b [ ]', 'Psisat [m]', 'Ksat [m/s]'};

end

%==========================EOF=====================================
