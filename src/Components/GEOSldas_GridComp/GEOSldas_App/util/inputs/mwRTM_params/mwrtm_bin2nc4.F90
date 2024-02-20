!This program converts original mwrtm_param.bin to nc4.
!The original mwrtm_param.bin is the same order as the tile file and encode with big endian 
!for example /gpfsm/dnb31/gdelanno/input/RTM_parms/EASEv2/L4SM_v001_Lit4_CalD0/SMAP_EASEv2_M36/mwRTM_param.bin
!
! Note: "mwrtm_param" file only includes time-invariant mwRTM parameters (i.e., excl. vegopacity)
!       -reichle, 21 Feb 2022
!
PROGRAM mwrtm_bin2nc4
  use mwRTM_types, only: mwRTM_param_type
  use mwRTM_types, only: mwRTM_param_nodata_check 

  implicit none
  INCLUDE 'netcdf.inc'

  integer       :: i,k,  n, command_argument_count, NTILES
  integer       :: NCFOutID, Vid, STATUS, CellID, TimID, nVars
  character(512):: Usage="mwrtm_bin2nc4.x mwrtm_BINFILE  mwRTM_param.nc4"
  character(512):: BINFILE, MWRTMNC4, arg(3)
  real, allocatable, dimension (:)    :: var
  integer, allocatable,dimension (:)  :: NT
  character(len=:),allocatable :: shnms(:)
  type(mwRTM_param_type), allocatable :: mwp(:)
  integer :: unitnum
  logical :: mwp_nodata


  nVars = 18
  shnms = [              &
       'MWRTM_VEGCLS   ',&
       'MWRTM_SOILCLS  ',&
       'MWRTM_SAND     ',&
       'MWRTM_CLAY     ',&
       'MWRTM_POROS    ',&
       'MWRTM_WANGWT   ',&
       'MWRTM_WANGWP   ',&
       'MWRTM_RGHHMIN  ',&
       'MWRTM_RGHHMAX  ',&
       'MWRTM_RGHWMIN  ',&
       'MWRTM_RGHWMAX  ',&
       'MWRTM_RGHNRH   ',&
       'MWRTM_RGHNRV   ',&
       'MWRTM_RGHPOLMIX',&
       'MWRTM_OMEGA    ',&
       'MWRTM_BH       ',&
       'MWRTM_BV       ',&
       'MWRTM_LEWT     ']
  
  ! processing command line agruments
  I = command_argument_count()
  
  if( I /=2 ) then
     print *, "Wrong Number of arguments: ", i
     print *, trim(Usage)
     stop
  end if

  do n=1,I
     call get_command_argument(n,arg(n))
  enddo

  read(arg(1),'(a)') BINFILE
  read(arg(2),'(a)') MWRTMNC4

  print *,trim(BINFILE)
  print *,trim(MWRTMNC4)

 ! reading mwrtm_bin

  unitnum = 10
  open (unitnum, file = trim(BINFILE), CONVERT='BIG_ENDIAN',form = 'unformatted', action ='read')
  read (unitnum) NTILES
  print*, "Ntiles", NTILES
  allocate(var(NTILES))
  allocate(NT(NTILES))

  status = NF_CREATE (trim(MWRTMNC4), NF_NETCDF4, NCFOutID)
  status = NF_DEF_DIM(NCFOutID, 'tile' , NTILES, CellID)

  do n = 1, nVars

     status = NF_DEF_VAR(NCFOutID,trim(shnms(n)) , NF_FLOAT, 1 ,CellID, vid)
     status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'long_name',               &
          LEN_TRIM(getAttribute(shnms(n), LNAME = 1)),  &
          getAttribute(shnms(n), LNAME = 1) )
     status = NF_PUT_ATT_TEXT(NCFOutID, vid, 'units',                   &
          LEN_TRIM(getAttribute(shnms(n), UNT = 1)),    &
          getAttribute(shnms(n), UNT = 1))
  end do

  ! reading and writing
  ! tile id 
  read(unitnum) NT ! read off the tile id
  do i= 1, NTILES
     if (i /= NT(i)) stop "not original one"
  enddo

  allocate(mwp(NTILES))

  read (unitnum) NT;   mwp(1:NTILES)%vegcls     =  NT(1:NTILES)
  read (unitnum) NT;   mwp(1:NTILES)%soilcls    =  NT(1:NTILES)

  read (unitnum) VAR;  mwp(1:NTILES)%sand       =  VAR(1:NTILES)
  read (unitnum) VAR;  mwp(1:NTILES)%clay       =  VAR(1:NTILES)

  read (unitnum) VAR;  mwp(1:NTILES)%poros      =  VAR(1:NTILES)

  read (unitnum) VAR;  mwp(1:NTILES)%wang_wt    =  VAR(1:NTILES)
  read (unitnum) VAR;  mwp(1:NTILES)%wang_wp    =  VAR(1:NTILES)

  read (unitnum) VAR;  mwp(1:NTILES)%rgh_hmin   =  VAR(1:NTILES)
  read (unitnum) VAR;  mwp(1:NTILES)%rgh_hmax   =  VAR(1:NTILES)
  read (unitnum) VAR;  mwp(1:NTILES)%rgh_wmin   =  VAR(1:NTILES)
  read (unitnum) VAR;  mwp(1:NTILES)%rgh_wmax   =  VAR(1:NTILES)
  read (unitnum) VAR;  mwp(1:NTILES)%rgh_Nrh    =  VAR(1:NTILES)
  read (unitnum) VAR;  mwp(1:NTILES)%rgh_Nrv    =  VAR(1:NTILES)
  read (unitnum) VAR;  mwp(1:NTILES)%rgh_polmix =  VAR(1:NTILES)

  read (unitnum) VAR;  mwp(1:NTILES)%omega      =  VAR(1:NTILES)

  read (unitnum) VAR;  mwp(1:NTILES)%bh         =  VAR(1:NTILES)
  read (unitnum) VAR;  mwp(1:NTILES)%bv         =  VAR(1:NTILES)
  read (unitnum) VAR;  mwp(1:NTILES)%lewt       =  VAR(1:NTILES)

  do i = 1, NTILES
    call mwRTM_param_nodata_check( mwp(i), mwp_nodata )
  enddo

  VAR = real(mwp(1:NTILES)%vegcls)
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(1))) ,(/1/),(/NTILES/),var )
  VAR = real(mwp(1:NTILES)%soilcls)
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(2))) ,(/1/),(/NTILES/),var )

  VAR =  mwp(1:NTILES)%sand
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(3))) ,(/1/),(/NTILES/),var )
  VAR =  mwp(1:NTILES)%clay
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(4))) ,(/1/),(/NTILES/),var )

  VAR =  mwp(1:NTILES)%poros
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(5))) ,(/1/),(/NTILES/),var )

  VAR =  mwp(1:NTILES)%wang_wt
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(6))) ,(/1/),(/NTILES/),var )
  VAR =  mwp(1:NTILES)%wang_wp 
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(7))) ,(/1/),(/NTILES/),var )

  VAR =  mwp(1:NTILES)%rgh_hmin
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(8))) ,(/1/),(/NTILES/),var )
  VAR =  mwp(1:NTILES)%rgh_hmax
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(9))) ,(/1/),(/NTILES/),var )
  VAR =  mwp(1:NTILES)%rgh_wmin
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(10))) ,(/1/),(/NTILES/),var )
  VAR =  mwp(1:NTILES)%rgh_wmax
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(11))) ,(/1/),(/NTILES/),var )
  VAR =  mwp(1:NTILES)%rgh_Nrh
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(12))) ,(/1/),(/NTILES/),var )
  VAR =  mwp(1:NTILES)%rgh_Nrv 
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(13))) ,(/1/),(/NTILES/),var )
  VAR =  mwp(1:NTILES)%rgh_polmix
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(14))) ,(/1/),(/NTILES/),var )

  VAR =  mwp(1:NTILES)%omega
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(15))) ,(/1/),(/NTILES/),var )

  VAR =  mwp(1:NTILES)%bh  
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(16))) ,(/1/),(/NTILES/),var )
  VAR =  mwp(1:NTILES)%bv
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(17))) ,(/1/),(/NTILES/),var )
  VAR =  mwp(1:NTILES)%lewt  
  status = NF_PUT_VARA_REAL(NCFOutID,VarID(NCFOutID,trim(shnms(18))) ,(/1/),(/NTILES/),var )

  
  STATUS   = NF_CLOSE (NCFOutID)
  close (10)

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
    case('MWRTM_VEGCLS');    LONG_NAME = 'L-band RTM model: Vegetation class. Type is Unsigned32';             UNITS = '1'
    case('MWRTM_SOILCLS');   LONG_NAME = 'L-band RTM model: Soil class. Type is Unsigned32';                   UNITS = '1'
    case('MWRTM_SAND');      LONG_NAME = 'L-band RTM model: Sand fraction';                                    UNITS = '1'
    case('MWRTM_CLAY');      LONG_NAME = 'L-band RTM model: Clay fraction';                                    UNITS = '1'
    case('MWRTM_POROS');     LONG_NAME = 'L-band RTM model: Porosity';                                         UNITS = 'm3 m-3'
    case('MWRTM_WANGWT');    LONG_NAME = 'L-band RTM model: Wang dielectric model transition soil moisture';   UNITS = 'm3 m-3'
    case('MWRTM_WANGWP');    LONG_NAME = 'L-band RTM model: Wang dielectric model wilting point soil moisture';UNITS = 'm3 m-3'
    case('MWRTM_RGHHMIN');   LONG_NAME = 'L-band RTM model: Minimum microwave roughness parameter'; UNITS = '1'
    case('MWRTM_RGHHMAX');   LONG_NAME = 'L-band RTM model: Maximum microwave roughness parameter'; UNITS = '1'
    case('MWRTM_RGHWMIN');   LONG_NAME = 'L-band RTM model: Soil moisture value below which maximum microwave roughness parameter is used'; UNITS = 'm3 m-3'
    case('MWRTM_RGHWMAX');   LONG_NAME = 'L-band RTM model: Soil moisture value above which minimum microwave roughness parameter is used'; UNITS = 'm3 m-3'
    case('MWRTM_RGHNRH');    LONG_NAME = 'L-band RTM model: H-pol. Exponent for rough reflectivity parameterization'; UNITS = '1'
    case('MWRTM_RGHNRV');    LONG_NAME = 'L-band RTM model: V-pol. Exponent for rough reflectivity parameterization'; UNITS = '1'
    case('MWRTM_RGHPOLMIX'); LONG_NAME = 'L-band RTM model: Polarization mixing parameter'; UNITS = '1'
    case('MWRTM_OMEGA');     LONG_NAME = 'L-band RTM model: Scattering albedo';             UNITS = '1'
    case('MWRTM_BH');        LONG_NAME = 'L-band RTM model: H-pol. Vegetation b parameter'; UNITS = '1'
    case('MWRTM_BV');        LONG_NAME = 'L-band RTM model: V-pol. Vegetation b parameter'; UNITS = '1'
    case('MWRTM_LEWT');      LONG_NAME = 'L-band RTM model: Parameter to transform leaf area index into vegetation water content'; UNITS = 'kg m-2'

    case default;        LONG_NAME = 'Checck_GridComp';                                                  UNITS = 'Checck_GridComp';
    end select

    if (present(LNAME)) str_atr = trim (LONG_NAME)
    if (present(UNT))   str_atr = trim (UNITS    )

  END FUNCTION getAttribute

END PROGRAM mwrtm_bin2nc4
