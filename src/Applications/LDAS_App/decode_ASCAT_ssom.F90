program decode_ASCAT_ssom

    implicit none
    
    real*8, dimension(15) :: tmp_vdata
    real*8, dimension(:,:), allocatable :: tmp_data

    integer, parameter :: lnbufr = 50
    integer, parameter :: max_obs = 250000
    integer :: idate,iret
    integer :: ireadmg,ireadsb
    integer :: N_obs

    character(8)    :: subset
    character(300)  :: fname, mastertable_path

! -------------------------------------------------------------------------

    fname = '/home/amfox/smap/SMAP_Nature/ASCAT_EUMETSAT/Metop_C/Y2023/M03/' // &
            'M03-ASCA-ASCSMO02-NA-5.0-20230301090900.000000000Z-20230301105557-4839070.bfr'
    mastertable_path = '/home/amfox/smap/SMAP_Nature/ASCAT_EUMETSAT'

!   Allocate the tmp_data array
    allocate(tmp_data(max_obs, 15))

    open(lnbufr, file=trim(fname), action='read',form='unformatted')

    call openbf(lnbufr,'SEC3', lnbufr)
    call mtinfo( trim(mastertable_path) // '/BUFR_mastertable/', 51, 52)
    call datelen(10)

    N_obs = 0        
    msg_report: do while(ireadmg(lnbufr,subset,idate) ==0)
        loop_report: do while(ireadsb(lnbufr) == 0)
            call ufbint(lnbufr,tmp_vdata,15,1,iret,'YEAR MNTH DAYS HOUR MINU SECO SSOM DOMO SMPF SMCF ALFR TPCX IWFR CLATH CLONH')
            N_obs = N_obs + 1
            tmp_data(N_obs,:) = tmp_vdata
        end do loop_report
    end do msg_report

    write(*,*) 'N_obs = ', N_obs
    write(*,*) tmp_vdata

    call closbf(lnbufr)
    close(lnbufr)

end program decode_ASCAT_ssom
