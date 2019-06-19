MODULE LDAS_HashTable

  IMPLICIT NONE ! Use strong typing
  private
  INTEGER, PARAMETER :: tbl_size = 50

  TYPE nodelist
     TYPE(nodelist), POINTER :: child => NULL()
     CHARACTER(len=:), ALLOCATABLE     :: key
     integer                           :: fid
   CONTAINS
     PROCEDURE :: put  => put_nodeinfo
     PROCEDURE :: get  => get_nodeinfo
     PROCEDURE :: free => free_nodeinfo
  END TYPE nodelist

  TYPE hash_table
     TYPE(nodelist), DIMENSION(:), ALLOCATABLE :: vec
     INTEGER                                 :: vec_len = 0
     LOGICAL                                 :: is_init = .FALSE.
     integer                                 :: max_file_N = 20
     integer                                 :: counter = 0
   CONTAINS
     PROCEDURE :: init => init_hash_table
     PROCEDURE :: put  => put_hash_table
     PROCEDURE :: get  => get_hash_table
     PROCEDURE :: free => free_hash_table
  END TYPE hash_table
  PUBLIC :: nodelist
  PUBLIC :: hash_table
CONTAINS

  RECURSIVE SUBROUTINE put_nodeinfo(list,key,fid)
    CLASS(nodelist),target, INTENT(inout) :: list
    CHARACTER(len=*), INTENT(in)    :: key
    integer, INTENT(in)             :: fid
    ! local
    INTEGER                         :: keylen

    keylen = LEN(key)
    IF (ALLOCATED(list%key)) THEN
       IF (list%key /= key) THEN
          IF ( .NOT. ASSOCIATED(list%child) ) then
             ALLOCATE(list%child)
          ENDIF
          CALL put_nodeinfo(list%child,key,fid)
       END IF
    ELSE
       IF (.NOT. ALLOCATED(list%key)) &
            ALLOCATE(CHARACTER(len=keylen) :: list%key)

       list%key = key 
       list%fid = fid

    END IF
  END SUBROUTINE put_nodeinfo
  
  RECURSIVE SUBROUTINE get_nodeinfo(list,key,fid)
    CLASS(nodelist),   target, INTENT(in)    :: list
    CHARACTER(len=*),          INTENT(in)    :: key
    integer,                   INTENT(out)   :: fid

    IF (ALLOCATED(list%key) .AND. (list%key == key)) THEN
       fid = list%fid
    ELSE IF(ASSOCIATED(list%child)) THEN ! keep going
       CALL get_nodeinfo(list%child,key,fid)
    ELSE ! At the end of the list, no key found
       fid = -9999
       RETURN
    END IF
  END SUBROUTINE get_nodeinfo

  RECURSIVE SUBROUTINE free_nodeinfo(list,closefile)
    CLASS(nodelist), INTENT(inout) :: list
    external :: closefile
    integer  :: rc
    IF (ASSOCIATED(list%child)) THEN
       CALL free_nodeinfo(list%child, closefile )
       DEALLOCATE(list%child)
    END IF
    list%child => NULL()
    if (list%fid > 0) call closefile(list%fid)
    IF (ALLOCATED(list%key)) then
        DEALLOCATE(list%key)
    ENDIF

  END SUBROUTINE free_nodeinfo

  SUBROUTINE init_hash_table(tbl,tbl_len)
    CLASS(hash_table),   INTENT(inout) :: tbl
    INTEGER,     OPTIONAL, INTENT(in)    :: tbl_len

    if(tbl%is_init) return

    IF (allocated(tbl%vec)) DEALLOCATE(tbl%vec)
    IF (PRESENT(tbl_len)) THEN
       ALLOCATE(tbl%vec(0:tbl_len-1))
       tbl%vec_len = tbl_len
    ELSE
       ALLOCATE(tbl%vec(0:tbl_size-1))
       tbl%vec_len = tbl_size
    END IF
    tbl%is_init = .TRUE.
  END SUBROUTINE init_hash_table

  ! The first part of the hashing procedure using the string
  ! collating sequence
  ELEMENTAL FUNCTION sum_string(str) RESULT(sig)
    CHARACTER(len=*), INTENT(in)   :: str
    INTEGER                        :: sig
    CHARACTER, DIMENSION(LEN(str)) :: tmp
    INTEGER :: i

    FORALL (i=1:LEN(str))
       tmp(i) = str(i:i)
    END FORALL
    sig = SUM(ICHAR(tmp))
  END FUNCTION sum_string

  SUBROUTINE put_hash_table(tbl,key,fid)
    CLASS(hash_table), INTENT(inout) :: tbl
    CHARACTER(len=*),    INTENT(in)    :: key
    integer, INTENT(in)             :: fid
    !local
    INTEGER                            :: hash

    hash = MOD(sum_string(key),tbl%vec_len)
    CALL tbl%vec(hash)%put(key,fid)
    tbl%counter = tbl%counter + 1
  END SUBROUTINE put_hash_table


  SUBROUTINE get_hash_table(tbl,key,fid)
    CLASS(hash_table),INTENT(in)    :: tbl
    CHARACTER(len=*), INTENT(in)    :: key
    integer,          INTENT(out)   :: fid

    ! local
    INTEGER :: hash

    hash = MOD(sum_string(key),tbl%vec_len)
    CALL tbl%vec(hash)%get(key,fid)
  END SUBROUTINE get_hash_table


  SUBROUTINE free_hash_table(tbl,closefile,forced)
    CLASS(hash_table), INTENT(inout) :: tbl    
    INTEGER     :: i, low, high
    external :: closefile
    logical,intent(in) :: forced ! force to clean up

    if((.not. forced) .and. (tbl%counter < tbl%max_file_N)) return

    low  = LBOUND(tbl%vec,dim=1)
    high = UBOUND(tbl%vec,dim=1) 
    IF (allocated(tbl%vec)) THEN
       DO i=low,high
          CALL tbl%vec(i)%free(closefile)
       END DO
       DEALLOCATE(tbl%vec)
    END IF

    tbl%vec_len = 0
    tbl%counter = 0
    tbl%is_init = .FALSE.

  END SUBROUTINE free_hash_table

END MODULE LDAS_HashTable
