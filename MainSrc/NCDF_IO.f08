MODULE NCDF_IO

  USE netcdf

  IMPLICIT NONE

CONTAINS

!==================================================================================================================================!
!Subroutine HANDLE_ERR
!
! Checks for errors in NetCDF routines and aborts program if any are detected
!
! NOTES::
!
! WARNINGS::
!
! OUTPUTS::
!
! INPUTS::
!
!==================================================================================================================================!
SUBROUTINE HANDLE_ERR(Status,Location)
  INTEGER,INTENT(IN):: Status
  CHARACTER(*),INTENT(IN),OPTIONAL:: Location

  IF (Status .NE. nf90_noerr) THEN
    WRITE(*,*) '***   NETCDF ERROR ENCOUNTERED   ***'
    IF (PRESENT(Location)) WRITE(*,*) '***   Error occured in ',Location
    WRITE(*,*) '***   ',Trim(nf90_strerror(Status))
    STOP 'Program Aborted'
  END IF

END SUBROUTINE HANDLE_ERR

!==================================================================================================================================!
!Subroutine NF_OPEN_FILE
!
! Opens a NetCDF dataset - can open an old dataset or create a new one
!
! NOTES::
!   If a new file is created (Ftype = 'new') then it will be opened in define mode with read/write access
!   New files are always created as an HDF5 type
!   If an old file is opened (Ftype = 'old') then it will be opened in data mode
!
! OUTPUTS::
!   ncID - the NetCDF dataset ID number (integer)
!
! INPUTS::
!   Fname - name of (or path to) the file/dataset
!   Ftype - specifies if the file to open is new or old: Ftype = 'new' will create a new file and overwrite any existing files ->
!           -> of the same name. Ftype = 'old' expects a file named Fname to already exist to open
!   rw - an optional input to declare whether to give read-only or read/write access to an opened file (only for Ftype = 'old')
!
!==================================================================================================================================!
SUBROUTINE NF_OPEN_FILE(ncID,Fname,Ftype,rw)
  INTEGER,INTENT(OUT):: ncID
  CHARACTER(*),INTENT(IN):: Fname, Ftype
  CHARACTER(*),INTENT(IN),OPTIONAL:: rw
  INTEGER:: Status, rw2
  CHARACTER(42):: Location = 'MODULE: NCDF_IO / SUBROUTINE: NF_OPEN_FILE'

  !--------------------------------------------------!
  !         New File -- call nf90_create             !
  !--------------------------------------------------!
  IF ( ANY(Ftype .EQ. (/'new','New','NEW'/)) ) THEN
    Status = nf90_create(Fname,NF90_NETCDF4,ncID)

  !--------------------------------------------------!
  !         Old File -- call nf90_open               !
  !--------------------------------------------------!
  ELSE IF ( ANY(Ftype .EQ. (/'old','Old','OLD'/)) ) THEN
    IF (PRESENT(rw)) THEN !if user input rw, check
      IF ( ANY(rw .EQ. (/'r   ','R   ','read','Read','READ'/)) ) THEN !read only access
        rw2 = NF90_NOWRITE

      ELSE IF ( ANY(rw .EQ. (/'w    ','W    ','rw   ','RW   ','write','Write','WRITE'/)) ) THEN !read/write access
        rw2 = NF90_WRITE

      ELSE !unrecognized rw
        WRITE(*,*) 'Bad rw detected in NF_OPEN_FILE (',rw,')'
        STOP

      END IF

    ELSE !if no rw input, default to read/write access
      rw2 = NF90_WRITE

    END IF

    Status = nf90_open(Fname,rw2,ncID)

  !--------------------------------------------------!
  !         Unrecognized Ftype -- abort              !
  !--------------------------------------------------!
  ELSE
    WRITE(*,*) 'Bad Ftype detected in NF_OPEN_FILE (',Ftype,')'
    STOP

  END IF

  !Checking for any errors that may have occured
  CALL HANDLE_ERR(Status,Location)

END SUBROUTINE NF_OPEN_FILE

!==================================================================================================================================!
!Subroutine NF_CLOSE_FILE
!
! Closes an open NetCDF dataset
!
!==================================================================================================================================!
SUBROUTINE NF_CLOSE_FILE(ncID)
  INTEGER,INTENT(IN):: ncID
  INTEGER:: Status
  CHARACTER(43):: Location = 'MODULE: NCDF_IO / SUBROUTINE: NF_CLOSE_FILE'

  !calling NetCDF function
  Status = nf90_close(ncID)

  !Checking for any errors that may have occured
  CALL HANDLE_ERR(Status,Location)

END SUBROUTINE NF_CLOSE_FILE

!==================================================================================================================================!
!Subroutine NF_REDEF
!
! Puts a NetCDF dataset into define mode
!
!==================================================================================================================================!
SUBROUTINE NF_REDEF(ncID)
  INTEGER,INTENT(IN):: ncID
  INTEGER:: Status
  CHARACTER(38):: Location = 'MODULE: NCDF_IO / SUBROUTINE: NF_REDEF'

  !calling NetCDF function
  Status = nf90_redef(ncID)

  !Checking for any errors that may have occured
  CALL HANDLE_ERR(Status,Location)

END SUBROUTINE NF_REDEF

!==================================================================================================================================!
!Subroutine NF_ENDDEF
!
! Puts a NetCDF dataset into data mode
!
! NOTES::
!   Changes made previously to the file in question while in define mode will be committed to disk
!
!==================================================================================================================================!
SUBROUTINE NF_ENDDEF(ncID)
  INTEGER,INTENT(IN):: ncID
  INTEGER:: Status
  CHARACTER(39):: Location = 'MODULE: NCDF_IO / SUBROUTINE: NF_ENDDEF'

  !calling NetCDF function
  Status = nf90_enddef(ncID)

  !Checking for any errors that may have occured
  CALL HANDLE_ERR(Status,Location)

END SUBROUTINE NF_ENDDEF

!==================================================================================================================================!
!Subroutine NF_DEF_DIM
!
! Defines a new dimension in a NetCDF dataset
!
! WARNINGS::
!   This routine should only be called to act on a NetCDF dataset that is in define mode
!   Attempting to use this routine on a NetCDF dataset in data mode will fail and abort the program
!
! OUTPUTS::
!   dimID(int) - the dimension ID number
!
! INPUTS::
!   ncID(int) - the NetCDF dataset ID number
!   Dname(str) - name of the dimension
!   Len(int) - length of the dimension, Len=0 defines an unlimited length dimension
!
!==================================================================================================================================!
SUBROUTINE NF_DEF_DIM(dimID,ncID,Dname,Len)
  INTEGER,INTENT(OUT):: dimID
  INTEGER,INTENT(IN):: ncID, Len
  CHARACTER(*),INTENT(IN):: Dname
  INTEGER:: Status, Len2
  CHARACTER(40):: Location = 'MODULE: NCDF_IO / SUBROUTINE: NF_DEF_DIM'

  !Determining if defining an unlimited dimension
  IF (Len .EQ. 0) THEN !if len=0, dimension has unlimited length
    Len2 = NF90_UNLIMITED
  ELSE !if len/=0, dimension has set length
    Len2 = Len
  END IF

  !calling NetCDF function
  Status = nf90_def_dim(ncID,Dname,Len2,dimID)

  !Checking for any errors that may have occured
  CALL HANDLE_ERR(Status,Location)

END SUBROUTINE NF_DEF_DIM

!==================================================================================================================================!
!Subroutine NF_DEF_VAR
!
! Defines a new variable in a NetCDF dataset
!
! NOTES::
!   To create a scalar variable, dimID must be an array of size 0
!
! WARNINGS::
!   This routine should only be called to act on a NetCDF dataset that is in define mode
!   Attempting to use this routine on a NetCDF dataset in data mode will fail and abort the program
!
! OUTPUTS::
!   varID(int) - the variable ID number
!
! INPUTS::
!   ncID(int) - the NetCDF dataset ID number
!   dimID(int) - an array of the dimension ID numbers corresponding to the dimensions of the variable
!   Vname(str) - name of the variable
!   Vtype(str) - the type of the variable
!
!==================================================================================================================================!
SUBROUTINE NF_DEF_VAR(varID,ncID,dimID,Vname,Vtype)
  INTEGER,INTENT(OUT):: varID
  INTEGER,INTENT(IN):: ncID, dimID(:)
  CHARACTER(*),INTENT(IN):: Vname, Vtype
  INTEGER:: Status, type
  CHARACTER(40):: Location = 'MODULE: NCDF_IO / SUBROUTINE: NF_DEF_VAR'

  IF ( ANY(Vtype .EQ. (/'double','Double','DOUBLE'/)) ) THEN
    type = NF90_DOUBLE
  ELSE IF ( ANY(Vtype .EQ. (/'real','Real','REAL'/)) ) THEN
    type = NF90_FLOAT
  ELSE IF ( ANY(Vtype .EQ. (/'int','Int','INT'/)) ) THEN
    type = NF90_INT
  ELSE IF ( ANY(Vtype .EQ. (/'str   ','Str   ','STR   ','string','String','STRING'/)) ) THEN
    type = NF90_CHAR
  ELSE
    WRITE(*,*) 'Unrecognized Vtype detected in ',Location,' (',Vtype,')'
    STOP
  END IF

  !calling NetCDF function
  Status = nf90_def_var(ncID,Vname,type,dimID,varID)

  !Checking for any errors that may have occured
  CALL HANDLE_ERR(Status,Location)

END SUBROUTINE NF_DEF_VAR

!==================================================================================================================================!
!Subroutine NF_DEF_UNIT
!
! Defines a new variable in a NetCDF dataset
!
! NOTES::
!
! WARNINGS::
!   This routine should only be called to act on a NetCDF dataset that is in define mode
!   Attempting to use this routine on a NetCDF dataset in data mode will fail and abort the program
!
! OUTPUTS::
!
! INPUTS::
!   ncID(int) - the NetCDF dataset ID number
!   varID(int) - the ID number of the variable being given units
!   Unit(str) - the units
!
!==================================================================================================================================!
SUBROUTINE NF_DEF_UNIT(ncID,varID,Unit)
  INTEGER,INTENT(IN):: ncID, varID
  CHARACTER(*),INTENT(IN):: Unit
  INTEGER:: Status
  CHARACTER(41):: Location = 'MODULE: NCDF_IO / SUBROUTINE: NF_DEF_UNIT'

  !calling NetCDF function
  Status = nf90_put_att(ncID,varID,'Units',Unit)

  !Checking for any errors that may have occured
  CALL HANDLE_ERR(Status,Location)

END SUBROUTINE NF_DEF_UNIT

!==================================================================================================================================!
!
!==================================================================================================================================!

END MODULE NCDF_IO
