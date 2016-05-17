PROGRAM GEOSChem_Helen_he5_create

  USE ISO_C_BINDING
  IMPLICIT NONE

  ! -----------------------------------
  ! Dimensions and values of AMF arrays
  ! -----------------------------------
  INTEGER(KIND=C_LONG), PARAMETER :: nlon = 576, nlat = 361, nlev = 47
  
  ! --------------------
  ! HE5 Output variables
  ! --------------------
  CHARACTER (LEN=132)          :: he5file
  CHARACTER (LEN=3), PARAMETER :: sensor = 'OMI'
  CHARACTER (LEN=4), PARAMETER :: year   = '2007'

  ! ----------------------------
  ! Set the HE5 output file name
  ! ----------------------------
  he5file = TRIM(ADJUSTL(sensor))//'_H2O_MERRA.he5'

  ! ------------------------------------------------------
  ! Initialize output quantity and call HE5 output routine
  ! ------------------------------------------------------
  CALL create_he5_file ( he5file, TRIM(ADJUSTL(sensor)), TRIM(ADJUSTL(year)), nlon, nlat, nlev )

  STOP
END PROGRAM GEOSChem_Helen_he5_create

SUBROUTINE create_he5_file ( he5file, sensor, year, nlon, nlat, nlev )

  !------------------------------------------------------------------------------
  ! This subroutine creates an HE5 file with BrO Wavelength-Dependent AMFs
  !
  ! Input:
  !   he5file    - Name of HE5 output file
  !   sensor     - Name of satellite sensor
  !   nlon       - Number of longitudes
  !   nlat       - Number of latitudes
  !   nlev       - Number of vertical levels
  !------------------------------------------------------------------------------

  USE GEOSChem_Helen_he5_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  CHARACTER (LEN=*),    INTENT (IN) :: he5file, sensor, year
  INTEGER(KIND=C_LONG), INTENT (IN) :: nlon, nlat, nlev

  ! ---------------
  ! Local variables
  ! ---------------
  INTEGER(KIND=C_LONG)                       :: nle1
  INTEGER                                    :: swath_file_id  
  ! ID number for HE5 output file (required for closing it)
  INTEGER                                    :: swath_id       
  ! ID number for swath (required for writing to swath)
  REAL (KIND=r4), DIMENSION (nlon)           :: longitudes
  REAL (KIND=r4), DIMENSION (nlat)           :: latitudes
  REAL (KIND=r8), DIMENSION (nlon,nlat,nlev) :: H2O8
  REAL (KIND=r8), DIMENSION (nlon,nlat,nlev) :: Temperature8
  REAL (KIND=r8), DIMENSION (nlon,nlat,nlev) :: Psurface8
  REAL (KIND=r4), DIMENSION (nlon,nlat,nlev) :: gcdatar4
  INTEGER :: he5stat, isw, idf
  
  CHARACTER (LEN=132) :: swath_name, ascii_name


  he5stat = 0

  ! ---------------------------
  ! Initialize output variables
  ! ---------------------------
  longitudes   = 0.0_r4
  latitudes    = 0.0_r4
  H2O8         = 0.0_r8
  Temperature8 = 0.0_r8
  Psurface8    = 0.0_r8
  gcdatar4     = 0.0_r4
  nle1         = nlev+1
  
  ! ---------------------------------------------------------------
  ! Open HE5 output file and check AMF_SWATH_FILE_ID ( -1 if error)
  ! ---------------------------------------------------------------
  swath_file_id = HE5_SWopen ( TRIM(ADJUSTL(he5file)), he5f_acc_trunc )
  IF ( swath_file_id == -1 ) THEN
     WRITE (*,*) 'ERROR: HE5_SWopen failed!'; STOP 1
  END IF

  ! ------------------------
  ! Loop over Monthly Swaths
  ! ------------------------
  swathloop: DO isw = 1, nmonth
     
     swath_name = TRIM(ADJUSTL(sensor))//'_Climatology_'//TRIM(ADJUSTL(SwathFieldsMonth(isw)))

     ! ------------------------------------------------------
     ! Create HE5 swath and check AMF_SWATH_ID ( -1 if error)
     ! ------------------------------------------------------
     swath_id = HE5_SWcreate ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
     IF ( swath_id == -1 ) CALL  he5_error_stop ( 'HE5_SWcreate', '<'//   &
                                 TRIM(ADJUSTL(swath_name))//'>' )
  
     ! ----------------------------------
     ! Define new dimensions in HE5 swath
     ! ----------------------------------
     he5stat = HE5_SWdefdim ( swath_id, nLonDim, nlon ) 
     IF ( he5stat /= 0 ) CALL he5_error_stop ( 'HE5_SWdefdim', nLonDim )
     he5stat = HE5_SWdefdim ( swath_id, nLatDim, nlat )
     IF ( he5stat /= 0 ) CALL he5_error_stop ( 'HE5_SWdefdim', nLatDim )
     he5stat = HE5_SWdefdim ( swath_id, nETADim, nlev )
     IF ( he5stat /= 0 ) CALL he5_error_stop ( 'HE5_SWdefdim', nETADim )
     he5stat = HE5_SWdefdim ( swath_id, nETADp1, nle1 )
     IF ( he5stat /= 0 ) CALL he5_error_stop ( 'HE5_SWdefdim', nETADp1 )

     ! --------------------------
     ! Define Geolocation Fields
     ! --------------------------
     defgeofields: DO idf = 1, ngf
        he5stat = HE5_SWdefgfld ( swath_id, TRIM(ADJUSTL(GeoFieldNames(idf))), &
                                  TRIM(ADJUSTL(GeoFieldDims(idf))), " ",       &
                                  HE5T_Native_FLOAT, he5_hdfe_nomerge )
        IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWdefgfld', &
                                  TRIM(ADJUSTL(GeoFieldNames(idf))) )
     END DO defgeofields

     ! ------------------
     ! Define Data Fields
     ! ------------------
     defdatafields: DO idf = 1, ndf
        CALL define_data_fields ( nlon, nlat, nlev, idf, swath_id )
     END DO defdatafields
  
     ! -------------------------------------------------------------------------------
     ! Detach from and re-attach to created swath (recommended before adding to swath)
     ! -------------------------------------------------------------------------------
     he5stat   = HE5_SWdetach ( swath_id )
     swath_id  = HE5_SWattach ( swath_file_id, TRIM(ADJUSTL(swath_name)) )
     IF ( swath_id == -1 ) CALL  he5_error_stop ( 'HE5_SWattach', TRIM(ADJUSTL(swath_name)) )

     ! ---------------------------
     ! Read and Write Data Fields
     ! ---------------------------
     ascii_name = &
          './Helen_MERRA2/'//'ctm.'// &
          TRIM(ADJUSTL(year))//''//         &
          ShortNameMonths(isw)//'.dat'

     PRINT *, 'reading from file '//TRIM(ADJUSTL(ascii_name))//' & writing to '// &
              TRIM(ADJUSTL(he5file))
     CALL read_gchem_ascii_file  ( &
          TRIM(ADJUSTL(ascii_name)), nlon, nlat, nlev, &
          longitudes(1:nlon), latitudes(1:nlat),               &
          H2O8(1:nlon,1:nlat,1:nlev), Temperature8(1:nlon,1:nlat,1:nlev),         &
          Psurface8(1:nlon,1:nlat,1:nlev))
     
     wrtdatafields: DO idf = 1, ndf
        gcdatar4 = 0.0_r4
        SELECT CASE (TRIM(ADJUSTL(DataFieldNames(idf))))
        CASE ("VMRH2O             ")
           gcdatar4(1:nlon,1:nlat,1:nlev) = &
                REAL( H2O8, KIND=r4 ) ! VMR
        CASE ("SurfacePressure    ")
           gcdatar4(1:nlon,1:nlat,1:nlev) = &
                REAL( Psurface8, KIND=r4 ) ! hPa
        CASE ("TemperatureProfile ")
           gcdatar4(1:nlon,1:nlat,1:nlev) = REAL( Temperature8, KIND=r4 )
        CASE DEFAULT
           WRITE(*,*) 'Warning!! Field '//TRIM(ADJUSTL(DataFieldNames(idf)))//' not found.' 
        END SELECT        
        gcdatar4 = gcdatar4 / ScaleFactors(idf)
        CALL write_data_fields ( nlon, nlat, nlev, idf, swath_id, &
        gcdatar4(1:nlon,1:nlat,1:nlev) )

        CALL write_field_attributes ( 'data', isw, idf, swath_id )
     END DO wrtdatafields

     ! --------------------------
     ! Write Geolocation Fields
     ! --------------------------
     CALL write_geo_field ( nlon, swath_id, 'Longitudes', longitudes(1:nlon) )
     CALL write_geo_field ( nlat, swath_id, 'Latitudes' , latitudes (1:nlat) )
     wrtgeofields: DO idf = 1, ngf
        CALL write_field_attributes ( 'geo', isw, idf, swath_id )
     END DO wrtgeofields

     ! -----------------------
     ! Write Global Attributes
     ! -----------------------
     CALL write_swath_attributes ( swath_id )

     ! -----------------------------------------------
     ! Detach from HE5 swath
     ! -----------------------------------------------
     he5stat = HE5_SWdetach ( swath_id )
     IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWdetach', TRIM(ADJUSTL(swath_name)) )
  
  END DO swathloop

  ! -----------------------------------------------
  ! Close HE5 output file
  ! -----------------------------------------------
  he5stat = HE5_SWclose  ( swath_file_id )
  IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWclose', TRIM(ADJUSTL(he5file)) )

  
  RETURN
END SUBROUTINE create_he5_file 

SUBROUTINE  read_gchem_ascii_file  ( ascii_name, nlon, nlat, nlev, &
                                     longitudes, latitudes,        &
                                     H2O, Temperature, Psurface )

  USE GEOSChem_Helen_he5_module
  IMPLICIT NONE

  ! ---------------
  ! Input variables
  ! ---------------
  INTEGER(KIND=C_LONG), INTENT (IN) :: nlon, nlat, nlev
  CHARACTER (LEN=*),    INTENT (IN) :: ascii_name
  ! ---------------
  ! Output variables
  ! ---------------
  REAL (KIND=r4), INTENT (OUT), DIMENSION (nlon)           :: longitudes
  REAL (KIND=r4), INTENT (OUT), DIMENSION (nlat)           :: latitudes
  REAL (KIND=r8), INTENT (OUT), DIMENSION (nlon,nlat)      :: Psurface
  REAL (KIND=r8), INTENT (OUT), DIMENSION (nlon,nlat,nlev) :: H2O, Temperature
  ! ------------------------------
  ! Local variables and parameters
  ! ------------------------------
  INTEGER, PARAMETER :: aunit = 99
  INTEGER            :: ilon, ilat, ios

  OPEN (UNIT=aunit, FILE=TRIM(ADJUSTL(ascii_name)), STATUS='OLD', ACTION='READ', IOSTAT=ios)
  IF ( ios /= 0 ) THEN
     WRITE (*, '(A)') 'ERROR reading from ASCII file '//TRIM(ADJUSTL(ascii_name))
     STOP 1
  END IF

  DO ilon = 1, nlon
     DO ilat = 1, nlat
        READ (UNIT=aunit, FMT=*) longitudes(ilon), latitudes(ilat), Psurface(ilon,ilat)
        READ (UNIT=aunit, FMT=*) H2O(ilon,ilat,1:nlev)
        READ (UNIT=aunit, FMT=*) Temperature(ilon,ilat,1:nlev)
     END DO
  END DO
  CLOSE (aunit)

  RETURN
END SUBROUTINE read_gchem_ascii_file
