MODULE GEOSChem_Helen_he5_module

  USE ISO_C_BINDING

  IMPLICIT NONE
  INCLUDE 'hdfeos5.inc'


  ! --------------
  ! Precision KIND
  ! --------------
  INTEGER, PARAMETER :: r8 = KIND(1.0D0)
  INTEGER, PARAMETER :: r4 = KIND(1.0  )

  ! ------------
  ! HE5 Routines
  ! ------------
  INTEGER (KIND = 4), EXTERNAL ::                                             &
       he5_swattach, he5_swclose, he5_swcreate, he5_swdefdfld, he5_swdefgfld, &
       he5_swdefdim, he5_swdetach,  he5_swopen, he5_swwrfld,  he5_swwrlattr,  &
       he5_swdefcomch, he5_swwrattr

  ! ----------------------------------
  ! Start, Stride, and Endge Variables
  ! ----------------------------------
  INTEGER(C_LONG)               :: he5_start_1, he5_stride_1, he5_edge_1
  INTEGER(C_LONG), DIMENSION(2) :: he5_start_2, he5_stride_2, he5_edge_2
  INTEGER(C_LONG), DIMENSION(3) :: he5_start_3, he5_stride_3, he5_edge_3
  
  ! -----------------------------------------------------------------------
  ! Parameters and variables associated with field compression and chunking
  ! -----------------------------------------------------------------------
  INTEGER(C_LONG), DIMENSION (4) :: cchunk_dim
  INTEGER                        :: comp_par, comp_type, ncchunkdim

  ! ----------------------------------------
  ! General Swath Fields, one for each month
  ! ----------------------------------------
  INTEGER,                               PARAMETER :: nmonth = 12
  !CHARACTER (LEN=14), DIMENSION (nmonth), PARAMETER :: SwathFieldsMonth = (/    &
  !"01 (January)  ", "02 (February) ", "03 (March)    ", "04 (April)    ", "05 (May)      ", "06 (June)     ",&
  !"07 (July)     ", "08 (August)   ", "09 (September)", "10 (October)  ", "11 (November) ", "12 (December) " /)

  CHARACTER (LEN=12), DIMENSION (nmonth), PARAMETER :: SwathFieldsMonth = (/    &
  "01_January  ", "02_February ", "03_March    ", "04_April    ", "05_May      ", "06_June     ",&
  "07_July     ", "08_August   ", "09_September", "10_October  ", "11_November ", "12_December " /)

  CHARACTER (LEN=2), DIMENSION (nmonth), PARAMETER :: ShortNameMonths = (/ &
       "01",  "02",  "03",  "04",  "05",  "06",  "07",  "08",  "09", "10",  "11",  "12"  /)

  ! ------------------------------
  ! Names of Dimensions in Swaths
  ! ------------------------------
  CHARACTER (LEN= 4), PARAMETER :: nLonDim  = 'nLon'
  CHARACTER (LEN= 4), PARAMETER :: nLatDim  = 'nLat'
  CHARACTER (LEN= 4), PARAMETER :: nETADim  = 'nETA'
  CHARACTER (LEN= 4), PARAMETER :: nETADp1  = 'nEp1'
  CHARACTER (LEN= 4), PARAMETER :: n1Dim = 'nETA'
  CHARACTER (LEN= 9), PARAMETER :: n2Dim = 'nLon,nLat'
  CHARACTER (LEN=14), PARAMETER :: n3Dim = 'nLon,nLat,nETA'

  ! ----------------------------------------------
  ! Number of data and geolocation fields in swath
  ! ----------------------------------------------
  INTEGER, PARAMETER :: ndf = 3, ngf = 2

  ! ------------------------------
  ! Names of Data Fields in Swaths
  ! ------------------------------
  CHARACTER (LEN=10), DIMENSION (ngf), PARAMETER :: GeoFieldNames = (/ "Longitudes",  "Latitudes " /)
  CHARACTER (LEN= 4), DIMENSION (ngf), PARAMETER :: GeoFieldDims  = (/ "nLon",        "nLat"       /)


  ! ---------------------------------------
  ! Titles of Geolocation Fields in Swaths
  ! --------------------------------------
  CHARACTER (LEN=40), DIMENSION (ngf), PARAMETER :: GeoFieldTitles = (/ &
       "Geodetic Longitude at Center of Grid Box",  &
       "Geodetic Latitude at Center of Grid Box "    /)

  ! -------------------------------------
  ! Units of Geolocation Fields in Swaths
  ! -------------------------------------
  CHARACTER (LEN=8), DIMENSION (ngf), PARAMETER :: GeoFieldUnits = (/ &
       "deg     ",  &
       "deg     "   /) 


  ! -----------------------------------------------
  ! Short-Names of Data Fields (part of file names)
  ! -----------------------------------------------
  CHARACTER (LEN=4), DIMENSION (ndf), PARAMETER :: ShortDataFieldNames = (/ &
       "H2O_", &
       "PSrf", &
       "Temp"   /)

  ! ------------------------------
  ! Names of Data Fields in Swaths
  ! ------------------------------
  CHARACTER (LEN=19), DIMENSION (ndf), PARAMETER :: DataFieldNames = (/ &
       "VMRH2O             ",  &
       "SurfacePressure    ",  &
       "TemperatureProfile "    /)

  ! ------------------------------
  ! Units of Data Fields in Swaths
  ! ------------------------------
  CHARACTER (LEN=3), DIMENSION (ndf), PARAMETER :: DataFieldUnits = (/ &
       "ppb",  &
       "hPa",  &
       "K  "    /) 
 
  ! ------------------------------
  ! Units of Data Fields in Swaths
  ! ------------------------------
  REAL (KIND=r4), DIMENSION (ndf), PARAMETER :: ScaleFactors = (/ &
       1.0E+04_r4,  &
       1.0E-00_r4,  &
       1.0E-00_r4   /) 
 
  ! ------------------------------
  ! Titles of Data Fields in Swaths
  ! ------------------------------
  CHARACTER (LEN=25), DIMENSION (ndf), PARAMETER :: DataFieldTitles = (/ &
       "Water Volume Mixing Ratio",  &
       "Surface Pressure         ",  &
       "Temperature Profile      "   /)

  ! ------------------------------
  ! Names of Data Provider
  ! ------------------------------
  CHARACTER (LEN=10), DIMENSION (ndf), PARAMETER :: DataProvider = (/ &
       "Helen Wang",   &
       "Helen Wang",   &
       "Helen Wang"   /)

  ! ----------------------------------------------------
  ! Temporal Range of GChem Simulation (MERRA)
  ! ----------------------------------------------------
  CHARACTER (LEN=12), DIMENSION (ndf) :: TimeOfDay = (/ &
       "1300h-1400h",  &
       "1300h-1400h",  &
       "1300h-1400h"    /)

  ! ---------------------------
  ! Version of GChem Simulation
  ! ---------------------------
  CHARACTER (LEN=2), DIMENSION (ndf), PARAMETER :: MERRAVersion = (/ &
       "02",  &
       "02",  &
       "02"    /)

  ! -----------------------------------
  ! Data of production GChem Simulation
  ! -----------------------------------
  CHARACTER (LEN=11), DIMENSION (ndf), PARAMETER :: DateOfProduction = (/ &
       "2015-06-30",  &
       "2015-06-30",  &
       "2015-06-30"    /)


CONTAINS

  SUBROUTINE define_data_fields ( nlon, nlat, nlev, ifield, swath_id )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,               INTENT (IN) :: ifield, swath_id
    INTEGER (KIND=C_LONG), INTENT (IN) :: nlon, nlat, nlev

    ! --------------
    ! Local variable
    ! --------------
    INTEGER         :: he5stat

    ! -----------------------
    ! Define compression type
    ! -----------------------
    comp_par    = 9 ; comp_type = HE5_HDFE_COMP_SHUF_DEFLATE

    ! ---------------------------------------------------------------------
    ! Define array value data fields in HE5 swath. SurfacePressure is the
    ! only 2D variable, everything else is 3D/
    ! ---------------------------------------------------------------------
    IF ( TRIM(ADJUSTL(DataFieldNames(ifield))) == 'SurfacePressure' ) THEN
       !ncchunkdim  = 2 ; cchunk_dim(1:ncchunkdim) = (/ nlon, nlat /)
       !he5stat = HE5_SWdefcomch ( &
       !     swath_id, comp_type, comp_par, ncchunkdim, cchunk_dim(1:ncchunkdim) )
       !IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWdefcomch', TRIM(ADJUSTL(DataFieldNames(ifield))) )

       he5stat = HE5_SWdefdfld ( swath_id, &
            TRIM(ADJUSTL(DataFieldNames(ifield))),  n2Dim, " ", HE5T_Native_FLOAT, he5_hdfe_nomerge )

    ELSE
       !ncchunkdim  = 3 ; cchunk_dim(1:ncchunkdim) = (/ nlon, nlat, nlev /)
       !he5stat = HE5_SWdefcomch ( &
       !     swath_id, comp_type, comp_par, ncchunkdim, cchunk_dim(1:ncchunkdim) )
       !IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWdefcomch', TRIM(ADJUSTL(DataFieldNames(ifield))) )

       he5stat = HE5_SWdefdfld ( swath_id, &
            TRIM(ADJUSTL(DataFieldNames(ifield))),  n3Dim, " ", HE5T_Native_FLOAT, he5_hdfe_nomerge )
    END IF
    IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWdefdfld', TRIM(ADJUSTL(DataFieldNames(ifield))) )

    RETURN
  END SUBROUTINE define_data_fields

  SUBROUTINE write_data_fields ( nlon, nlat, nlev, ifield, swath_id, gcdata)

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,                                       INTENT (IN) :: ifield, swath_id
    INTEGER (KIND=C_LONG),                         INTENT (IN) :: nlon, nlat, nlev
    REAL    (KIND=r4), DIMENSION (nlon,nlat,nlev), INTENT (IN) :: gcdata

    ! --------------
    ! Local variable
    ! --------------
    INTEGER :: he5stat

    ! ---------------------------------------------------------------------
    ! Write array value data fields in HE5 swath. SurfacePressure is the
    ! only 2D variable, everything else is 3D/
    ! ---------------------------------------------------------------------
    IF ( TRIM(ADJUSTL(DataFieldNames(ifield))) == 'SurfacePressure' ) THEN
       he5_start_2  = (/    0,    0 /)
       he5_stride_2 = (/    1,    1 /)
       he5_edge_2   = (/ nlon, nlat /)

       he5stat = HE5_SWwrfld ( swath_id, &
            TRIM(ADJUSTL(DataFieldNames(ifield))), he5_start_2, he5_stride_2, he5_edge_2, &
            gcdata(1:nlon,1:nlat,1) )
    ELSE

       he5_start_3  = (/    0,    0,    0 /)
       he5_stride_3 = (/    1,    1,    1 /)
       he5_edge_3   = (/ nlon, nlat, nlev /)

       he5stat = HE5_SWwrfld ( swath_id, &
            TRIM(ADJUSTL(DataFieldNames(ifield))), he5_start_3, he5_stride_3, he5_edge_3, &
            gcdata(1:nlon,1:nlat,1:nlev) )
    END IF
    IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWwrfld', TRIM(ADJUSTL(DataFieldNames(ifield))) )


    RETURN
  END SUBROUTINE write_data_fields

  SUBROUTINE write_geo_field ( ndim, swath_id, geofield, gfdata )

    IMPLICIT NONE

    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,                             INTENT (IN) :: swath_id
    INTEGER (KIND=C_LONG),               INTENT (IN) :: ndim
    REAL    (KIND=r4), DIMENSION (ndim), INTENT (IN) :: gfdata
    CHARACTER (LEN=*),                   INTENT (IN) :: geofield

    ! --------------
    ! Local variable
    ! --------------
    INTEGER :: he5stat

    he5_start_1  = 0 ; he5_stride_1 = 1 ; he5_edge_1   = ndim

    he5stat = HE5_SWwrfld ( swath_id, &
            TRIM(ADJUSTL(geofield)), he5_start_1, he5_stride_1, he5_edge_1, gfdata(1:ndim) )
    IF ( he5stat /= 0 ) CALL  he5_error_stop ( 'HE5_SWwrfld', TRIM(ADJUSTL(geofield)) )

    RETURN
  END SUBROUTINE write_geo_field

  SUBROUTINE write_field_attributes ( geodat, imonth, ifield, swath_id )
  
    IMPLICIT NONE
  
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER,           INTENT (IN) :: imonth, ifield, swath_id
    CHARACTER (LEN=*), INTENT (IN) :: geodat
    ! --------------
    ! Local variable
    ! --------------
    INTEGER(C_LONG) :: n_attr, n_one = 1
    INTEGER         :: he5stat
  
  
    ! -------------------------------------------------
    ! Attributes for data fields but geolocation fields
    ! -------------------------------------------------
    IF ( INDEX ( TRIM(ADJUSTL(geodat)), 'geo') > 0 ) THEN

       n_attr = LEN_TRIM(ADJUSTL(GeoFieldTitles(ifield)))
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(GeoFieldNames(ifield))), "Title", &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(GeoFieldTitles(ifield))) )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'Title')

       he5stat = HE5_SWwrlattr (                                            &
            swath_id, TRIM(ADJUSTL(GeoFieldNames(ifield))), "ScaleFactor",  &
            HE5T_NATIVE_FLOAT, n_one, 1.0E+00                                 )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'ScaleFactor')
  
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(GeoFieldNames(ifield))), "MissingValue", &
            HE5T_NATIVE_FLOAT, n_one, 1.0E-30                                 )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'MissingValue')
  
       n_attr = LEN_TRIM(ADJUSTL(GeoFieldUnits(ifield)))
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(GeoFieldNames(ifield))), "Units",         &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(GeoFieldUnits(ifield))) )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'Units')

    ELSE

       n_attr = LEN_TRIM(ADJUSTL(DataFieldTitles(ifield)))
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "Title", &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(DataFieldTitles(ifield))) )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'Title')

       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "ScaleFactor",  &
            HE5T_NATIVE_FLOAT, n_one, ScaleFactors(ifield)                    )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'ScaleFactor')
  
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "MissingValue", &
            HE5T_NATIVE_FLOAT, n_one, 1.0E-30                                 )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'MissingValue')
  
       n_attr = LEN_TRIM(ADJUSTL(DataFieldUnits(ifield)))
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "Units", &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(DataFieldUnits(ifield))) )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'Units')

       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "YearOfRun", &
            HE5T_NATIVE_INT, n_one, 2007 )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'YearOfRun')

       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "MonthOfRun", &
            HE5T_NATIVE_INT, n_one, imonth )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'MonthOfRun')

       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "DayOfRun", &
            HE5T_NATIVE_INT, n_one, 0 )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'DayOfRun')

       n_attr = LEN_TRIM(ADJUSTL(TimeOfDay(ifield)))
       he5stat = HE5_SWwrlattr (                                          &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "TimeOfDay", &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(TimeOfDay(ifield)))   )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'TimeOfDay')

       n_attr = LEN_TRIM(ADJUSTL(MERRAVersion(ifield)))
       he5stat = HE5_SWwrlattr (                                                &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "MERRAVersion", &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(MERRAVersion(ifield)))   )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'MERRAVersion')

       n_attr = LEN_TRIM(ADJUSTL(DateOfProduction(ifield)))
       he5stat = HE5_SWwrlattr (                                                &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "DateOfProduction", &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(DateOfProduction(ifield)))   )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'DateOfProduction')

       n_attr = LEN_TRIM(ADJUSTL(DataProvider(ifield)))
       he5stat = HE5_SWwrlattr (                                             &
            swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "DataProvider", &
            HE5T_NATIVE_CHAR, n_attr, TRIM(ADJUSTL(DataProvider(ifield))) )
       IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'DataProvider')
  
    END IF

    he5stat = HE5_SWwrlattr (                                             &
         swath_id, TRIM(ADJUSTL(DataFieldNames(ifield))), "GEOSVersion", &
         HE5T_NATIVE_INT, n_one, 5 )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrlattr', 'GEOSVersion')
  

    RETURN
  END SUBROUTINE write_field_attributes


  SUBROUTINE write_swath_attributes ( swath_id )
  
    IMPLICIT NONE
  
    ! ---------------
    ! Input variables
    ! ---------------
    INTEGER, INTENT (IN) :: swath_id

    ! --------------
    ! Local variable
    ! --------------
    INTEGER(C_LONG) :: n_attr
    INTEGER         :: he5stat
  
  
    n_attr = LEN_TRIM(ADJUSTL("Gonzalo Gonzalez Abad"))
    he5stat = HE5_SWwrattr (                                             &
         swath_id, "Author", &
         HE5T_NATIVE_CHAR, n_attr, "gonzalo gonzalez abad" )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrgattr', 'Author')

    n_attr = LEN_TRIM(ADJUSTL("Harvard-Smithsonian Center for Astrophysics"))
    he5stat = HE5_SWwrattr (                                             &
         swath_id, "AuthorAffiliation", &
         HE5T_NATIVE_CHAR, n_attr, "Harvard-Smithsonian Center for Astrophysics" )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrgattr', 'AuthorAffiliation')
  
    n_attr = LEN_TRIM(ADJUSTL("ggonzale@cfa.harvard.edu"))
    he5stat = HE5_SWwrattr (                                             &
         swath_id, "AuthorContact", &
         HE5T_NATIVE_CHAR, n_attr, "ggonzale@cfa.harvard.edu" )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrgattr', 'AuthorContact')
  
    n_attr = LEN_TRIM(ADJUSTL("2016-05-17"))
    he5stat = HE5_SWwrattr (                                             &
         swath_id, "DateOfProduction", &
         HE5T_NATIVE_CHAR, n_attr, "2016-05-17" )
    IF ( he5stat /= 0 ) CALL he5_error_stop ('HE5_SWwrgattr', 'DateOfProduction')

    RETURN
  END SUBROUTINE write_swath_attributes


  SUBROUTINE he5_error_stop ( he5routine, fieldname )

    IMPLICIT NONE

    CHARACTER (LEN=*), INTENT (IN) :: he5routine, fieldname

    WRITE (*,'(A)') 'ERROR: '//TRIM(ADJUSTL(he5routine))//' failed for '//TRIM(ADJUSTL(fieldname))
    STOP 1

    RETURN
  END SUBROUTINE he5_error_stop

END MODULE GEOSChem_Helen_he5_module
