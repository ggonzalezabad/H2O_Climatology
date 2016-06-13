close, /all

; Available data (year & month)
years = ['2007']
month = ['01','02','03','04','05','06','07','08','09','10','11','12']

RESTORE, 'MERRA2_H2O.sav'
; Loop for years
FOR iyear = 0, N_ELEMENTS(years)-1 DO BEGIN

    ; Loop for months
    FOR imonth = 0, 11 DO BEGIN
        
        ; Generate filename for the output ASCII
        file = './Helen_MERRA2/ctm.'+years[iyear]+month[imonth]+'.dat'
        print, 'Writing to '+file

        OPENW, LUN, file, /GET_LUN

        ; Loop for longitudes
        FOR ilon = 0, nlon-1 DO BEGIN

            ; Loop for latitudes
            FOR ilat = 0, nlat-1 DO BEGIN

                ; Print longitudes, latitudes and surface pressure to file
                PRINTF, LUN, lon[ilon], lat[ilat], psurface[imonth,ilon,ilat]
                ; Print H2O to file
                PRINTF, LUN, h2ovmr[imonth,ilon,ilat,0:46],      FORMAT = '(47(1x,E13.5))'
                ; Print H2O standard deviation to file
                PRINTF, LUN, h2ovmr_sd[imonth,ilon,ilat,0:46],   FORMAT = '(47(1x,E13.5))'
                ; Print Temperature to file
                PRINTF, LUN, TEMPERATURE[imonth,ilon,ilat,0:46], FORMAT = '(47(1x,E13.5))'

            ENDFOR ; End latitudes loop

        ENDFOR ; End longitudes loop

        CLOSE, /ALL

    ENDFOR ; End months loop

ENDFOR ; End years loop

CLOSE, /all
endall:

END
