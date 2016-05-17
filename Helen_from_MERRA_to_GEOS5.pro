; Create file name strings:
fss = '/data/tempo2/hwang/MERRA2/MERRA2_300.instM_3d_ana_Np.2007'
ffs = '.nc4'
fcs =['01','02','03','04','05','06','07','08','09','10','11','12']

; Create lon & lat grid
nlat = 361 & lat = FINDGEN(nlat)*1.0/2.0-90.00
nlon = 576 & lon = FINDGEN(nlon)*5.0/8.0-180.00

; Missing value
miss_val = -1.0e-30

; Molecular weights
miuh2o = 0.018 & miudry = 0.029

; Create variables to hold data
psurface = FLTARR(12,nlon,nlat)
temperature = FLTARR(12,nlon,nlat,47)
h2ovmr = FLTARR(12,nlon,nlat,47)

; Loop over each one of the files
FOR ifile = 0, 11 DO BEGIN

   filename = fss+fcs[ifile]+ffs
   print, filename
   ; Read MERRA fields
   ncdf_get_1field,filename,'lon',lons_merra
   ncdf_get_1field,filename,'lat',lats_merra
   ncdf_get_1field,filename,'lev',levs_merra
   ncdf_get_1field,filename,'PS',ps_merra  ; Pa
   ncdf_get_1field,filename,'T',t_merra ; K
   ncdf_get_1field,filename,'H',h_merra ; geopotential height [m]
   ncdf_get_1field,filename,'QV',qv_merra  ; Specific humidity [kg/kg]
   ncdf_get_1field,filename,'Var_QV',var_qv_merra ; Variance of qv [kg/kg kg/kg]

   ; Convert ps_merra to hPa
   ps_merra = ps_merra / 100.0
   ; Values 1e15 are converted to NAN
   dummy = WHERE(t_merra EQ 1E15, count) & IF (count GE 1) THEN t_merra[dummy] = !Values.F_nan
   dummy = WHERE(h_merra EQ 1E15, count) & IF (count GE 1) THEN h_merra[dummy] = !Values.F_nan
   dummy = WHERE(qv_merra EQ 1E15, count) & IF (count GE 1) THEN qv_merra[dummy] = !Values.F_nan
   dummy = WHERE(var_qv_merra EQ 1E15, count) & IF (count GE 1) THEN var_qv_merra[dummy] = !Values.F_nan

   ; Convert qv_merra to vmr_merra
   dummy = (miuh2o/miudry) * (1.0d0 - qv_merra) / qv_merra
   vmr_merra = 1.0d / ( 1.0d + dummy) * 1e9 ; [ppb]

   ; Loop over latitudes
   FOR ilat = 0, nlat-1 DO BEGIN
      ; Loop over longitudes
      FOR ilon = 0, nlon-1 DO BEGIN
         
         ; Find out the closest MERRA pixel to lon[ilon] & lat[ilat]
         dummy = MIN(ABS(lons_merra - lon[ilon]),jj)
         dummy = MIN(ABS(lats_merra - lat[ilat]),ii)

         ; Save surface for pixel
         psurface[ifile,ilon,ilat] = ps_merra[jj,ii]

         ; Use internal GAMAP function to work out pixel pressure grid
         dummy = GetEta('GEOS5', 47, Psurf=psurface[ifile,ilon,ilat], pressure=pcenter)
         dummy = GetEta('GEOS5', 47, Psurf=psurface[ifile,ilon,ilat], /Edge, pressure=pedge)

         ; Interpolate MERRA temperature to GEOS-Chem
         ; Question to answer (how or where is defined
         ; the MERRA vertical grid and what is it's
         ; relationship with surface pressure). Using
         ; pressure levels is not a very good idea. Why
         ; not to use the hybrid vertical grid.
         temperature[ifile,ilon,ilat,*] = INTERPOL(t_merra[jj,ii,*],levs_merra,pcenter)
         ; If they are /Nan values in h2ovmr[ifile,ilon,ilat,*]
         dummy = WHERE(FINITE(temperature[ifile,ilon,ilat,*], /Nan), count)
         IF (count GE 1) THEN BEGIN
            ; Find non NAN water VMR closest to surface
            FOR ilev = 0, 46 DO BEGIN
               IF (FINITE(temperature[ifile,ilon,ilat,ilev], /NAN)) THEN CONTINUE
               lower = temperature[ifile,ilon,ilat,ilev]
               BREAK
            ENDFOR
            temperature[ifile,ilon,ilat,dummy] = lower
         ENDIF
         
         ; Interpolate vmr_merra.         
         h2ovmr[ifile,ilon,ilat,*] = INTERPOL(vmr_merra[jj,ii,*],levs_merra,pcenter)
         ; If they are /Nan values in h2ovmr[ifile,ilon,ilat,*]
         dummy = WHERE(FINITE(h2ovmr[ifile,ilon,ilat,*], /Nan), count)
         IF (count GE 1) THEN BEGIN
            ; Find non NAN water VMR closest to surface
            FOR ilev = 0, 46 DO BEGIN
               IF (FINITE(h2ovmr[ifile,ilon,ilat,ilev], /NAN)) THEN CONTINUE
               lower = h2ovmr[ifile,ilon,ilat,ilev]
               BREAK
            ENDFOR
            h2ovmr[ifile,ilon,ilat,dummy] = lower
         ENDIF
      ENDFOR ; End longitudes loop
   ENDFOR ; End latitudes loop

ENDFOR ; End files (month loop)

dummy = WHERE(FINITE(temperature, /NAN), COUNT) & IF (count GE 1) THEN temperature[dummy] = miss_val
dummy = WHERE(FINITE(h2ovmr, /NAN), COUNT) & IF (count GE 1) THEN h2ovmr[dummy] = miss_val

; Convert lon & lat grid to center
lon = lon + (0.5*5.0/8.0) & lat = lat + ( 0.5*1.0/2.0)

SAVE, lon, lat, nlon, nlat, psurface, temperature, h2ovmr, filename = 'MERRA2_H2O.sav'

endall:
END
