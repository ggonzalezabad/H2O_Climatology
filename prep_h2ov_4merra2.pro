pro prep_h2ov_4merra2,fin,lat,lon,ps,pcenter,h2onum_gc,roul,droulroul,temperature_gc, $
                      fout=fout,saveit=saveit

if not keyword_set(saveit) then saveit = 'y'
if not keyword_set(fout) then fout = 'h2ov_clim.txt'

miuh2o = 0.018
miudry = 0.029
missing = 1.e+15
;------------------------------------------------------------------
modelinfo = ctm_type('GEOS5_47L',res=2)

;------------------------------------------------------------------
ncdf_get_1field,fin,'lon',lons
ncdf_get_1field,fin,'lat',lats
ncdf_get_1field,fin,'lev',levs

trash = min(abs(lons - lon),ii)
trash = min(abs(lats - lat),jj)
print,'use profile at lat/lon = ',lats[jj],lons[ii]
;-------------------------------------------------------------------
ncdf_get_1field,fin,'PS',temp ; Pa
surfp  = temp / 100. ; mbar
ps = surfp [ii,jj]
psfc_pgrid,ps,pcenter,pedge,dp,modelinfo=modelinfo
undefine,temp

;-------------------------------------------------------------------
ncdf_get_1field,fin,'T',temp; K
temperature = reform(temp[ii,jj,*])
rounum = levs * 100. / temperature / 8.314 * 6.022e+23 ; #_wetair / m^3
undefine,temp
;-------------------------------------------------------------------
ncdf_get_1field,fin,'H',temp ; geopotential height [m]
gph = reform(temp[ii,jj,*])

;-------------------------------------------------------------------
ncdf_get_1field,fin,'QV',qv ; kg / kg
temp = (miuh2o/miudry) * (1.d0 - qv) / qv 
nmxr =  1.d0 / (1.d0 + temp)
ind = where(qv gt 1.e+14,cc)
if (cc gt 0) then nmxr[ind] = 0.d0
undefine,temp

roumxr = reform(nmxr[ii,jj,*])
h2onum = roumxr * rounum ; #_h2ov / m^3

;------------------------------------------------------------------
ncdf_get_1field,fin,'Var_QV',temp
trash = sqrt(temp)
output = trash / qv
if (cc gt 0) then output[ind] = 0.

deltamxr = reform(output[ii,jj,*])
;-------------------------------------------------------------------
h2onum_gc = interpol(h2onum,levs,pcenter)
temperature_gc = interpol(temperature,levs,pcenter)
nz = n_elements(pcenter)


gph_gc = interpol(gph,levs,pedge)
dh = shift(gph_gc,-1) - gph_gc
dh = dh[0:nz-1] ; meters

droulroul = interpol(deltamxr,levs,pcenter)
;roumxr_gc = interpol(roumxr,levs,pcenter)
;===================================================================

roul = h2onum_gc * dh * 1.e-4 ; # / cm^2
roul = roul > 0.d0

if (saveit eq 'y') then begin
close,11
openw,11,fout
for i = 0, nz-1 do begin
   k = nz - 1 - i
   printf,11,pcenter[k],roul[k],droulroul[k],format='(E15.2,1x,E21.7,E15.2)'
endfor
close,11
endif

;-------------------------------------------------------------------
return
end
