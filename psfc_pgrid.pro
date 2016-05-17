pro psfc_pgrid,psfc,pcenter,pedge,dp,modelinfo=modelinfo

nargs = n_params()
if nargs lt 1 then begin
  print,' psfc_pgrid,psfc,pcenter,pedge,dp,modelinfo=modelinfo'
  print,' calculate pedge and pcenter , dp from psfc'
  return
endif
if not keyword_set(modelinfo) then $
modelinfo=ctm_type('GEOS5_47L',res=2)

gridinfo = ctm_grid(modelinfo)
AP = gridinfo.AP
BP = gridinfo.BP
nn = n_elements(AP)
;---------------------------------------------------------
ss = size(psfc)
ndim = ss(0)

if ndim eq 1 then begin
   pedge = AP + BP * psfc
   pcenter = (shift(pedge,-1)+pedge) / 2.0
   pcenter = pcenter (0:nn-2)
   dp = -(shift(pedge,-1)-pedge)
   dp = dp (0:nn-2)
endif else if ndim eq 2 then begin
   dim1 = ss(1) & dim2 = ss(2)
   pedge = fltarr(dim1,dim2,nn)
   pcenter = fltarr(dim1,dim2,nn-1)
   dp = fltarr(dim1,dim2,nn-1)
   for i = 0, dim1-1 do begin
      for j = 0, dim2-1 do begin
         tmp = AP(*) + BP(*) * psfc(i,j)
         tmpc = (shift(tmp,-1)+tmp)/2.
         detp = -(shift(tmp,-1)-tmp)
         pedge(i,j,*) = tmp(*)
         pcenter(i,j,*) = tmpc(0:nn-2)
         dp(i,j,*) = detp(0:nn-2)
      endfor
   endfor
endif else if (ndim eq 0) then begin
   pedge = AP + BP * psfc
   pcenter = (shift(pedge,-1)+pedge) / 2.0
   pcenter = pcenter[0:nn-2]
   dp = pedge - shift(pedge,-1)
   dp = dp [0: nn -2]
endif else begin
   message,'dimension not supported.'
endelse
;--------------------------------------------------------
return
end   
