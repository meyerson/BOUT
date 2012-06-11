function GInterpol, y, x, xnew, order=order
;
; Global polynomial interpolation in 1D
;-----------------------------------------------------------------;

 if not keyword_set(ORDER) then ORDER=1 ;-linear fit by default

 n=N_ELEMENTS(x)
 if (n ne N_ELEMENTS(y)) then STOP, 'x and y dimensions must be same in GInterpol'
 ORDER=MIN([ORDER,n-1]) ;-otherwise ill-conditioned

 nnew=N_ELEMENTS(xnew)

 coef=POLY_FIT(x,y,order)
 ynew=fltarr(nnew)

 for iord=0,order do begin
     ynew=ynew+coef[iord]*(xnew^iord)
 endfor

return, ynew
end


function int2x, xarr, farr, narr, SORT=SORT, X0=X0, Y0=Y0, ORDER=ORDER, DEBUG=DEBUG
;
; Integrate f(x)dx from x[0] to x[i]
; where x and f are given by arrays xarr and farr 
;
; ORDER is the order of polynomial used to find the reference point
;------------------------------------------------------------------;

 yarr=fltarr(narr)

  for i=1,narr-1 do begin
   xarrNow=xarr[0:i]
   farrNow=farr[0:i]
   yarr[i]=0.5*TOTAL((xarrNow[1:i]-xarrNow[0:i-1])*(farrNow[1:i]+farrNow[0:i-1])) ;-trapezoid rule
  endfor


  ;-interpolate result to given value X0
  if (N_ElEMENTS(X0) gt 0) then begin

      if (order eq -1) then begin 
          ;-local interpolation to 4 nearby points
          Y0=INTERPOL(yarr,xarr,X0,/Q)
      endif else begin
          ;-global interpolation to full set
          Y0=GInterpol(yarr,xarr,X0,order=order)
      endelse

  endif

if keyword_set(DEBUG) then STOP
return, yarr
end



function Gindgen, i0, j0
;-generalized indgen, [i0,...,j0]

 if (j0 gt i0) then return, i0+indgen(j0-i0+1)

 if (j0 lt i0) then return, j0+REVERSE(indgen(i0-j0+1))

end



pro ReflectProfile, prof, ix0
;
; Make a symmetric profile by reflecting values from the outer part to
; the inner part
; Inputs: vector prof
;         index ix0
;
; inner part [0,ix0], outer part [ix0+1:*]
;---------------------------------------------------------------------;

nx=n_elements(prof)

if (ix0 gt 0 AND ix0 lt (nx-1)) then begin

 i_out=Gindgen(ix0+1,nx-1)
 n_out=(nx-1)-(ix0+1)+1

 i_in=Gindgen(0, ix0)
 n_in=ix0+1

 if (n_in gt n_out) then begin

     i_in1=Gindgen(0,n_in-n_out-1)
     i_in2=Gindgen(n_in-n_out,ix0)
     prof[i_in2]=REVERSE(prof[i_out])
     prof[i_in1]=prof[i_in2[0]]

 endif else begin

     prof[i_in]=REVERSE(prof[i_out[0:n_in-1]])

 endelse

endif

end



function sideval, f, ipue, jrue 
;
; calculate values on UEDGE grid cell sides from
; known values on the four corners
;
; Inputs: ipue-UEDGE poloidal index
;         jrue-UEDGE radial index
;-------------------------------------------------;

 fSW=f[ipue,jrue,1] 
 fSE=f[ipue,jrue,2] 
 fNW=f[ipue,jrue,3] 
 fNE=f[ipue,jrue,4] 

  fEast=0.5*(fNE+fSE)
   fWest=0.5*(fNW+fSW)
    fNorth=0.5*(fNE+fNW)
     fSouth=0.5*(fSE+fSW)
      fCntr=f[ipue,jrue,0]

return, {east:fEast, west:fWest, north:fNorth, south:fSouth, cntr:fcntr}
end



function remlim_arr, arr, sz, ix_lim, nx, ixdim, info=info, debug=debug
;
; Remove limiter cells from first dimension of an array
;------------------------------------------------------------;

if not keyword_set(INFO) then info=0

CASE sz[0] of

    1:  begin
        if (info) then print, '...Reformatting 1D array'
        indx=[gindgen(0,ix_lim-1), gindgen(ix_lim+2,nx-1)]
        arrnew=arr[indx, *]
        return, arrnew
    end


    2:  begin
        if (info) then print, '...Reformatting 2D array'

        indx=[gindgen(0,ix_lim-1), gindgen(ix_lim+2,nx-1)]
        arrnew=arr[indx,*]
        ;stop
        return, arrnew 
    end


    3:  begin
        if (info) then print, '...Reformatting 3D array'

        indx=[gindgen(0,ix_lim-1), gindgen(ix_lim+2,nx-1)]
        arrnew=arr[indx,*,*]
        ;stop
        return, arrnew
    end

    else: STOP, '...unexpected dimensionality !!!'


ENDCASE

end




function remlim_str, d, ix_lim, nx, info=info
;
; Remove guard cells around limiter for all arrays
; members of the input structure
;------------------------------------------------------------;

 print, 'Removing limiter cells from structure'

 dnew=CREATE_STRUCT('blank','blank') ;-initialize new structure
 dnames=TAG_NAMES(d)
 nn=n_elements(dnames)
  
 for i=0,nn-1 do begin

     if keyword_set(info) then print, 'testing ', dnames[i], size(d.(i))

     sz=(size(d.(i)))         ;-field dimension

     if (sz[0] gt 0) then begin ;-all arrays

         arr=d.(i)

                                ;-find what dimension is x
         ixdim=where(sz eq nx) 

         if (ixdim[0] gt -1) then begin
                                ;-found such dimension in array

             arrNew=RemLim_arr(arr, sz, ix_lim, nx, ixdim[0], info=info)

         endif else begin
                                ;-array but has no such dimension
             arrNew=arr
         endelse

     endif else begin

                                ;-not an array
         arrNew=d.(i)
     
     endelse

                                ;-add new field
     dnew=CREATE_STRUCT(dnew, dnames[i], arrNew)

        
 endfor

;
;
;
return, dnew
end



pro showDomain, d1, ixmin, ixmax, jymin, jymax, color=color, FILL=FILL, var=var, $
EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR
;
; Show variable on grid or just show the domain itself
;------------------------------------------------------------------------------;

   if not keyword_set(EPS) then EPS=1e-3

   if keyword_set(var) then begin
       if not keyword_set(MAXVAR) then maxVar=MAX(var)
       if not keyword_set(MINVAR) then minVar=MIN(var)
   endif


   for ix=ixMin,ixMax do begin
    for jy=jyMin,jyMax do begin

	rr=REFORM(d1.RM_COM[jy,ix,[1,2,4,3,1]])
	zz=REFORM(d1.ZM_COM[jy,ix,[1,2,4,3,1]])

            if keyword_set(var) then begin
             color_now=3+250*(var[jy,ix]-minVar)/MAX([(maxVar-minVar),1e-10])

             ;-avoid over- and under-shots
             if (color_now lt 3 OR color_now gt 253) then color_now=255

             ;-mark zero level by black
             if (ABS(var[jy,ix]) lt EPS) then color_now=0 ;-this is for debugging
            endif else begin
             color_now=color
            endelse


            if keyword_set(fill) then begin
             POLYFILL, rr, zz, color=color_now
            endif else begin 
             oplot, rr, zz, color=color_now
            endelse

       plots, d1.RM_COM[jy,ix,[0]], d1.ZM_COM[jy,ix,[0]], psym=3

    endfor
   endfor
;
end



pro showGrid, d1, d2, d3, FILL=FILL, XRANGE=XRANGE, YRANGE=YRANGE, var=var, title=title,$
EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR
;
; Plot a given variable on the grid domain, or just the domain itself
;-------------------------------------------------------------------------;

  if keyword_set(VAR) then begin
      loadct,39
      if not keyword_set(MAXVAR) then maxVar=MAX(var)
      if not keyword_set(MINVAR) then minVar=MIN(var)
      print, "MIN=", minVar, ", MAX=", maxVar
  endif else begin
      tek_color
  endelse


  subtitle=d1.RUNIDG_GRD

  if keyword_set(XRANGE) then XSTYLE=1
  if keyword_set(YRANGE) then YSTYLE=1

  plot, d1.RM_COM[*,*,0], d1.ZM_COM[*,*,0], /iso, psym=3, subtit=subtitle,/NOD, $
    XRANGE=XRANGE, YRANGE=YRANGE, XSTYLE=XSTYLE, YSTYLE=YSTYLE, title=title


   IF (d1.nxpt_com lt 2) then begin

    ;-core
    ShowDomain, d1, 0, d1.IYSPTRX1_COM[0], d1.IXPT1_COM[0]+1, d1.IXPT2_COM[0],$
	 col=2, FILL=FILL, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-sol
    ShowDomain, d1, d1.IYSPTRX1_COM[0]+1, d1.NY_COM+1, d1.IXLB_COM[0], d1.IXRB_COM[0]+1, $
	col=3, FILL=FILL, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-PF, outer
    ShowDomain, d1, 0, d1.IYSPTRX1_COM[0], d1.IXPT2_COM[0]+1, d1.IXRB_COM[0]+1, $
	col=4, FILL=FILL, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-PF, inner
    ShowDomain, d1, 0, d1.IYSPTRX1_COM[0], d1.IXLB_COM[0], d1.IXPT1_COM[0],$ 
	col=5, FILL=FILL, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

  ENDIF ELSE BEGIN

    ;-core inboard side
    ShowDomain, d1, 0, d1.IYSPTRX1_COM[0], d1.IXPT1_COM[0]+1, d1.IXPT2_COM[0], $
                col=4,/f, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-core outboard side
    ShowDomain, d1, 0, d1.IYSPTRX1_COM[0], d1.IXPT1_COM[1]+1, d1.IXPT2_COM[1], $
                col=2,/f, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-internal sol inboard side
    ShowDomain, d1, d1.IYSPTRX1_COM[0]+1,d1.IYSPTRX2_COM[0], 0, d1.IXPT2_COM[0], $
                col=5,/f, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-internal sol outboard side
    ShowDomain, d1, d1.IYSPTRX1_COM[0]+1,d1.IYSPTRX2_COM[0], d1.IXPT1_COM[1]+1, $
                d1.IXRB_COM[1]+1, col=3,/f, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-external sol inboard side
    ShowDomain, d1, d1.IYSPTRX2_COM[0]+1,d1.ny_com+1, d1.IXLB_COM[0], d1.IXRB_COM[0]+1, $
                col=7,/f, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-external sol outboard side
    ShowDomain, d1, d1.IYSPTRX2_COM[0]+1,d1.ny_com+1, d1.IXLB_COM[1], d1.IXRB_COM[1]+1, $
                col=8,/f, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-internal pf inboard side
    ShowDomain, d1, 0, d1.IYSPTRX1_COM[0], d1.IXLB_COM[0], d1.IXPT1_COM[0], $
                col=6,/f, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-internal pf outboard side
    ShowDomain, d1, 0, d1.IYSPTRX1_COM[0], d1.IXPT2_COM[1]+1, d1.IXRB_com[1]+1, $
                col=9,/f, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-external pf inboard side
    ShowDomain, d1, 0, d1.IYSPTRX2_COM[0], d1.IXPT2_COM[0]+1, d1.IXRB_COM[0]+1, $
                col=11,/f, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

    ;-external pf outboard side
    ShowDomain, d1, 0, d1.IYSPTRX2_COM[0], d1.IXLB_COM[1], d1.IXPT1_COM[1],  $
                col=12,/f, var=var, EPS=EPS, MINVAR=MINVAR, MAXVAR=MAXVAR

  ENDELSE 

  ;;CONTOUR, var, d1.RM_COM[*,*,0], d1.ZM_COM[*,*,0], lev=[-1,0,1],/over

end



pro ImportUE, d1, d2, d3, PATH=PATH, NOREFORMAT=NOREFORMAT,no
;
; Import data from UEDGE data files
;-------------------------------------------------;

  if not keyword_set(PATH) then PATH='.'

  f1=PATH+'/gridue.pdb' 
  f2=PATH+'/uedgegrd.pdb' 
  f3=PATH+'/uedgeout.pdb' 

       print, ''
       print, 'Reading data from files:'
         print, f1
         print, f2
         print, f3

       ;-read UEDGE variables from PDB files
       d1=pd_import(f1)
       d2=pd_import(f2)
       d3=pd_import(f3)


       if not keyword_set(NOREFORMAT) then begin
           d1=Reformat_struc(d1)
           d2=Reformat_struc(d2)
           d3=Reformat_struc(d3)
       endif


       if (d1.IX_LIM_COM gt 0) then begin
                                ;remove limiter cells from all structures
           d1=REMLIM_STR(d1,d1.IX_LIM_COM,d1.nx_com+2,info=0)
           d2=REMLIM_STR(d2,d1.IX_LIM_COM,d1.nx_com+2,info=0)
           d3=REMLIM_STR(d3,d1.IX_LIM_COM,d1.nx_com+2,info=0)
           d1.nx_com=d1.nx_com-2
           d1.IXPT2_COM[0]=d1.IXPT2_COM[0]-2
           d1.IXRB_COM[0]=d1.IXRB_COM[0]-2
       endif

end



pro get_theta, d1, d2, d3, theta, dTheta, theta_all
;
; Define Theta coordinate on the grid
;--------------------------------------------------------------------;

;;-Theta uniform in poloidal index, jyseps1+1 <=> Theta=0., jyseps2 <=> Theta=2*PI,
;;-from x-point to x-point the full poloidal angle is 2*PI

print, ''
print, 'Calculating poloidal angle theta'


        Theta=fltarr(d1.nx_com+2, d1.ny_com+2) 


        IF (d1.nxpt_com gt 1) THEN BEGIN ;-double-null

          dTheta = 2*!PI/float((d1.IXPT2_COM[0] - d1.IXPT1_COM[0])+$
                               (d1.IXPT2_COM[1] - d1.IXPT1_COM[1])+1)
          print, "dTheta=", dTheta



          for ip=0,d1.nx_com+1 do begin
           for jr=0,d1.ny_com+1 do begin


             CASE 1 OF

                 ;-core outboard side
                 ((jr ge 0) AND (jr le d1.IYSPTRX1_COM[0]) AND $ 
                  (ip ge (d1.IXPT1_COM[1]+1)) AND (ip le d1.IXPT2_COM[1])) : $
                 BEGIN
                       Theta[ip,jr]=2*!PI+float(ip-((d1.IXPT2_COM[1]+1))-1)*dTheta	
                 END


                 ;-core inboard side
	         ((jr ge 0) AND (jr le d1.IYSPTRX1_COM[0]) AND $
                  (ip ge (d1.IXPT1_COM[0]+1)) AND (ip le d1.IXPT2_COM[0])) : $
                 BEGIN
                       Theta[ip,jr]=float(ip-((d1.IXPT1_COM[0]+1)))*dTheta	
                 END

            
                 ;-internal sol outboard side
                 ((jr ge (d1.IYSPTRX1_COM[0]+1)) AND (jr le d1.IYSPTRX2_COM[0]) AND $ 
                   (ip ge (d1.IXPT1_COM[1]+1)) AND (ip le (d1.IXRB_COM[1]+1))) : $
                 BEGIN
                       Theta[ip,jr]=2*!PI+float(ip-(d1.IXPT2_COM[1]+1)-1)*dTheta	
                 END


                 ;-internal sol inboard side
                 ((jr ge (d1.IYSPTRX1_COM[0]+1)) AND (jr le d1.IYSPTRX2_COM[0]) AND $ 
                  (ip ge 0) AND (ip le d1.IXPT2_COM[0])) : $
                 BEGIN
                       Theta[ip,jr]=float(ip-((d1.IXPT1_COM[0]+1)))*dTheta	
                 END


                 ;-external sol outboard side
                 ((jr ge (d1.IYSPTRX2_COM[0]+1)) AND (jr le (d1.ny_com+1)) AND $ 
                  (ip ge d1.IXLB_COM[1]) AND (ip le (d1.IXRB_COM[1]+1))) : $
                 BEGIN
                       Theta[ip,jr]=2*!PI+float(ip-((d1.IXPT2_COM[1]+1))-1)*dTheta	
                 END


                 ;-external sol inboard side
                 ((jr ge (d1.IYSPTRX2_COM[0]+1)) AND (jr le (d1.ny_com+1)) AND $ 
                  (ip ge d1.IXLB_COM[0]) AND (ip le (d1.IXRB_COM[0]+1))) : $
                 BEGIN
                       Theta[ip,jr]=float(ip-((d1.IXPT1_COM[0]+1)))*dTheta	
                 END

                 ;-internal pf outboard side
                 ((jr ge 0) AND (jr le d1.IYSPTRX1_COM[0]) AND $
                  (ip ge (d1.IXPT2_COM[1]+1)) AND (ip le (d1.IXRB_com[1]+1))) : $
                 BEGIN
                       Theta[ip,jr]=2*!PI+float(ip-((d1.IXPT2_COM[1]+1))-1)*dTheta	
                 END


                 ;-internal pf inboard side
                 ((jr ge 0) AND (jr le d1.IYSPTRX1_COM[0]) AND $ 
                  (ip ge d1.IXLB_COM[0]) AND (ip le d1.IXPT1_COM[0])) : $
                 BEGIN
                       Theta[ip,jr]=float(ip-((d1.IXPT1_COM[0]+1)))*dTheta	
                 END


                 ;-external pf outboard side
                 ((jr ge 0) AND (jr le d1.IYSPTRX2_COM[0]) AND $ 
                  (ip ge d1.IXLB_COM[1]) AND (ip le d1.IXPT1_COM[1])) : $
                 BEGIN
                       Theta[ip,jr]=2*!PI+float(ip-((d1.IXPT2_COM[1]+1))-1)*dTheta	
                 END


                 ;-external pf inboard side
                 ((jr ge 0) AND (jr le d1.IYSPTRX2_COM[0]) AND $
                  (ip ge (d1.IXPT2_COM[0]+1)) AND (ip le (d1.IXRB_COM[0]+1))) : $
                 BEGIN
                       Theta[ip,jr]=float(ip-((d1.IXPT1_COM[0]+1)))*dTheta	
                 END


                 ;;ELSE: 

             ENDCASE

           endfor
          endfor

        ENDIF ELSE BEGIN ;single null geometry

          dTheta = 2*!PI/float(d1.IXPT2_COM[0] - d1.IXPT1_COM[0])
          print, "dTheta=", dTheta

          for ip=0,d1.nx_com+1 do begin
           for jr=0,d1.ny_com+1 do begin
	     Theta[ip,jr]=float(ip-((d1.IXPT1_COM[0]+1)))*dTheta	
           endfor
          endfor

        ENDELSE


;-now define Theta at each grid corner
        Theta_all=fltarr(d1.nx_com+2, d1.ny_com+2, 5)

        FOR i=0,d1.nx_com+1 DO BEGIN
            FOR j=0,d1.ny_com+1 DO BEGIN
                ;
                Theta0=Theta[i,j]
                Theta_all[i,j,0]=Theta0 ;-center
                Theta_all[i,j,1]=Theta0-0.5*dTheta ;-SW
                Theta_all[i,j,3]=Theta0-0.5*dTheta ;-NW
                Theta_all[i,j,2]=Theta0+0.5*dTheta ;-SE
                Theta_all[i,j,4]=Theta0+0.5*dTheta ;-NE
                ;
            ENDFOR
        ENDFOR

end


pro Get_Jacobian, DEBUG=DEBUG,$
                  d1, d2, d3, dTheta, Rthe, Zthe, Rpsi, Zpsi, $
                  dpsi, dlthe, hthe, qsafe, kapagc1, kapagc2, kapanc
;
; Get components of the Jacobian matrix and related quantities. The
; poloidal grid cell size and partial derivatives w/respect to Theta and Psi
; are calculated using the four corners of the grid
;---------------------------------------------------------------------;

print, ''
print, 'Calculating elements of the Jacobian matrix'
if keyword_set(DEBUG) then STOP

       dlthe=fltarr(d1.nx_com+2, d1.ny_com+2) 
       dpsi=fltarr(d1.nx_com+2, d1.ny_com+2) 

       rthe=fltarr(d1.nx_com+2, d1.ny_com+2) 
       zthe=fltarr(d1.nx_com+2, d1.ny_com+2) 

       rpsi=fltarr(d1.nx_com+2, d1.ny_com+2) 
       zpsi=fltarr(d1.nx_com+2, d1.ny_com+2) 

       kapagc1=fltarr(d1.nx_com+2, d1.ny_com+2) 
       kapanc=fltarr(d1.nx_com+2, d1.ny_com+2) 


          FOR ip=0,d1.nx_com+1 DO BEGIN
           FOR jr=0,d1.ny_com+1 DO BEGIN

               rr=Sideval(d1.rm_com,ip,jr)
               zz=Sideval(d1.zm_com,ip,jr)
               bb=Sideval(d1.b_com,ip,jr)
               pp=Sideval(d1.psi_com,ip,jr)

               dlthe[ip,jr]=MAX([sqrt((rr.East-rr.West)^2+(zz.East-zz.West)^2),1e-10])
               dpsi[ip,jr]=pp.North-pp.South
               if (abs(dpsi[ip,jr]) lt 1e-12) then dpsi[ip,jr]=1e-12

               ;-poloidal derivs   
               Rthe[ip,jr]=(rr.East-rr.West)/dTheta
               Zthe[ip,jr]=(zz.East-zz.West)/dTheta
               kapagc1[ip,jr]=(1./(8.*!PI))*(bb.East^2-bb.West^2)/dTheta

               ;-radial derivs   
               Rpsi[ip,jr]=(rr.North-rr.South)/dpsi[ip,jr]
               Zpsi[ip,jr]=(zz.North-zz.South)/dpsi[ip,jr]
               kapanc[ip,jr]=(1./(8.*!PI))*(bb.North^2-bb.South^2)/dpsi[ip,jr]

           ENDFOR
         ENDFOR

if keyword_set(DEBUG) then STOP

        ;-metric coefficient
        hthe=dlthe/dTheta

        ;-define locally, for convenience
        Btxy=d1.bphi_com[*,*,0]
        Bpxy=d1.bpol_com[*,*,0]
        Bxy=d1.b_com[*,*,0]
        Rxy=d1.rm_com[*,*,0]

        ;-local pitch factor
        qsafe=hthe*Btxy/(Rxy*Bpxy)

        ;-put proper normalization coefficients
        kapagc1 = kapagc1*(4.*!pi*Btxy/(hthe*Bxy^3.))
        kapagc2 = kapagc1*(Bpxy/(Btxy*Rxy*hthe))*(Zthe*Zpsi + Rthe*Rpsi)
        kapanc = (4.*!PI/Bxy)*kapanc

;
;
if keyword_set(DEBUG) then STOP
end



pro Get_Refs, d1, d2, d3, hthe, gjy0, gjy1, R0, hthe0, Ni_x, Te_x, Ti_x, Vi_x, $
              bmag, iNixnorm, jNixnorm
;
; Find normalization values for geometry and plasma params
;---------------------------------------------------------------------;

print, ''
print, 'Calculating normalization constants'

       ;-select reference values R0 and hthe0, ir=1 remains in BOUT domain
       R0=max(d1.rm_com[1:d1.NX_COM,1,0], gjy0) ;-skip guard cells where values can be crazy
       dummy=min(d1.rm_com[1:d1.NX_COM,1,0], gjy1)

       gjy0=gjy0+1 ;-compensate for skipping the guard cells
       gjy1=gjy1+1 

       hthe0=hthe[gjy0,1]
       bmag=d1.b_com[gjy0,1,0]


       ;print, "UEDGE outer midplane index ", gjy0 ;;-adjust for guard cells
       print, "R0=", R0
       print, "hthe0=", hthe0
       print, "bmag=", bmag
       print, ""


           ;-select reference values Ni_x, Te_x etc.
           
           ;-old way of choosing Ni_x as a maximum value of Ni (may be in divertor leg)
           ;-new way of choosing Ni_x at outer midplane
           iNixnorm = 0 ;-radial
           jNixnorm = gjy0-1 ;-poloidal

            Ni_x = d3.NI___1__VALUE[gjy0,1]/1e20
            Vi_x = d3.UP___1__VALUE[gjy0,1]
            Te_x = d3.TE____EV_VALUE[gjy0,1]
            Ti_x = d3.TI____EV_VALUE[gjy0,1]
            phi_x = d3.PHI____VALUE[gjy0,1]
    

        print, "Ni_x=", Ni_x, " x 1e20 1/m^3"
        print, "Vi_x=", Vi_x, " m/s"
        print, "Te_x=", Te_x, " eV"
        print, "Ti_x=", Ti_x, " eV"
        print, "phi_x=", phi_x, " V"

end



pro Get_Psi_Derivs, d1, d2, d3, qsafe, dqdpsi, kappaN, kappaV, kappaTe, kappaTi, VE0
;
; Calculate derivatives d/dPsi for quantities that are defined on cell
; centers
;---------------------------------------------------------------------;

print, ''
print, 'Calculating radial derivatives'

    dqdpsi = fltarr(d1.nx_com+2, d1.ny_com+2)
    kappaN = fltarr(d1.nx_com+2, d1.ny_com+2)
    kappaV = fltarr(d1.nx_com+2, d1.ny_com+2)
    kappaTe = fltarr(d1.nx_com+2, d1.ny_com+2)
    kappaTi = fltarr(d1.nx_com+2, d1.ny_com+2)
    VE0 = fltarr(d1.nx_com+2, d1.ny_com+2)


    ;-since dpsi was taken by absolute value we need to do same with psixy
    if (min(DERIV(d1.psi_com[5,*,0])) lt 0.) then begin
        sign=-1 
        print, 'Poloidal field negative'
    endif else begin
        sign=1
        print, 'Poloidal field positive'
    endelse



    FOR ip=0,d1.nx_com+1 DO BEGIN

        kappaN[ip,*]  = -DERIV(d1.psi_com[ip,*,0], d3.NI___1__VALUE[ip,*])/1e20 
        kappaV[ip,*]  = -DERIV(d1.psi_com[ip,*,0], d3.UP___1__VALUE[ip,*])
        kappaTi[ip,*] = -DERIV(d1.psi_com[ip,*,0], d3.TI____EV_VALUE[ip,*])
        kappaTe[ip,*] = -DERIV(d1.psi_com[ip,*,0], d3.TE____EV_VALUE[ip,*])
        VE0[ip,*]     = -DERIV(d1.psi_com[ip,*,0], d3.PHI____VALUE[ip,*])
        dqdpsi[ip,*]  =  DERIV(d1.psi_com[ip,*,0], qsafe[ip,*])

    ENDFOR

;
;
;
end



function Int_subd, x, f, x0, ipMin, ipMax, irMin, irMax, order=order,$
                   offset=offset, DEBUG=DEBUG
;
; Calculate Int f dX _x0^x over a logically 
; rectangular subdomain [ipmin, ipmax, irmin, irmax]
;------------------------------------------------------------------------;

 if not keyword_set(OFFSET) then OFFSET=fltarr(irMax-irMin+1)

 res=fltarr(ipMax-ipMin+1,irMax-irMin+1)

 for ir=irMin, irMax do begin
            
     ip_line=Gindgen(ipMin, ipMax) ;-set of ip indices
     n_line=N_ELEMENTS(ip_line)
     x_line=x[ip_line,ir]
     
     f_line=f[ip_line,ir]
     res_line=Int2x(x_line, f_line, n_line, x0=x0, y0=y0, order=order, DEBUG=DEBUG)
     res[*,ir-irMin]=res_line-y0+offset[ir-irMin]
     
 endfor

if keyword_set(DEBUG) then STOP
return, res
end



pro Int_dTheta, d1, d2, d3, gjy0, gjy1, theta, var_in, var_out
;
; Calculate integral Int f dTheta over full domain
;
; gjy0 is the reference outboard location 
; gjy1 is the reference inboard location
; var_in represents integrand function f
; var_out is the integral
;-------------------------------------------------------------------------;

print, ''
print, 'Calculating poloidal integral Int{ f dTheta}'

 var_out = fltarr(d1.nx_com+2, d1.ny_com+2)


 IF (d1.nxpt_com gt 1) THEN BEGIN ;-double-null case


;-Begin core domain-


       ;-core outboard
       ipMinCoreOut=d1.IXPT1_COM[1]+1
       ipMaxCoreOut=d1.IXPT2_COM[1]
       irMin=0
       irMax=d1.IYSPTRX1_COM[0]
       tOffset=theta[gjy0,0]
       ;
       var_out[ipMinCoreOut:ipMaxCoreOut, irMin:irMax]=$
         Int_subd(theta, var_in, tOffset, $
                  ipMinCoreOut, ipMaxCoreOut, irMin, irMax, order=-1)  


       ;-core inboard, use offset from outboard part
       ipMinCoreIn=d1.IXPT1_COM[0]+1
       ipMaxCoreIn=d1.IXPT2_COM[0]
       irMin=0
       irMax=d1.IYSPTRX1_COM[0]
       tOffset=theta[ipMinCoreOut,0]
       yOffset=var_out[ipMinCoreOut, irMin:irMax]
       ;
       var_out[ipMinCoreIn:ipMaxCoreIn, irMin:irMax]=$
         Int_subd(theta, var_in, tOffset, $
                  ipMinCoreIn, ipMaxCoreIn, irMin, irMax, offset=yOffset, order=-1)
       

;-Finished with the core domain-


;-Begin internal sol domain-

       ;-internal sol outboard
       ipMinIsolOut=d1.IXPT1_COM[1]+1
       ipMaxIsolOut=d1.IXRB_COM[1]+1
       irMin=d1.IYSPTRX1_COM[0]+1
       irMax=d1.IYSPTRX2_COM[0]
       tOffset=theta[gjy0,0]
       ;
       var_out[ipMinIsolOut:ipMaxIsolOut, irMin:irMax]=$
         Int_subd(theta, var_in, tOffset, $
                  ipMinIsolOut,ipMaxIsolOut, irMin, irMax, order=-1)


       ;-internal sol inboard, use offset from outboard part
       ipMinIn=d1.IXLB_COM[0]
       ipMaxIn=d1.IXPT2_COM[0]
       irMin=d1.IYSPTRX1_COM[0]+1
       irMax=d1.IYSPTRX2_COM[0]
       tOffset=theta[ipMinIsolOut,0]
       yOffset=var_out[ipMinIsolOut, irMin:irMax]
       ;
       var_out[ipMinIn:ipMaxIn, irMin:irMax]=$
         Int_subd(theta, var_in, tOffset, $
                  ipMinIn,ipMaxIn, irMin, irMax, order=-1, offset=yOffset)

;-Finished internal sol-


;-Begin external sol outboard side-

       ;-external sol outboard side
       ipMin = d1.IXLB_COM[1]
       ipMax = d1.IXRB_COM[1]+1
       irMin = d1.IYSPTRX2_COM[0]+1
       irMax = d1.ny_com+1
       tOffset=theta[gjy0,0]
       ;
       var_out[ipMin:ipMax, irMin:irMax]=$
         Int_subd(theta, var_in, tOffset, $
                  ipMin, ipMax, irMin, irMax, order=-1) 

;-Finished external sol outboard side-


;-Begin inboard external SOL-

       ipMin=d1.IXLB_COM[0]
       ipMax=d1.IXRB_COM[0]+1
       irMin=d1.IYSPTRX2_COM[0]+1
       irMax=d1.NY_COM+1

       ;-for offset values extrapolate radially from internal sol
       tOffset=theta[gjy1,0] ;-inboard midplane
       yset1=var_out[gjy1, 0:irMin-1]
       pset1=d1.psi_com[gjy1,0:irMin-1,0]
       pset2=d1.psi_com[gjy1,irMin:irMax,0]
       yset2=INTERPOL(yset1,pset1,pset2,/q)
       yOffset=yset2
       ;
       var_out[ipMin:ipMax, irMin:irMax]=$
         Int_subd(theta, var_in, tOffset, $
                  ipMin, ipMax, irMin, irMax, order=1, offset=yOffset)

;-Finish inboard external SOL-


;-Begin PF regions-

      ;-extrapolate from external inboard SOL to upper inboard PF
       ipMin=d1.IXPT2_COM[0]+1
       ipMax=d1.IXRB_COM[0]+1
       irMin=0
       irMax=d1.IYSPTRX2_COM[0]+1 ;-take extra point here
       ;
       for ip=ipMin,ipMax do begin
           qarr=var_out[ip, irMax:*]
           parr=d1.psi_com[ip, irMax:*,0]
           parr2=d1.psi_com[ip, irMin:irMax,0]
           qarr2=INTERPOL(qarr,parr,parr2)
           qarr2[*]=qarr[1] ;-extend with constant value
           var_out[ip, irMin:irMax]=qarr2
       endfor


       ;-extrapolate to outboard upper PF from SOL
       ipMin=d1.IXLB_COM[1]
       ipMax=d1.IXPT1_COM[1]
       irMin=0
       irMax=d1.IYSPTRX2_COM[0]+1 ;-take extra point here
       ;
       for ip=ipMin,ipMax do begin
           qarr=var_out[ip, irMax:*]
           parr=d1.psi_com[ip, irMax:*,0]
           parr2=d1.psi_com[ip, irMin:irMax,0]
           qarr2=INTERPOL(qarr,parr,parr2)
           qarr2[*]=qarr[1] ;-extend with constant value
           var_out[ip, irMin:irMax]=qarr2
       endfor


       ;-extrapolate to outboard lower PF from SOL
       ipMin=d1.IXPT2_COM[1]
       ipMax=d1.IXRB_COM[1]+1
       irMin=0
       irMax=d1.IYSPTRX2_COM[0]+1 ;-take extra point here
       ;
       for ip=ipMin,ipMax do begin
           qarr=var_out[ip, irMax:*]
           parr=d1.psi_com[ip, irMax:*,0]
           parr2=d1.psi_com[ip, irMin:irMax,0]
           qarr2=INTERPOL(qarr,parr,parr2)
           qarr2[*]=qarr[1] ;-extend with constant value
           var_out[ip, irMin:irMax]=qarr2
       endfor


       ;-extrapolate to inboard lower PF from SOL
       ipMin=d1.IXLB_COM[0]
       ipMax=d1.IXPT1_COM[0]
       irMin=0
       irMax=d1.IYSPTRX1_COM[0]+1 ;-take extra point here
       ;
       for ip=ipMin,ipMax do begin
           qarr=var_out[ip, irMax:*]
           parr=d1.psi_com[ip, irMax:*,0]
           parr2=d1.psi_com[ip, irMin:irMax,0]
           qarr2=INTERPOL(qarr,parr,parr2) ;-extend linearly
           qarr2[*]=qarr[1] ;-extend with constant value
           var_out[ip, irMin:irMax]=qarr2
       endfor

;-Finished PF regions


 ENDIF ELSE BEGIN   
;-Single-null case

       ;-core
       ipMin=d1.IXPT1_COM[0]+1
       ipMax=d1.IXPT2_COM[0]
       irMin=0
       irMax=d1.IYSPTRX1_COM[0]
       ;
       var_out[ipMin:ipMax, irMin:irMax]=$
         Int_subd(theta, var_in, theta[gjy0,0], $
                  ipMin, ipMax, irMin, irMax, order=-1)


       ;-sol
       ipMin=d1.IXLB_COM[0]
       ipMax=d1.IXRB_COM[0]+1
       irMin=d1.IYSPTRX1_COM[0]+1
       irMax=d1.NY_COM+1
       ;
       var_out[ipMin:ipMax, irMin:irMax]=$
         Int_subd(theta, var_in, theta[gjy0,0], $
                  ipMin, ipMax, irMin, irMax, order=-1)       


       ;-PF, outer
       ipMin=d1.IXPT2_COM[0]+1
       ipMax=d1.IXRB_COM[0]+1
       irMin=0
       irMax=d1.IYSPTRX1_COM[0]
       

       if (d1.IYSPTRX1_COM[0] lt d1.NY_COM) then begin
           ;-extrapolate to PF from SOL
           print, 'Extrapolating to inner PF region'
           for ip=ipMin,ipMax do begin
               qarr=var_out[ip,*]
               ReflectProfile, qarr, irMax
               var_out[ip,*]=qarr
           endfor
       endif


       ;-PF, inner
       ipMin=d1.IXLB_COM[0]
       ipMax=d1.IXPT1_COM[0]
       irMin=0
       irMax=d1.IYSPTRX1_COM[0]
       

       if (d1.IYSPTRX1_COM[0] lt d1.NY_COM) then begin
           ;-extrapolate to PF from SOL
           print, 'Extrapolating to inner PF region'
           for ip=ipMin,ipMax do begin
               qarr=var_out[ip,*]
               ReflectProfile, qarr, irMax
               var_out[ip,*]=qarr
           endfor
       endif

   ENDELSE

;
;
;
end



pro Calc_Curv, d1, qsafe, qinty, sinty, theta_all, bxcvx, bxcvy, bxcvz, DEBUG=DEBUG
;
; Calculate curvature-related quantities
;--------------------------------------------------------------------;


    CURVATURE, d1.nx_com+2, d1.ny_com+2, $
               d1.rm_com, d1.zm_com, d1.br_com, d1.bz_com, d1.bphi_com, d1.psi_com,$
               Theta_all,BXCV=BXCV

    ;-since dpsi was taken by absolute value we need to do same with psixy
    if (min(DERIV(d1.psi_com[5,*,0])) lt 0.) then sign=-1 else sign=1

    ;-calculate coefficients in front of partial derivatives in the curvature terms
    bxcvx=sign*BXCV.psi
    bxcvy=BXCV.theta
    bxcvz=sign*(BXCV.phi - sinty*BXCV.psi - qsafe*BXCV.theta)

if keyword_set(DEBUG) then STOP
end



function SetStructure, $
  d1, d2, d3, gjy0, R0, hthe0,iNixnorm, jNixnorm, Te_x, Ti_x, Ni_x, Vi_x, $
  bmag, simagxg, sibdryg, Rxy, Zxy, Bpxy, Btxy, Bxy, psixy,$
  hthe, dlthe, qsafe, kapanc, kapagc1, kapagc2, Rthe, Zthe, qinty, sinty,$
  Theta, dTheta, dx, Zpsi,Rpsi, dpsi,bxcvx, bxcvy, bxcvz,$
  Ni0, Vi0, Te0, Ti0, Nn0, kappaN, kappaV, kappaTe, kappaTi,$
  phi0, Ve0, yyc, Jpar0, INFO=INFO

if not keyword_set(INFO) then INFO='no comment available'

IF (d1.nxpt_com lt 2) then begin
;-Single-null case-;

    print, ''
    print, 'Setting data structure for single-null case'

    midpoint = LONG(d1.nx_com/2)       ;-use symmetric point for a single null
    nyue=d1.ny_com 
    nxue=d1.nx_com 

    data={$
         AREADME:info,$
         $
         nx:nyue,$ ;-number of radial points (without guards)        
         ny:nxue,$ ;-number of poloidal points (without guards)        
         $
         ixseps1:d1.IYSPTRX1_COM[0],$ ;-radial location right outside of the separatrix
         ixseps2:d1.IYSPTRX1_COM[0],$ ;-radial location right outside of the separatrix
         $
         jyseps1_1:LONG(d1.ixpt1_com[0]-1),$ ;-poloidal location right below the inner branch cut
         jyseps2_2:LONG(d1.ixpt2_com[0]-1),$ ;-poloidal location right below the inner branch cut
         jyseps1_2:midpoint,$
         jyseps2_1:midpoint,$
         $
         ixlb2:LONG(midpoint + 2),$
         $
         gjy0:gjy0-1,$  ;-shift because guard cells are removed
         R0:R0,$
         hthe0:hthe0,$
         $
         iNixnorm:0,$      ;-radial
         jNixnorm:gjy0-1,$ ;-poloidal
         $
         Te_x:Te_x,$
         Ti_x:Ti_x,$
         Ni_x:Ni_x,$
         Vi_x:Vi_x,$
         $
         bmag:bmag,$
         simagxg:d2.SIMAGXG_GRD,$
         sibdryg:d2.SIBDRYG_GRD,$
         $
         Rxy:TRANSPOSE(d1.rm_com[1:nxue,1:nyue,0]),$
         Zxy:TRANSPOSE(d1.zm_com[1:nxue,1:nyue,0]),$
         Bpxy:TRANSPOSE(d1.bpol_com[1:nxue,1:nyue,0]),$
         Btxy:TRANSPOSE(d1.bphi_com[1:nxue,1:nyue,0]),$
         Bxy:TRANSPOSE(d1.b_com[1:nxue,1:nyue,0]),$
         psixy:TRANSPOSE(d1.psi_com[1:nxue,1:nyue,0]),$
         $
         hthe:TRANSPOSE(hthe[1:nxue,1:nyue]),$
         dlthe:TRANSPOSE(dlthe[1:nxue,1:nyue]),$
         q_safe:TRANSPOSE(qsafe[1:nxue,1:nyue]),$
         kappa_n:TRANSPOSE(kapanc[1:nxue,1:nyue]),$
         kappa_g1:TRANSPOSE(kapagc1[1:nxue,1:nyue]),$
         kappa_g2:TRANSPOSE(kapagc2[1:nxue,1:nyue]),$
         $
         Rthe:TRANSPOSE(Rthe[1:nxue,1:nyue]),$
         Zthe:TRANSPOSE(Zthe[1:nxue,1:nyue]),$
         qinty:TRANSPOSE(qinty[1:nxue,1:nyue]),$
         sinty:TRANSPOSE(sinty[1:nxue,1:nyue]),$
         Theta:TRANSPOSE(Theta[1:nxue,1:nyue]),$
         $
         dy:REPLICATE(dTheta,nyue,nxue),$
         dx:TRANSPOSE(1./d2.gy_com[1:nxue,1:nyue]),$
         Zpsi:TRANSPOSE(Zpsi[1:nxue,1:nyue]),$
         Rpsi:TRANSPOSE(Rpsi[1:nxue,1:nyue]),$
         dpsi:TRANSPOSE(dpsi[1:nxue,1:nyue]),$
         $
         bxcvx:TRANSPOSE(bxcvx[1:nxue,1:nyue]),$
         bxcvy:TRANSPOSE(bxcvy[1:nxue,1:nyue]),$
         bxcvz:TRANSPOSE(bxcvz[1:nxue,1:nyue]),$
         $
         Ni0:TRANSPOSE(d3.NI___1__VALUE[1:nxue,1:nyue])/1e20,$
         Vi0:TRANSPOSE(d3.UP___1__VALUE[1:nxue,1:nyue]),$
         Te0:TRANSPOSE(d3.TE____EV_VALUE[1:nxue,1:nyue]),$
         Ti0:TRANSPOSE(d3.TI____EV_VALUE[1:nxue,1:nyue]),$
         Nn0:TRANSPOSE(d3.NN___1__VALUE[1:nxue,1:nyue])/1e20,$
         $
         kappaN:TRANSPOSE(kappaN[1:nxue,1:nyue]),$
         kappaV:TRANSPOSE(kappaV[1:nxue,1:nyue]),$
         kappaTe:TRANSPOSE(kappaTe[1:nxue,1:nyue]),$
         kappaTi:TRANSPOSE(kappaTi[1:nxue,1:nyue]),$
         phi0:TRANSPOSE(d3.PHI____VALUE[1:nxue,1:nyue]),$
         VE0:TRANSPOSE(Ve0[1:nxue,1:nyue]),$
         $
         x_array:d2.yyc_com[1:nyue],$
         Jpar0:fltarr(nyue,nxue)}

ENDIF



return, data
end



pro read_uedata3, data, path=path, save=save, noreformat=noreformat, info=info,$
debug=debug, filename=filename, NOPLOTS=NOPLOTS,nopdb = nopdb
;
;
; Calculate data for uedge.grd file
;-------------------------------------------------------------------;

  if not keyword_set(filename) then filename='uedge.grd.pdb'


;-import data from UEDGE data files
  if keyword_set(nopdb) then begin
     if not keyword_set(PATH) then PATH='.' 
     file1=PATH + '/gridue.pdb' 
     file2=PATH + '/uedgegrd.pdb'
     file3=PATH + '/uedgeout.pdb'
     restore,file1
     restore,file2
     restore,file3
  endif else begin
     ImportUE, d1, d2, d3, PATH=PATH, noreformat=noreformat
  endelse

   if not keyword_set(NOPLOTS) then begin
;-show the grid layout
      ShowGrid, d1, d2, d3, /FILL, title='Grid layout'
  wait, 3
endif


;-define Theta coordinate on the grid
  Get_Theta, d1, d2, d3, theta, dTheta, theta_all
if not keyword_set(NOPLOTS) then begin
  ShowGrid, d1, d2, d3, /FILL, var=theta, title='Theta'
  wait, 3
endif

;-calculate partial derivs dR/dTheta, dR/dPsi, dZ/dTheta, dZ/dPsi
  Get_Jacobian, $
                d1, d2, d3, dTheta, Rthe, Zthe, Rpsi, Zpsi, $
                dpsi, dlthe, hthe, qsafe, kapagc1, kapagc2, kapanc


;-select reference location and plasma parameters
  Get_Refs, d1, d2, d3, hthe, gjy0, gjy1, R0, hthe0, Ni_x, Te_x, Ti_x, Vi_x, $
            bmag, iNixnorm, jNixnorm


;-calculate derivatives d/dPsi
  Get_Psi_Derivs, d1, d2, d3, qsafe, dqdpsi, kappaN, kappaV, kappaTe, kappaTi, VE0
  

;-calculate qinty
  Int_dTheta, d1, d2, d3, gjy0, gjy1, theta, qsafe, qinty
  ;-adjust integration constant to have zero at inner branch cut
  qq=qinty
  for ip=0,d1.IXRB_COM[0]+1 do begin 
      qq[ip,*]=qinty[ip,*]-qinty[d1.IXPT1_COM[0]+1,*]
  endfor
  qinty=qq
if not keyword_set(NOPLOTS) then begin
  ShowGrid, d1, d2, d3, /FILL, var=qinty, title='qinty'
  wait, 3
endif


;-calculate sinty
  Int_dTheta, d1, d2, d3, gjy0, gjy1, theta, dqdpsi, sinty
if not keyword_set(NOPLOTS) then begin
  ShowGrid, d1, d2, d3, /FILL, var=sinty, title='sinty'
  wait, 3
endif


;-calculate curvature-related quantities
  Calc_Curv, d1, qsafe, qinty, sinty, theta_all, bxcvx, bxcvy, bxcvz


;-show contours of plasma background profiles
if not keyword_set(NOPLOTS) then begin
  ShowGrid, d1, d2, d3, /FILL, var=d3.NI___1__VALUE, tit="Ni, 1/m^3"
  wait,3

  ShowGrid, d1, d2, d3, /FILL, var=d3.TE____EV_VALUE, tit="Te, eV"
  wait,3

  ShowGrid, d1, d2, d3, /FILL, var=d3.TI____EV_VALUE, tit="Ti, eV"
  wait,3

  ShowGrid, d1, d2, d3, /FILL, var=d3.Nn___1__VALUE, tit="Nn, 1/m^3"
  wait,3

  ShowGrid, d1, d2, d3, /FILL, var=d3.UP___1__VALUE, tit="Vi||, m/s"
  wait,3

  ShowGrid, d1, d2, d3, /FILL, var=d3.PHI____VALUE, tit="Phi, V"
  wait,3   
endif


;-form a structure containing all data
  data=SetStructure($
       d1, d2, d3, gjy0, R0, hthe0,iNixnorm, jNixnorm, Te_x, Ti_x, Ni_x, Vi_x, $
       bmag, simagxg, sibdryg, Rxy, Zxy, Bpxy, Btxy, Bxy, psixy,$
       hthe, dlthe, qsafe, kapanc, kapagc1, kapagc2, Rthe, Zthe, qinty, sinty,$
       Theta, dTheta, dx, Zpsi,Rpsi, dpsi,bxcvx, bxcvy, bxcvz,$
       Ni0, Vi0, Te0, Ti0, Nn0, kappaN, kappaV, kappaTe, kappaTi,$
       phi0, Ve0, yyc, Jpar0, INFO=INFO)



;-save the data in a PDB file
  if keyword_set(SAVE) then begin
    ;Grd2pdb, data=data, output=filename

     

     
     names = ["AREADME","nx",$
         "ny","ixseps1","ixseps2","jyseps1_1","jyseps2_2","jyseps1_2","jyseps2_1","ixlb2",$
         "gjy0","R0","hthe0","iNixnorm","jNixnorm","Te_x","Ti_x","Ni_x","Vi_x",$
         "bmag","simagxg","sibdryg","Rxy","Zxy","Bpxy","Btxy","Bxy","psixy","hthe","dlthe",$
            "q_safe",$
            "kappa_n","kappa_g1","kappa_g2","Rthe","Zthe","qinty","sinty","Theta","dy","dx",$
            "Zpsi","Rpsi","dpsi","bxcvx","bxcvy","bxcvz","Ni0","Vi0","Te0","Ti0", "Nn0",$
            "kappaN","kappaV","kappaTe","kappaTi","phi0","VE0","x_array","Jpar0"]
     
     grdpdb = file_export(filename,data,custom_name= names)

  

     
  endif

;
;
;
if keyword_set(DEBUG) then STOP
end
