;hlmk_grids without paramaersl

;hlmk_grids,/simple,/narrow,/small_y_res

;for example to create  grid to analyze radially localized modes
;hlmk_grids,/simple,/narrow,/local_r

;to create simple local grid with only Bz field
;hlmk_grids,/simple,/narrow,/local_r,Bz0=.1,bphi0=0.0,gridname =
;'Helimak_Bz'

;create a grid with a super weak density gradient
;hlmk_grids,/simple,/narrow,/local_r,Bz0=.1,bphi0=0.0,gridname
;='Helimak_Bz'

;build a 5x64 grid
;hlmk_grids,/simple,/narrow,/local_r,Bz0=.1,bphi0=0.0,gridname
;='Helimak_Bz',grid_size = 6 

;hlmk_grids,/simple,/narrow,/local_r,Bz0=.1,bphi0=.1,gridname='Helimak_bz_1_1',grid_size=6 

;hlmk_grids,/simple,/narrow,/local_r,Bz0=.1,bphi0=0.0,gridname='Helimak_bz',grid_size=5

;messing with radial comp, local r but not narrow just to see where
;Rxy is used in bOUT++
;hlmk_grids,/simple,/local_r,Bz0=.1,bphi0=0.0,gridname='Helimak_bz_local_wide',grid_size=5



;create a cold dense plasma
;hlmk_grids,/simple,/narrow,/local_r,Bz0=.1,bphi0=0.0,Te=2,Ni0 =
;1e18,grid_size = 5,gridname ='Helimak_Bz_COLD+DENSE'

;hlmk_grids,/simple,/local_r,Bz0=.1,bphi0=0.0,Te=2,Ni0 = 1e18,grid_size = 5,gridname ='Helimak_Bz_COLD+DENSE+WIDE2'

;create a cold plasma
;hlmk_grids,/simple,/narrow,/local_r,Bz0=.1,bphi0=0.0,Te=2,Ni0 =
;1e16,grid_size = 5,gridname ='Helimak_Bz_COLD'

;create a cold plasma with radial resolve
;hlmk_grids,/simple,Bz0=.1,bphi0=0.0,Te=2,Ni0 = 1e18,grid_size = 5,gridname ='Helimak_Bz_COLD+DENSE+NARROW',/narrow 




pro hlmk_grids,full=full,Lc = Lc, $
               Ln = Ln, Lphi = Lphi,Lte = Lte,$
                 grid_size= N,Te0 = Te0,Ti0 = Ti0,$                 
               narrow= narrow,simple = simple,$
               small_y_res=small_y_res,$
               name= name,local_r = local_r,gridname=gridname,$
               bphi0 = bphi0, Bz0 = Bz0,Ni0 = Ni0
  
  ;to avoid making this thing too general I will myself to 
  ;assuming that this script will simply allow the user to tweak
  ;Bz(connection length), phi0V, ni and Te gradients,  
  

  ;we need a way to append some metadata to the grids rather than
  ;just creating long ugly filenames
  Zmax = 2.0
  
                                ;simple exp profiles, constant ExB
  if keyword_set(simple)then begin
     ni_profile_type=2 ;simple exp profile, lam_n ignored
     ti_profile_type=0
     te_profile_type=0
     phi_profile_type = 0
     if not(keyword_set(N))then N = 5
  endif
  
  if keyword_set(narrow) then begin
     rMin = .95
     rMax = 1.05
  endif else begin
     rMin= .90
     rMax = 1.1
  endelse


  
                                ;bphi is typically fixed
  
  if not keyword_set(bphi0) then begin
     print, 'seting bphi0 = 0.0'
     bphi0 = 0.0
  endif
  
  if not keyword_set(gridname) then gridname = 'Helimak' 
  
                                ;Bz0 specification overides Lc
  if keyword_set(Bz0) then begin
     if Bz0 EQ 'zero' then Bz0 = 0
     Btot = sqrt(Bz0^2 + bphi0^2)
     Lc = Zmax*Btot/Bz0
  endif else begin
     if not keyword_set(Lc)then begin ;no Bz0 or Lc
        Bz0 = bphi0/10.
        
        Btot = sqrt(Bz0^2 + bphi0^2)
        Lc = Zmax*Btot/Bz0
     endif else begin           ;Lc set but Bz0 is not
        
        Bz0 = bphi0/sqrt((Lc/Zmax)^2 -1) 
     endelse
  endelse            ;Bz0 set but no LC
  
  Btot = sqrt(Bz0^2 + bphi0^2)
  Lc = Zmax*Btot/Bz0
 


  if (not keyword_set(N)) then N = 4 
  
  

;Nz here is in fact the number of grid point ALONG the field lines,
;this code was originally used to create grids for machines with
;simple cylindrical geometries, where Bz = Bpar


  if keyword_set(small_y_res)then begin
      Nz = 2^(N-1)
      Nr = 2^(N) + 4
  endif else begin
     Nz = 2^N
     Nr = 2^N + 4
  endelse

  ;local r modes only
  if keyword_set(local_r) then begin
     Nr =  1 + 4
     Nz = 2^(N)
  endif
  temp = [gridname+"_",string(Nr-4),"x",string(Nz),"_",string(Lc),".nc"]

  if (not keyword_set(filename)) then filename = strcompress(strjoin(temp),/remove_all)

                                ;if the slopes are not indicated set
                                ;them to be on par with the system
                                ;size, but small enough to keep all
                                ;quantaties positive


  if not keyword_set(slope_te) then slope_te_amp = 1 /(rMax - rMin)
  if not keyword_set(slope_n) then slope_n_amp = .5 /(rMax - rMin)
  if not keyword_set(slope_ti) then slope_ti_amp = 1 /(rMax - rMin)
  
  
  ;ni0 = 5e16 ;units? set_mesh_cyl will want m^-3
  
  if not keyword_set(Te0) then Te0 = 10.0 
  if not keyword_set(Ni0) then Ni0 = 5e16
  
  ;we we create 10 grid with different Te gradient
  for i=0,8 do begin
     slope_n = (i-4.)/4. * slope_n_amp 
     slope_te = 0.0
     slope_ti = 0.0
     lam_n = (10+(i-4.)/1.0)/100
     ;lam_n = (10+(i-4.)/1.0)/1000
     te0 = Te0
     ni0 = Ni0
     print,"lam_n: ",lam_n
     
     
     temp = [gridname,"_",string(Nr-4),"x",string(Nz),"_",sigfig(lam_n,2,sci = sci),"_lam_n.nc"]

     filename = strcompress(strjoin(temp),/remove_all)

     
     set_mesh_cyl,/export,Nr = Nr, Nz = Nz,rMin = rMin, rMax = rMax,ni0 =ni0 $
                  ,te0=te0,Bz0 = Bz0,bphi0 = bphi0,Zmax=Zmax,$
                  ni_profile_type = ni_profile_type,ti0 =ti0,$
                  te_profile_type = te_profile_type,$
                  ti_profile_type = ti_profile_type,phi_profile_type = phi_profile_type,$            
                  slope_te = slope_te, slope_ti = slope_ti,$
                  slope_n = slope_n,lam_n = lam_n
     read_uedata3, /s, d, /noref, /NOPLOTS,/nopdb,filename = filename
     spawn,"rm *.pdb"

  endfor
;this script with generate a helimak grid witha 
  
  

;;   plot,(shot_data.set3.vfloat)[*,0]
;;   oplot,(shot_data.set3.vfloat)[*,1]
  
;;   oplot,(shot_data.set3.vfloat)[*,2]
;;   oplot,(shot_data.set3.vfloat)[*,3]
;;   oplot,(shot_data.set3.vfloat)[*,4]
;;   oplot,(shot_data.set3.vfloat)[*,5]
  
;;   oplot,(shot_data.set3.vfloat)[*,6]
  
  
end
