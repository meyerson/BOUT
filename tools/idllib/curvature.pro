function pdiff_rz, r, z, f
;
; Calculate partial derivatives df/dr and df/dz
; for function f given on a set of 5 points (r,z)
; using singular-value decomposition
;
; Inputs: arrays r[5],z[5],f[5]
; Output: structure {r:df/dr, z:df/dz}
;-------------------------------------------------

   A=TRANSPOSE([[fltarr(5)+1],[r-r(0)],[z-z(0)]])

   SVDC, A,W,U,V

   res=SVSOL(U,W,V,f)

   pdiff={r:res[1],z:res[2],phi:0.0}

;
;
;
return, pdiff
end


function curlcyl, vecR, vecV, gradVr, gradVphi, gradVz
;
; Calculate curl of a axisymmetric vector field V
; in cylindrical coords
;
; Inputs: 
;        vecR - location vector in cylindrical components {r:r,z:z}
;        vecV - vector V in cylindrical components {r:Vr,phi:Vphi,z:Vz} 
;        gradVr - gradient of the r-component,     {dVr/dr,dVr/dz}
;        gradVphi - gradient of the phi-component, {dVphi/dr,dVphi/dz}
;        gradVz - gradient of the z-component,     {dVphi/dr,dVphi/dz}
;
; Output: curl in cylindrical coordinates
;-------------------------------------------------


  curl={r:-gradVphi.z, phi:gradVr.z-gradVz.r, z:vecV.phi/vecR.r-gradVphi.r}
;
;
;
return, curl
end


function xprod, v1, v2, minus=minus
;
; Calculate cross-product of two vectors
; in cylindrical coordinates
;
; Inputs:
;        v1={r,phi,z}
;        v2={r,phi,z}
;
; Output: v1xv2 {r,phi,z}
;---------------------------------------


    r = v1.phi*v2.z-v1.z*v2.phi
    phi = v1.z*v2.r-v1.r*v2.z
    z = v1.r*v2.phi-v1.phi*v2.r

;
 if keyword_set(MINUS) then begin
   res={r:-r,phi:-phi,z:-z} 
 endif else begin
   res={r:r,phi:phi,z:z}
 endelse

return, res
end


function dotprod, v1, v2
;
; Calculate dot-product of two vectors
; in cylindrical coordinates
;
; Inputs:
;        v1={r,phi,z}
;        v2={r,phi,z}
;
; Output: (v1,v2)
;---------------------------------------

    res=v1.r*v2.r + v1.phi*v2.phi + v1.z*v2.z

return, res
end


pro curvature, nx, ny, Rxy, Zxy, BRxy, BZxy, BPHIxy, PSIxy, THETAxy,$
CURLB=CURLB, JXB=JXB, CURVEC=CURVEC, BXCURVEC=BXCURVEC, BXCV=BXCV,$
DEBUG=DEBUG
;
; Calculate the magnetic field curvature and other related quantities
;--------------------------------------------------------------------

print, 'Calculating curvature-related quantities...'

;;-vector quantities are stored as 2D arrays of structures {r,phi,z}
   vec={r:0.,phi:0.,z:0.}
   curlb=REPLICATE(vec,nx,ny) 
   jxb=REPLICATE(vec,nx,ny) 
   curvec=REPLICATE(vec,nx,ny) 
   bxcurvec=REPLICATE(vec,nx,ny)

   vec2={psi:0.,theta:0.,phi:0.}
   bxcv=REPLICATE(vec2,nx,ny)



          FOR i=0,nx-1 DO BEGIN
           FOR j=0,ny-1 DO BEGIN
             ;
              grad_Br=pdiff_rz(REFORM(Rxy[i,j,*]), REFORM(Zxy[i,j,*]), REFORM(BRxy[i,j,*]))
              grad_Bz=pdiff_rz(REFORM(Rxy[i,j,*]), REFORM(Zxy[i,j,*]), REFORM(BZxy[i,j,*]))
              grad_Bphi=pdiff_rz(REFORM(Rxy[i,j,*]), REFORM(Zxy[i,j,*]), REFORM(BPHIxy[i,j,*]))

              grad_Psi=pdiff_rz(REFORM(Rxy[i,j,*]), REFORM(Zxy[i,j,*]), REFORM(PSIxy[i,j,*]))
              grad_Theta=pdiff_rz(REFORM(Rxy[i,j,*]), REFORM(Zxy[i,j,*]), REFORM(THETAxy[i,j,*]))
              grad_Phi={r:0.0,z:0.0,phi:1./Rxy[i,j,0]} ;-gradient of the toroidal angle

              vecR={r:Rxy[i,j,0],z:Zxy[i,j,0]}
              vecB={r:BRxy[i,j,0],z:BZxy[i,j,0],phi:BPHIxy[i,j,0]}

              curlb[i,j]=CurlCyl(vecR, vecB, grad_Br, grad_Bphi, grad_Bz)
              jxb[i,j]=Xprod(curlb[i,j], vecB)


              ;-magnitude of B at 5 locations in cell
              bStrength=SQRT(REFORM(BRxy[i,j,*])^2 + REFORM(BZxy[i,j,*])^2 + REFORM(BPHIxy[i,j,*])^2)

              ;-unit B vector at cell center
              vecB_unit={r:(BRxy[i,j,0])/bStrength[0], $
                          z:(BZxy[i,j,0])/bStrength[0], $
                           phi:(BPHIxy[i,j,0])/bStrength[0]}

                   ;-components of gradient of unit B vector at 5 locations in cell
                   grad_Br_unit=pdiff_rz(REFORM(Rxy[i,j,*]), $
                                         REFORM(Zxy[i,j,*]), $
                                          REFORM(BRxy[i,j,*]/bStrength))

                   grad_Bz_unit=pdiff_rz(REFORM(Rxy[i,j,*]), $
                                         REFORM(Zxy[i,j,*]), $
                                          REFORM(BZxy[i,j,*]/bStrength))

                   grad_Bphi_unit=pdiff_rz(REFORM(Rxy[i,j,*]), $
                                           REFORM(Zxy[i,j,*]), $
                                            REFORM(BPHIxy[i,j,*]/bStrength))

              ;-curl of unit B vector at cell center
              curlb_unit=CurlCyl(vecR, vecB_unit, grad_Br_unit, grad_Bphi_unit, grad_Bz_unit)

              ;-curvature vector at cell center
              curvec[i,j]=Xprod(vecB_unit,curlb_unit,/MINUS)

              ;-unit b cross curvature vector at cell center
              bxcurvec[i,j]=Xprod(vecB_unit,curvec[i,j])


              ;-calculate bxcurvec dotted with grad_psi, grad_theta, and grad_phi
              bxcv[i,j].psi=Dotprod(bxcurvec[i,j],grad_Psi)
              bxcv[i,j].theta=Dotprod(bxcurvec[i,j],grad_Theta)
              bxcv[i,j].phi=Dotprod(bxcurvec[i,j],grad_Phi)

             ;
           ENDFOR
          ENDFOR

if keyword_set(DEBUG) then STOP
print, '...done'
;
;
end
lb_unit,/MINUS)

              ;-unit b cross curvature vector at cell center
              bxcurvec[i,j]=Xprod(vecB_unit,curvec[i,j])


              ;-calculate bxcurvec dotted with grad_psi, grad_theta, and grad_phi
              bxcv[i,j].psi=Dotprod(bxcurvec[i,j],grad_Psi)
              bxcv[i,j].theta=Dotprod(bxcurvec[i,j],grad_Theta)
              bxcv[i,j].phi=Dotprod(bxcurvec[i,j],grad_Phi)

         