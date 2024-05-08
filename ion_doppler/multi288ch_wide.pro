PRO Multi288CH_wide
!P.multi=0 & device,decomposed=0 & loadct,39 &!P.color=0 & !P.background=255
a= findgen(4) * (!PI*2.0/4.0) & USERSYM, cos(A), sin(A), /FILL
tmax=100
EM_max=1.e8
lambda0=480.6;471.3;656.3;480.6;587.562;486.133;656.3;486.133;468.57;486.133;nm
mass=39.95;4.;12.01
resolution=-1.420387E-12*lambda0^3 - 2.156031E-09*lambda0^2 + 1.250038E-06*lambda0 + 3.830769E-03;0.0037714; -0.000000000001420387*lambda0^3 - 0.000000002156031*lambda0^2 + 0.000001250038*lambda0 + 0.003830769
z=read_ascii("C:\Users\haruaki_tanaka\Documents\ion-doppler\z_negative.txt") & z=reform(z.field1)*1.e-3
p=read_Ascii("C:\Users\haruaki_tanaka\Documents\ion-doppler\r.txt") & p=p.field1*1.e-3
edge=0.33     ; 2022/7　解析～　ポテンシャルの範囲より仮定
filename=dialog_pickfile(path="\\192.168.1.111\experiment\results\Doppler\Andor\320CH\20210924",filter="*.asc")
;filename=dialog_pickfile(path="C:\Users\haruaki_tanaka\Desktop\20210128",filter="*.asc")

d=read_ascii(filename[0]) & d=transpose(d.field0001(1:*,0:*))
bg=read_Ascii("\\192.168.1.111\experiment\results\Doppler\Andor\320CH\20210916\bg.asc") & bg=transpose(bg.field0001(1:*,0:*)) & D=D-Bg
calib=read_ascii("C:\Users\haruaki_tanaka\OneDrive\ドキュメント\ion-doppler\Ar_calibration.0916_remake.txt",data_start=1) & calib=calib.field1
ch=reform(calib[0,*])
Center=reform(calib[1,*])
Smile=reform(calib[2,*])
relative=reform(calib[4,*])
instru=reform(calib[5,*])
Ti_instru=1.69e8*mass*(2.*resolution*instru*sqrt(2.*alog(2.))/lambda0)^2

separation=where((CH-1) mod 16 eq 1)-1 & separation[3]=separation[3]+1
;Only since late in 2021 (Effective when CH 1 is deleted)***************************************************************************
;separation = [0, separation]                                                              
;***********************************************************************************************************************************
lambda=dblarr([1024,n_elements(CH)])
x=dindgen(1024)
for i=0,n_elements(CH)-1 do begin
lambda[*,i]=(x-smile[i])*resolution+lambda0;+0.13
endfor

spectra=fltarr([1024,n_elements(CH)])
!p.MULTI=[0,1,2]
window,0,xsize=800,ysize=1000
loadct,1
interval=7
contour,d,x,x,/fill,nlevels=16,xst=1,yst=1,zst=1,zr=[min(d),(max(d)-min(d))*0.1+min(d)],yr=[min(center)-interval-5,max(center)+interval+5],position=[0.1,0.4,0.95,0.95]
loadct,39
for i=0,n_elements(CH)-1 do begin
oplot,smile[i]+[-40,40],center[i]*[1,1],color=250 & oplot,smile[i]*[1.,1.],center[i]+[-7,7],color=100 &
oplot,smile[i]+[-40,40],center[i]*[1,1]+interval,color=200,linestyle=2 & oplot,smile[i]+[-40,40],center[i]*[1,1]-interval,color=200,linestyle=2
endfor

xr=[410,540]
lambdaA=(x-(xr[0]+xr[1])/2.)*resolution+lambda0 & lambdaB=lambdaA & lambdaA=lambdaA(xr[0]:xr[1])
ybin=TOTAL(d,2)
PLOT,x,ybin,XST=1,YST=1,PSYM=-1,position=[0.1,0.05,0.95,0.4]
oplot,x(xr[0]:xr[1]),gaussfit(x[xr[0]:xr[1]],ybin[xr[0]:xr[1]],coeff,nterm=5),color=250,thick=2
A=[coeff[0],200,coeff[2],coeff[0],350,coeff[2],coeff[0],500,coeff[2],coeff[0],650,coeff[2],coeff[0],800,coeff[2],coeff[3],0,0]
fita=A*0.+1.

;**************************
;*Line-integrated analysis*
;**************************

passive_Ti=CH*0.
passive_Timax=CH*0.
passive_Timin=CH*0.
passive_Em=CH*0.

!P.multi=[0,4,4]
window,2,xsize=1800,ysize=900
for i=0,n_elements(CH)-1 do begin
  for j=0,n_elements(X)-1 do begin
    spectra[j,i]=spectra[j,i]+total(d(j:j,center[i]-interval:center[i]+interval))*relative[i]         ;各チャンネルのガウシアン信号を積分
  endfor
  
  input=reform(spectra[*,i]) 
  fit=gaussfit(x(round(smile[i])-75:round(smile[i])+75),input(round(smile[i])-75:round(smile[i])+75),coeff,nterms=4,sigma=sigma)
  passive_Ti[i]=1.69e8*mass*(2.*resolution*(coeff[2])*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru[i]
  passive_Timax[i]=1.69e8*mass*(2.*resolution*(coeff[2]+sigma[2])*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru[i]
  passive_Timin[i]=1.69e8*mass*(2.*resolution*(coeff[2]-sigma[2])*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru[i]
  if i+1 mod 16 eq 0 then window,(i+1)/16+2,xsize=1000,ysize=900
  plot,x(round(smile[i])-75:round(smile[i])+75),input(round(smile[i])-75:round(smile[i])+75),xst=1,yst=1,title="CH::" + strcompress(CH[i]),charsize=2,psym=1
  oplot,x(round(smile[i])-75:round(smile[i])+75),fit,color=250
  passive_Em[i]=resolution*total(input(round(smile[i])-75:round(smile[i])+75)-min(smooth(reform(input(round(smile[i])-75:round(smile[i])+75)),20)))
  spectra[*,i]=spectra[*,i]-min(smooth(reform(input(round(smile[i])-75:round(smile[i])+75)),20))
  
endfor

;window,24 & !P.multi=0
;plot,CH,passive_Ti,psym=8,xst=1,yst=1,yr=[0,max(passive_Ti)]
;errplot,CH,passive_Timin,passive_Timax
;for i=0,n_elements(separation)-1 do oplot,CH[separation[i]]*[1,1],[-10000,10000],color=150

Ti2D=fltarr([n_elements(p),n_elements(z)])
Em2D=fltarr([n_elements(p),n_elements(z)])
for i=0,n_elements(z)-2 do begin
  if i mod 2 eq 0 then Ti2D[*,i]=spline(CH(separation[i]:separation[i+1]-1)-CH(separation[i]),passive_Ti(separation[i]:separation[i+1]-1),indgen(16))
  if i mod 2 eq 1 then Ti2D[*,i]=spline(CH(separation[i]:separation[i+1]-1)-CH(separation[i]),passive_Ti(separation[i]:separation[i+1]-1),indgen(16))
  if i mod 2 eq 0 then Em2D[*,i]=spline(CH(separation[i]:separation[i+1]-1)-CH(separation[i]),passive_Em(separation[i]:separation[i+1]-1),indgen(16))
  if i mod 2 eq 1 then Em2D[*,i]=spline(CH(separation[i]:separation[i+1]-1)-CH(separation[i]),passive_Em(separation[i]:separation[i+1]-1),indgen(16))
if i mod 2 eq 1 then Ti2D[*,i]=reverse(Ti2D[*,i]) 
if i mod 2 eq 1 then Em2D[*,i]=reverse(Em2D[*,i])

endfor

Ti2D[*,n_elements(z)-1]=spline(CH(separation[i]:*)-CH(separation[i]),passive_Ti(separation[i]:*),indgen(16)) & Ti2D[*,n_elements(z)-1]=reverse(Ti2D[*,n_elements(z)-1])
Em2D[*,n_elements(z)-1]=spline(CH(separation[i]:*)-CH(separation[i]),passive_Em(separation[i]:*),indgen(16)) & Em2D[*,n_elements(z)-1]=reverse(Em2D[*,n_elements(z)-1])
Em2D[*,3]=spline(CH(separation[3]:separation[4]-1)-CH(separation[3])+1,passive_Em(separation[3]:separation[4]-1),indgen(16)) & Em2D[*,n_elements(z)-1]=reverse(Em2D[*,n_elements(z)-1]) & Em2D[*,3]=reverse(Em2D[*,3])
window,25,xsize=1000, ysize=800 & !P.multi=[0,2,1]
contour,Ti2D/max(Ti2D),p,z,/fill,nlevels=256,xst=1,yst=1,zst=1,/isotropic,zr=[0,1]
loadct,1 & contour,Em2D/(max(Em2D)*0.75),p,z,/fill,nlevels=256,xst=1,yst=1,zst=1,/isotropic,zr=[0,1];max(Em2D)*0.75]
loadct,39
!P.multi=0

;==================relative check ==============================================================================
;i=12
;relative_1d = fltarr(16) & relative_1d=spline(CH[separation[i]:separation[i+1]-1]-CH[separation[i]],relative(separation[i]:separation[i+1]-1),indgen(16)) 
;if i mod 2 eq 1 then relative_1d=reverse(relative_1d)
;window,18,xsize=1000,ysize=800 & !p.multi=[0,1,2]
;plot,p,Em2D[*,i],psym=-1,title="Em2D"
;plot,p,relative_1d,psym=-1,title="relative"
;******************
;*Abel-inversion*
;******************
num=250
yy=findgen(num+1)/num*(edge-min(p))+min(p)
spectra_interp=fltarr([n_elements(lambdaA),n_elements(yy),n_elements(z)])

;#####16bundole.ver####################################
for i=0,n_elements(z)-2 do begin
window,i;,xsize=400,ysize=300
  if i mod 2 eq 0 then result=trigrid_interpor_for_r_lambda([[reform(spectra(0:*,separation[i]:separation[i+1]-1))],[fltarr(n_elements(lambdaB))]],[[reform(lambda[0:*,separation[i]:separation[i+1]-1])],[lambdaB]],[p(CH(separation[i]:separation[i+1]-1)-CH(separation[i])),edge],lambdaA,yy)
  if i mod 2 eq 1 then result=trigrid_interpor_for_r_lambda([[fltarr(n_elements(lambdaB))],[reform(spectra(0:*,separation[i]:separation[i+1]-1))]],[[lambdaB],[reform(lambda[0:*,separation[i]:separation[i+1]-1])]],[edge,reverse(p(CH(separation[i]:separation[i+1]-1)-CH(separation[i])))],lambdaA,yy)
  spectra_interp[*,*,i]=result.z 
endfor
window,17
result=trigrid_interpor_for_r_lambda([[fltarr(n_elements(lambdaB))],[reform(spectra(0:*,separation[17]:*))]],[[lambdaB],[reform(lambda[0:*,separation[17]:*])]],[edge,reverse(p(CH(separation[17]:N_ELEMENTS(ch)-1)-CH(separation[17])))],lambdaA,yy)
spectra_interp[*,*,17]=result.z 

;**********2020 January****************************
;for i=0,130 do spectra_interp[i,*,13]=reverse(spectra_interp[i,*,13],2)

;*************************************************

Local_spectra=fltarr([n_elements(lambdaA),num+1,n_elements(z)])
spectra_interp=smooth(spectra_interp,[5,num/16,1])
dy=yy[1]-yy[0]

for k=0,n_elements(z)-1 do begin
  for l=0,n_elements(lambdaA)-1 do begin
 derivative=deriv(yy,reform(spectra_interp[l,*,k]))
    for i=0,num do begin
      for j=i,num-1 do local_spectra[l,i,k]=local_spectra[l,i,k]-1./!pi*derivative[j]*alog(yy[j+1]*(1.+sqrt(1.-yy[i]^2/yy[j+1]^2))/(yy[j]*(1.+sqrt(1.-yy[i]^2/yy[j]^2))));Balandin's Abel inversion
    endfor
  endfor
endfor
Local_Spectra=smooth(Local_Spectra,[5,num/16,1])

emission=fltarr([num+1,n_elements(z)])
Ti_2D=fltarr([num+1,n_elements(z)])
Ti_max=fltarr([num+1,n_elements(z)])
Ti_min=fltarr([num+1,n_elements(z)])
Ti_instru2=total(Ti_instru(separation[0]:separation[1]))/(separation[1]-separation[0])
for i=1,n_elements(z)-2 do Ti_instru2=[Ti_instru2,total(Ti_instru(separation[i]:separation[i+1]))/(separation[i+1]-separation[i])]
Ti_instru2=[Ti_instru2,total(Ti_instru(separation[i]:n_elements(CH)-1))/(n_elements(CH)-1-separation[i])]

!P.multi=[0,6,6]
for i=0,num do begin
  ;if i mod 36 eq 0 then window, fix(i/36),xsize=1700,ysize=900;,/free
  for j=0,n_elements(z)-1 do begin
  input=reform(Local_spectra[*,i,j]) 
  for l=0,130 do begin 
      if input[l] lt 0 then input[l] = -input[l] -min(abs((smooth(reform(Local_spectra[*,i,j]),20))))
  endfor 
  fit = gaussfit(lambdaA,input,coeff,nterms=3,sigma=sigma)
  ;plot,lambdaA,input,xst=1,yst=1,title="(R,Z) ="+strcompress(yy[i],/remove_all)+","+strcompress(z[j],/remove_all),charsize=1.5
  ;oplot,lambdaA,fit,color=250,thick=2
  Ti_2D[i,j]=1.69e8*MASS*(2.*coeff[2]*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru2[j]
  Ti_max[i,j]=1.69e8*MASS*(2.*(abs(coeff[2])+3.*abs(sigma[2]))*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru2[j]
  Ti_min[i,j]=1.69e8*MASS*(2.*(abs(coeff[2])-3.*abs(sigma[2]))*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru2[j]
  emission[i,j]=total(input*resolution);coeff[0]
  checker=float(abs(coeff[1]-lambda0) lt 0.1)*(coeff[0] gt 0)*(emission[i,j] gt 100)*float(abs(Ti_max[i,j]-Ti_min[i,j]) lt Ti_2D[i,j]);*float(abs(Ti_max[i,j]-Ti_min[i,j]) lt Ti_2D[i,j]+Ti_instru2[j]);*(emission_local[i,j] gt EM_max*0.1)
  Ti_2D[i,j]=Ti_2D[i,j]*checker & Ti_max[i,j]=Ti_max[i,j]*checker & Ti_min[i,j]=Ti_min[i,j]*checker
  endfor
endfor

;remove "NAN" grid point
for i=0,n_elements(z)-1 do begin
  for j=0,n_elements(yy)-1 do begin
    if finite(Ti_2D[j,i]) ne 1 then begin
    Ti_2D[j,i]=0.
    Ti_max[j,i]=0.
    Ti_min[j,i]=0.
    print,i,j   
    endif 
  endfor
endfor

start=where((yy gt p[0])) & start=start[0]
;Ti_2D=Ti_2D*float(emission_local gt 0)
;Ti_2D=Ti_2D*float(emission_local gt EM_max*0.04)
;for i=0,n_elements(z)-1 do Ti_2D[*,i]=Ti_2D[*,i]*float(abs(Ti_max[*,i]-Ti_min[*,i]) lt Ti_2D[*,i]+Ti_instru2[i])
checker=Ti_2D le 0
for i=start,n_elements(yy)-2 do begin
    for j=0,n_elements(z)-1 do begin
      if checker[i,j] then begin
        Ti_2D[i,j]=sqrt((Ti_2D[i+1,j]+Ti_2D[i-1,j])/2.)
        Ti_max[i,j]=(Ti_max[i+1,j]+Ti_max[i-1,j])/2.  
        Ti_min[i,j]=(Ti_min[i+1,j]+Ti_min[i-1,j])/2.        
      endif
    endfor  
endfor
loadct,39
window,26,ysize=800 & !P.multi=[0,2,1]
contour,Ti_2D/tmax,yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,/isotropic,zr=[0,1],xr=[min(p),max(p)]
loadct,1 
contour,emission/max(emission*0.75),yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,/isotropic,zr=[0,1],xr=[min(p),max(p)];max(Em2D)*0.75]

window,27,ysize=500,xsize=1400
!P.multi=[0,4,1]
loadct,39 & contour,transpose(Ti2D),z,p,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,tmax],xtitle="Z[m]",ytitle="R[m]",title="Ti(Projection)";,/isotropic
plot,Ti2D[*,0],p,xst=1,yst=1,xr=[0,tmax*2],/nodata,xtitle="Ti [eV]",ytitle="R[m]",title="Ti(Projection)" & for i=0,5 do oplot,Ti2D[*,i],p,color=i*50,psym=-8
loadct,39 & contour,transpose(Ti_2D),z,yy,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,tmax],xtitle="Z[m]",ytitle="R[m]",title="Ti(Local)",yr=[min(p),max(p)];,/isotropic
plot,Ti_2D[*,0],yy,xst=1,yst=1,xr=[0,tmax*2],/nodata,xtitle="Ti [eV",ytitle="R[m]",title="Ti(Local)",yr=[min(p),max(p)] & for i=0,17 do oplot,Ti_2D[*,i],yy,color=i*15,psym=-8

print,"The average Ti in downstream is" + string(mean(Ti_2D[41:120,6:16])) + "eV"
save,filename=strmid(filename,0,70)+".sav",Ti_2D,yy,z,emission
;epson,filename='C:\Users\haruaki_tanaka\Documents\ion-doppler\Doppler288ch\20200116\'+'Shot10.eps',aspect=aspect_ratio
;!p.multi=[0,3,1]
;loadct,39 &contour,Ti_2D/tmax,yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,1],/isotropic,ytitle="Z [m]",xtitle="R [m]",title="Ti (Local)",color=0,charsize=1.0,CHARTHICK=2.0
;loadct,5 & contour,emission/max(emission*0.75),yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,1],/isotropic,xtitle="R [m]",ytitle="Z [m]",title="Emission (Local)",color=0,charsize=1.0,CHARTHICK=2.0
;loadct,39 &contour,smooth(Ti_2D,[num/16,2])/tmax,yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,1],/isotropic,ytitle="Z [m]",xtitle="R [m]",title="Ti (Local)",color=0,charsize=1.0,CHARTHICK=2.0
;epsoff
stop

data=mag_haru(dateshot=210305009)
!P.background=255;16777215L
!P.color=0
time=8
r=data.c[0:86] & z_m=data.b[3:94] & Psi=data.A[3:94,0:86,time] & Bp=data.I[3:94,0:86,time] 

window,30 & !P.multi=[0,2,1] & loadct,39
contour,Ti_2D/tmax,yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,1],/isotropic,xtitle="R [m]",ytitle="Z [m]",title="Ti (Local)"+strmid(filename,38,6),xr=[p[0],max(p)],color=0,charsize=1.0,CHARTHICK=2.0
;for i=0,n_elements(z)-1 do oplot,fltarr(n_elements(p))+Z[i],p,psym=8,color=250
contour,transpose(Psi/max(abs(Psi))),r,z_m,/overplot,nlevels=20,color=0,thick=2.0
contour,-transpose(Psi/max(abs(Psi))),r,z_m,/overplot,nlevels=20,color=0,thick=2.0
loadct,5 & contour,emission/max(emission*0.75),yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,1],/isotropic,xtitle="R [m]",ytitle="Z [m]",title="Emission (Local)",xr=[p[0],max(p)],color=0,charsize=1.0,CHARTHICK=2.0
loadct,39;& for i=0,n_elements(z)-1 do oplot,fltarr(n_elements(p))+Z[i],p,psym=8,color=150
contour,transpose(Psi/max(abs(Psi))),r,z_m,/overplot,nlevels=20,color=0,thick=2.0
contour,-transpose(Psi/max(abs(Psi))),r,z_m,/overplot,nlevels=20,color=0,thick=2.0

window,31 & !P.multi=[0,2,1] 
contour,smooth(Ti_2D,[num/16,2])/tmax,yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,1],/isotropic,ytitle="Z [m]",xtitle="R [m]",title="Ti (Local)",xr=[p[0],max(p)],color=0,charsize=1.0,CHARTHICK=2.0
;for i=0,n_elements(z)-1 do oplot,p,fltarr(n_elements(p))+Z[i],psym=8,color=250
contour,transpose(Psi/max(abs(Psi))),r,z_m,/overplot,nlevels=20,color=0,thick=2.0
contour,-transpose(Psi/max(abs(Psi))),r,z_m,/overplot,nlevels=20,color=0,thick=2.0
loadct,5 & contour,emission/max(emission*0.75),yy,z,/fill,nlevels=256,xst=1,yst=1,zst=1,zr=[0,1],/isotropic,xtitle="R [m]",ytitle="Z [m]",title="Emission (Local)",xr=[p[3],max(p)],color=0,charsize=1.0,CHARTHICK=2.0
loadct,39;& for i=0,n_elements(z)-1 do oplot,fltarr(n_elements(p))+Z[i],p,psym=8,color=150
contour,transpose(Psi/max(abs(Psi))),r,z_m,/overplot,nlevels=20,color=0,thick=2.0
contour,-transpose(Psi/max(abs(Psi))),r,z_m,/overplot,nlevels=20,color=0,thick=2.0
stop

shot=210422097
time=13
data=mag_haru(dateshot=shot)
!P.background=255;16777215L
!P.color=0
r=data.c[15:86] & z_m=data.b[9:87] & Psi=data.A[9:87,15:86,time] 
restore,"\\192.168.1.111\experiment\results\Doppler\Andor\320CH\20210422\shot97.sav"
tmax=20
num=250
p=read_Ascii("C:\Users\haruaki_tanaka\Documents\ion-doppler\r.txt") & p=p.field1*1.e-3
window,31,xsize=400,ysize=400 & !P.multi=0 ;,[float(num-41)/16,2])
loadct,39
contour,smooth(Ti_2D[37:*,0:17],[float(num-37)/16,2])/tmax,yy[37:*],z[0:17],/fill,nlevels=256,yr=[z[1]-0.01,z[16]+[0.01]],xtickv =[0.11,0.16, 0.21, 0.26],/isotropic,xticks=3,xst=1,yst=1,zst=1,zr=[0,1] $
,ytitle="Z [m]",xtitle="R [m]",title="Ti (Local)",xr=[p[3],max(p)],color=0
contour,transpose(smooth(Psi/max(abs(Psi)),5)),r,z_m,/overplot,nlevels=15,color=0,C_thick=2.0
contour,-transpose(smooth(Psi/max(abs(Psi)),5)),r,z_m,/overplot,nlevels=15,color=0,C_thick=2.0
Ti_trim = Ti_2D[41:*,1:16] & yy_trim = yy[41:*] & z_trim = z[1:16]
mx = Max(Ti_trim, location)
ind = array_indices(Ti_trim, location)
print,ind, Ti_trim[ind[0],ind[1]]
oplot,[yy_trim[ind[0]]],[z_trim[ind[1]]],psym=8, thick=2, symsize=2
for i=0,n_elements(z)-1 do oplot,p,fltarr(n_elements(p))+Z[i],psym=8,color=250;fltarr(n_elements(p))+

colorbar=[[findgen(256)],[findgen(256)]];transpose([[findgen(256)],[findgen(256)]])
x=findgen(256)/256*tmax
y=[0,1]
window,xsize=800
contour,colorbar,x,y,/fill,nlevel=256,yst=7,xst=1,zst=1,title="Ti[eV]",color=0,background=255,CHARTHICK=2.0,charsize=1.0,position=[0.1,0.1,0.9,0.15]
;cb=colorbar(target=smooth(Ti_2D,[num/16,2])/tmax,orientation=0
;colorbar =[[findgen(255)-127.)/127.],[findgen(255)-127.)/127.]]
;contour,colorbar,(findgen(255)-127.)/127.,[0,1],/fill,nlevels=32,/noerase,position=[0.1,0.85,0.9,0.9],title="Ti::"+filename,color=white,xst=1,yst=4

;図面用
Ti_m = [mean([13.9208,15.7185,12.9921,14.8409,14.2638]),mean([15.1807,15.5784,16.1286,12.7858,18.8456]),mean([20.4092,18.2873,19.2168,18.7018]),mean([24.5393,23.0986,20.8332,22.6381,28.8369]),mean([28.4685,31.7799,23.1856,19.5128,29.6242]),mean([31.3746,31.9483,26.3716,27.6756,26.4899])]
Br_m = [mean([0.0146687,0.0174026,0.0148782,0.0145523,0.015247]),mean([0.0161971,0.0165794,0.0170618,0.0150024,0.0161736]),mean([0.0200045,0.01984,0.019711,0.019711]),mean([0.0234491,0.022342,0.0225947,0.0242208,0.0245243]),mean([0.026898,0.0236132,0.0256574,0.0243157,0.0241727]),mean([0.0278776,0.0274316,0.027667,0.0260187,0.025784])]
Ti_std = [stddev([13.9208,15.7185,12.9921,14.8409,14.2638]),stddev([15.1807,15.5784,16.1286,12.7858,18.8456]),stddev([20.4092,18.2873,19.2168,18.7018]),stddev([24.5393,23.0986,20.8332,22.6381,28.8369]),stddev([28.4685,31.7799,23.1856,19.5128,29.6242]),stddev([31.3746,31.9483,26.3716,27.6756,26.4899])]
Br_std = [stddev([0.0146687,0.0174026,0.0148782,0.0145523,0.015247]),stddev([0.0161971,0.0165794,0.0170618,0.0150024,0.0161736]),stddev([0.0200045,0.01984,0.019711,0.019711]),stddev([0.0234491,0.022342,0.0225947,0.0242208,0.0245243]),stddev([0.026898,0.0236132,0.0256574,0.0243157,0.0241727]),stddev([0.0278776,0.0274316,0.027667,0.0260187,0.025784])]
window,xsize=800,ysize=600 & !p.multi=0
plot,Br_m, Ti_m,psym=2, thick=2, charsize=2, xst=1,yst=1,xr=[0.015,0.028], yr=[13, 33],xtitle="Br[T]",ytitle="Ti[eV]"
errplot, Br_m, Ti_m-Ti_std,Ti_m+Ti_std

Ti = [mean([19.7501,15.8,20.43]),mean([26.5483,15.6761]),mean([12.7735,13.3002,15.906]),mean([17.459,12.8982]),mean([16.0186,20.6209,14.44])$
  ,mean([24.0215,25.3256,20.2669]),mean([15.27,16.7,21.4451]),mean([21.76]),mean([13.4571,17.49,17.9]),mean([15.87,16.09,13.7415])]
Bratio = [mean([4.554054054,5.097102063,4.180456799]),mean([4.863343477,4.7218536]),mean([5.602239598,5.5812085,5.626833015]),mean([6.608224964,7.005097642]),mean([7.553289859,8.112647348,7.480498061])$
  ,mean([4.426741643,4.944316827,5.323098692]),mean([6.390561632,5.532326875,5.445026178]),mean([5.401871578]),mean([7.680111479,7.735750007,8.433843818]),mean([7.691254359,7.717901783,7.561639357])]
Ti_std = [stddev([19.7501,15.8,20.43]),stddev([26.5483,15.6761]),stddev([12.7735,13.3002,15.906]),stddev([17.459,12.8982]),stddev([16.0186,20.6209,14.44])$
  ,stddev([24.0215,25.3256,20.2669]),stddev([15.27,16.7,21.4451]),stddev([21.76]),stddev([13.4571,17.49,17.9]),stddev([15.87,16.09,13.7415])]
Bratio_std_1 = [stddev([4.554054054,5.097102063,4.180456799]),stddev([4.863343477,4.7218536]),stddev([5.602239598,5.5812085,5.626833015]),stddev([6.608224964,7.005097642]),stddev([7.553289859,8.112647348,7.480498061])$
  ,stddev([4.426741643,4.944316827,5.323098692]),stddev([6.390561632,5.532326875,5.445026178]),stddev([5.401871578]),stddev([7.680111479,7.735750007,8.433843818]),stddev([7.691254359,7.717901783,7.561639357])]
window,30,xsize=800,ysize=600 & !p.multi=0
plot,Bratio, Ti,psym=2, thick=2, charsize=2, xst=1,yst=1,xr=[4,9], yr=[10, 40],xtitle="Br[T]",ytitle="Ti[eV]"
errplot, Bratio, Ti-Ti_std,Ti+Ti_std

stop

END

