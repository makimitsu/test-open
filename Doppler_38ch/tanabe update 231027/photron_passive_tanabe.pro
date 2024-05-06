PRO Photron_passive_tanabe
a= findgen(17) * (!PI*2.0/16.0)
USERSYM, cos(A), sin(A), /FILL
device,decomposed=0 & loadct,39
;restore,"\\OBI-WAN\htanabe_backup\fast_camera\UVI\0530\shot09\shot09.tif_accumulate.sav"
;restore,"\\OBI-WAN\htanabe_backup\fast_camera\UVI\0719\shot02\shot02.tif_accumulate.sav"
;restore,dialog_pickfile(path="\\OBI-WAN\htanabe_backup\fast_camera\UVi_calibration\",filter="*.sav")
;fil=dialog_pickfile(path="\\OBI-WAN\htanabe_backup\fast_camera\UVi\",filter="*.tif")
;fil=dialog_pickfile(path="\\192.168.1.111\experiment\results\Doppler\Photron\",filter="*.tif")
fil=dialog_pickfile(path="\\Mac\Home\Downloads\Doppler\Photron\230914",filter="*.tif")
;fil="\\Mac\Home\Downloads\Doppler\Photron\230913\H_dial48635_shot20_300kfps_3us_gain20\H_dial48635_shot20_300kfps_3us_gain20.tif"
D=read_tiff(fil,image_index=0)
p=read_Ascii("\\Mac\Home\Downloads\Doppler\Photron\calibration_summary\2022_Funato_38CH\test2\radius.txt") & p=p.field1*1.e-3 & p=p
relative=read_Ascii("\\Mac\Home\Downloads\Doppler\Photron\calibration_summary\2022_Funato_38CH\test2\relative.txt") & relative=relative.field1
smile= read_Ascii("\\Mac\Home\Downloads\Doppler\Photron\calibration_summary\2022_Funato_38CH\test2\smile.txt") & smile=smile.field1
center= read_Ascii("\\Mac\Home\Downloads\Doppler\Photron\calibration_summary\2022_Funato_38CH\test2\center.txt") & center=center.field1
width=5
Instrument=read_ascii("\\Mac\Home\Downloads\Doppler\Photron\calibration_summary\2022_Funato_38CH\test2\instrument.txt") & Instrument=Instrument.field1

center=center-((1024.)-n_elements(d[*,0]))/2.

Tmax=50.
nframe=249;00;250;200;300
lambda0=486.133
offset=0.;21.;0.15;-0.325
;resolution=0.011941488;-0.00000003104171*lambda0^2 + 0.00002333564*lambda0 + 0.007991794
resolution=0.01393172
Data=fltarr([n_elements(d[*,0]),n_elements(d[0,*]),nframe])

Instrument=Instrument*resolution
M=1.
Ti_instru=1.69e8*M*(2.*instrument*sqrt(2.*alog(2.))/lambda0)^2
rate=300.e3
time=findgen(nframe)/rate*1.e6
lambda=-(findgen(n_elements(D[0,*]))-n_elements(D[0,*])/2.)*resolution+lambda0;-offset ;波長ひっくり返ってるので反映させる
spectra=fltarr([36,n_elements(d[0,*]),nframe])
delay=0.;80.
tend=nframe-1
;stop

V_pixel=n_elements(data[0,*,0])
smile=smile-512.+V_pixel/2.+offset
bg=max(where(time lt 400.))

for i=0,nframe-1 do begin
;for i=0,0 do begin
Data[*,*,i]=read_tiff(fil,image_index=i)
  for j=0,35 do begin
    subimage=Data(center[j]-width:center[j]+width,0:*,i:i)
    spectra[j,*,i]=total(subimage,1,/double)*relative[j]
   spectra[j,*,i]=spectra[j,*,i]-total(reform(spectra[j:j,0:*,0:bg]),2)/(bg+1.); バックグラウンド信号除去
;    subimage=D(center[j]-width:center[j]+width,0:*)
;    spectra[j,*,i]=total(subimage,1)
  endfor
endfor

!P.multi=[0,9,4] & window,0,xsize=1800,ysize=800
for i=0,35 do contour,spectra[i,*,*],lambda,time,/fill,nlevels=64,xst=1,yst=1,zst=1,yr=[450,600],xr=[lambda0-0.5,lambda0+0.5],title="spectrum_CH::"+strcompress(i+1)
;spectra=smooth(spectra,[2,2,1])
Ti=fltarr([36,nframe])
Ti_max=fltarr([36,nframe])
Ti_min=fltarr([36,nframe])
flow=fltarr([36,nframe])
flow_max=fltarr([36,nframe])
flow_min=fltarr([36,nframe])

lambda2=fltarr([n_elements(lambda),36])
for i=0,35 do lambda2[*,i]=-(findgen(n_elements(D[0,*]))-n_elements(D[0,*])/2.-smile[i])*resolution+lambda0;-offset*resolution)
for i=0,35 do contour,spectra[i,*,*],lambda2[*,i],time,/fill,nlevels=64,xst=1,yst=1,zst=1,yr=[400,600],title="spectrum_CH::"+strcompress(i+1)

for i=delay,tend do begin
;for i=0,0 do begin
  for j=0,35 do begin
    fitting=gaussfit(reform(lambda2[*,j]),spectra[j,*,i],coeff,nterms=3,sigma=sigma)
    Ti[j,i]=1.69e8*M*(2.*coeff[2]*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru[j]
    Ti_max[j,i]=1.69e8*M*(2.*(coeff[2]+3.*sigma[2])*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru[j]
    Ti_min[j,i]=1.69e8*M*(2.*(coeff[2]-3.*sigma[2])*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru[j]
;    spectra[j,*,i]=spectra[j,*,i]-coeff[3]
;    flow[j,i]=3.0e8/lambda0*(coeff[1]-(smile[j]*resolution+lambda0-n_elements(D[0,*])/2.*resolution))/1000.
    flow[j,i]=3.0e8/lambda0*(coeff[1]-lambda0)/1000.;*resolution
;    plot,findgen(n_elements(lambda)),spectra[j,*,0],xst=1,yst=1
;    fitting=gaussfit(findgen(n_elements(lambda)),spectra[j,*,i],coeff,nterms=6)
;    oplot,findgen(n_elements(lambda)),fitting,color=250
;    print,coeff
  endfor
endfor
checker=(Ti ge abs(Ti_max-Ti_min)) *(Ti lt 300) 
Ti=Ti*checker
Ti_max=Ti_max*checker
Ti_min=Ti_min*checker

loadct,39
window,2,ysize=700 & !P.multi=0
contour,transpose(Ti),time,p,/fill,nlevels=256,xst=1,yst=1,zst=1,xr=[450,550],zr=[0,tmax],title="Ti [eV]",charsize=2
for i=0,n_elements(time)-1 do oplot,time[i]*float(p gt 0.),p,psym=8
!P.background=16777215L & !P.color=0

window,3,ysize=700 & !P.multi=0
contour,transpose(flow),time,p,/fill,nlevels=256,xst=1,yst=1,zst=1,xr=[450,550],title="flow [km/s]",charsize=2,zr=[-50.,50.]
!P.background=16777215L & !P.color=0

window,4
contour,transpose(total(spectra,2,/double)),time,p,/fill,nlevels=256,xst=1,yst=1,zst=1,xr=[450,550],title="emission [a.u.]",charsize=2,zr=[0.,max(total(spectra,2,/double))]
!P.background=16777215L & !P.color=0
for i=0,n_elements(p)-1 do oplot,time,p[i]*(fltarr(n_elements(time))+1.),psym=8,color=50

;!P.background=16777215L & !P.color=0 & Tmax=100
window,5
plot,time,total(Ti,1)/36.,xst=1,yst=1,psym=-8,yr=[0,tmax],charsize=1.5,xtitle="time [us]",ytitle="Ti [eV]",title=fil,xr=[400,600]
ERRPLOT,time,total(Ti_min,1)/36.,total(Ti_max,1)/36.
;window,5,xsize=500,ysize=500
;plot,time-400.,total(Ti,1)/36.,xst=1,yst=1,psym=-8,yr=[0,tmax],charsize=1.5,xtitle="time [us]",ytitle="Ti [eV]",title=fil,xr=[50,100]
;ERRPLOT,time-400.,total(Ti_min,1)/36.,total(Ti_max,1)/36.

!P.background=16777215L & !P.color=0
window,6 & !P.multi=0
plot,time,total(flow,1)/36.,xst=1,yst=1,psym=-8,charsize=1.5,xtitle="time [us]",ytitle="Vi [km/s]",title=fil,xr=[400,600],yr=[-50,50]
ERRPLOT,time,total(flow_min,1)/36.,total(flow_max,1)/36.

!P.color=16777215L & !P.background=0
;stop

;Start_Abel_inversion
num=250.
emission=fltarr([num+1,n_elements(time)])
input=total(spectra,2,/nan,/double)
edge=0.35
yy=findgen(num+1)*(edge-min(p))/num+min(p)
dy=(edge-min(p))/num
interp2=trigrid_interpor3([input,transpose(fltarr(249))],[p,edge],time,yy,time)

for i=delay,tend  do begin
;  interp=spline([p,edge],[input[*,i],0.],yy)
interp=reform(interp2.z[*,i])
  dIdt=deriv(yy,interp)
  for j=0,num-1 do begin
;      emission[j,i]=total(-1./!PI*(dIdt/sqrt(yy^2-yy[j]^2)*dy),/nan,/double)
      emission[j,i]=total(-1./!PI*(dIdt[j+1:*]/sqrt(yy[j+1:*]^2-yy[j]^2)*dy),/nan,/double)
 endfor
endfor

loadct,39
;contour,smooth(transpose(emission),[2,num/16.]),time-400.,yy,xr=[0,200],/fill,nlevels=256,zst=1,xst=1,yst=1,zr=[0,max(emission)],charsize=2;,color=100,
contour,transpose(emission),time-400.,yy,xr=[0,200],/fill,nlevels=256,zst=1,xst=1,yst=1,zr=[0,max(emission)],charsize=2;,color=100,
;loadct,39
for i=0,n_elements(p)-1 do oplot,time,p[i]*(fltarr(n_elements(time))+1.),psym=8,color=50
;save,filename=fil+"_EM.sav",emission,time,p,yy
;stop

Ti_r=emission*0.
t_start=min(where(time ge 450.))
t_end=max(where(time le 550.))
lambda_out=lambda[3:n_elements(lambda)-5]
spectra_interp=fltarr([n_elements(lambda_out),n_elements(yy),t_end-t_start+1])
time2=time[t_start:t_end]
for i=0,t_end-t_start do begin
result=trigrid_interpor_for_r_lambda2([[transpose(reform(spectra[*,*,i+t_start]))],[fltarr(n_elements(lambda))]],[[lambda2],[lambda]],[p,edge],lambda_out,yy)
spectra_interp[*,*,i]=result.z
endfor
Local_spectra=spectra_interp*0.
Ti_local=transpose(total(spectra_interp,1)*0.)
Ti_local_max=Ti_local
Ti_local_min=Ti_local
peak=Ti_local
emission2=Ti_local
D_shift=Ti_local
Ti_instru2=total(ti_instru)/36.

spectra_interp=smooth(spectra_interp,[5,num/18.,1])
;stop
for i=0,t_end-t_start do begin
  for j=0,n_elements(lambda_out)-1 do begin
    dIdt=deriv(yy,reform(spectra_interp[j,*,i]))
    for k=1,num-1 do begin
      Local_spectra[j,k,i]=total(-1./!PI*(dIdt[k+1:*]/sqrt(yy[k+1:*]^2-yy[k]^2)*dy),/nan,/double)
    endfor
  endfor
  emission2[i,*]=total(Local_spectra[*,*,i],1,/nan,/double)
  for j=0,num do begin
    fit=gaussfit(lambda_out,reform(Local_spectra[*,j,i]),coeff,nterms=4,sigma=sigma)
    Ti_local[i,j]=    1.69e8*M*(2.* coeff[2]              *sqrt(2.*alog(2.))/lambda0)^2-Ti_instru2
    Ti_local_max[i,j]=1.69e8*M*(2.*(coeff[2]+abs(sigma[2]))*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru2
    Ti_local_min[i,j]=1.69e8*M*(2.*(coeff[2]-abs(sigma[2]))*sqrt(2.*alog(2.))/lambda0)^2-Ti_instru2
    peak[i,j]=coeff[0]
    D_shift[i,j]=coeff[1]
  endfor
endfor
checker=(peak gt 0) * (D_shift lt lambda0+0.2) * (D_shift gt lambda0-0.2) * (Ti_local_max-Ti_local_min lt Ti_local*2.)*(Ti_local gt 0)  
Ti_local=Ti_local*checker
Ti_local_max= Ti_local_max*checker
Ti_local_min= Ti_local_min*checker
checker2=where((finite(Ti_local) ne 1))
for i=0,n_elements(checker2)-1 do begin
  Ti_local[checker2[i]]=0.
  Ti_local_max[checker2[i]]=0.
  Ti_local_min[checker2[i]]=0.
endfor
checker=Ti_local gt 0

for kkk=0,100 do begin
  for i=0,t_end-t_start do begin
    for j=1,num-1 do begin
      if checker[i,j] eq 0 then begin
        Ti_local[i,j]=(Ti_local[i,j+1]+Ti_local[i,j-1])/2.
        Ti_local_max[i,j]=(Ti_local_max[i,j+1]+Ti_local_max[i,j-1])/2.
        Ti_local_min[i,j]=(Ti_local_min[i,j+1]+Ti_local_min[i,j-1])/2.
      endif
    endfor
  endfor
endfor
window,7 & loadct,39
contour,emission2,time2,yy,/fill,nlevels=256,zr=[0,max(emission2)]
for i=0,n_elements(p)-1 do oplot,time,p[i]*(fltarr(n_elements(time))+1.),psym=8,color=50

window,8
loadct,39 & contour,smooth(ti_local,[num/18,1]),time2,yy,/fill,nlevels=256,zr=[-0.1,25],color=255
;loadct,39 & contour,smooth(ti_local,[num/18,1]),time2,yy,/fill,nlevels=256,zr=[-0.1,15],color=255,/nodata
;loadct,25 & contour,smooth(ti_local,[num/18,1]),time2,yy,/fill,nlevels=256,zr=[-0.1,15],/overplot
;for i=0,n_elements(p)-1 do oplot,time,p[i]*(fltarr(n_elements(time))+1.),psym=8,color=50
stop
END
