Function Trigrid_interpor_for_r_lambda,z,x,y,xout,yout

;Nterms=N_elements(x)*N_elements(y)
Nterms=N_elements(x)
Xtemp=DBLarr(Nterms) & Ytemp=DBLarr(Nterms) & SigTemp=DBLarr(Nterms)

For i=0,Nterms-1 do begin
    index=i mod N_elements(Y)
    Xtemp[i]   = X[floor(i/N_elements(Y)),index]
    Ytemp[i]   = Y[index]
    SigTemp[i] = Z[floor(i/N_elements(Y)),index]
Endfor

;PLOT,Xtemp,Ytemp,psym=1,xstyle=1,ystyle=1;,charsize=2
       Triangulate,Xtemp,Ytemp,Tr,B1
       FOR i=0, N_ELEMENTS(tr)/3-1 DO BEGIN
            ; Subscripts of vertices [0,1,2,0]:
            t_index = [tr[*,i], tr[0,i]]
         ; Connect triangles:
            OPLOT, Xtemp[t_index], Ytemp[t_index],color=100
       ENDFOR
;A=Trigrid(Xtemp,Ytemp,SigTemp,Tr,xout=xout,yout=yout,Extra=B1);,/quintic)
A=Trigrid(Xtemp,Ytemp,SigTemp,Tr,xout=xout,yout=yout);,/quintic)
;contour,A,Xout,Yout,nlevels=32,/fill,xst=1,yst=1,zst=1
        FOR i=0, N_ELEMENTS(tr)/3-1 DO BEGIN
            ; Subscripts of vertices [0,1,2,0]:
            t_index = [tr[*,i], tr[0,i]]
            ; Connect triangles:
            OPLOT, Xtemp[t_index], Ytemp[t_index],color=100
        ENDFOR

result={Z:A,X:Xout,Y:Yout}
return, result
END
