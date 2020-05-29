!Ledicia Díaz Lago
		function random(idum)
        implicit none
        
        integer, parameter :: entero=SELECTED_INT_KIND(9)
        integer, parameter :: doblep=SELECTED_REAL_KIND(15,307)

        !random es real*8. idum integer*4
        !devuelve un numero alatorio entre (0,1). 
		!Si el parametro idum es negativo, o es la primera vez que se la llama inicializa el generador en base a idum.

        integer (kind=entero) :: idum
        real (kind=doblep), parameter :: mbig=4.d+06,mseed=1618033.d00
        real (kind=doblep), parameter :: mz=0.d00,fac=1.d00/mbig
        integer (kind=entero) :: i,iff,ii,inext,inextp,k
        real (kind=doblep) :: random,mj,mk,ma(55)
        save iff,inext,inextp,ma                                          !Guarda esto para la proxima llamada
        data iff/0/                                                       !Inicializa a cero cuando compila.
                                                                          !Solo asigna un valor a iff la primera c¡vez que se llama a la funcion
                                                                          
		!Comprueba se es la primeira vez que se llama a la función.
        !Si idum es negativo ou cero vuelve a inicializar la función; sino sigue con los datos de antes.
        if (idum<=0.or.iff==0) then
           iff=1
           mj=abs(mseed-abs(idum))
           mj=mod(mj,mbig)
           ma(55)=mj
           mk=1
           do i=1,54
              ii=mod(21*i,55)
              ma(ii)=mk
              mk=mj-mk
              if (mk.lt.mz) mk=mk+mbig
              mj=ma(ii)
           end do
           do k=1,4
              do i=1,55
                 ma(i)=ma(i)-ma(1+mod(i+30,55))
                 if (ma(i).lt.mz) ma(i)=ma(i)+mbig
              end do
           end do
           inext=0
           inextp=31
           idum=1
        end if
        inext=inext+1
        if (inext==56) inext=1
        inextp=inextp+1
        if (inextp==56) inextp=1
        mj=ma(inext)-ma(inextp)
        if(mj.lt.mz) mj=mj+mbig
        ma(inext)=mj
        random=mj*fac
        return
        end function random
    

