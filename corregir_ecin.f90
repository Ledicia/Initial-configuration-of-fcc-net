subroutine corr_ecin(np,vx,vy,vz,ecin,epot,etotq)
	
	use definir_precision
    use variables_comunes

    implicit none

!	Definimos las variables utilizadas en la subrutina indicando su formato: entrada ('in') o salida ('out'):
	integer(kind=entero), intent(in):: np
	real(kind=doblep), dimension(npmax), intent(inout):: vx,vy,vz
    real(kind=doblep), intent(inout):: ecin
	real(kind=doblep), intent(in):: etotq,epot

!	Definimos:
!	ecinq: energía cinética deseada
!	corr_v: factor para corregir las velocidades
!	px,py,pz: momentos lineales
	real(kind=doblep):: ecinq,corr_v,px,py,pz 
    !npd: np en doble precisión:
    real(kind=doblep)::npd 
    npd=dble(np)
 
!	Nos aseguramos de que el momento total sea nulo:
	px=sum(vx)
	py=sum(vy)
	pz=sum(vz)
    
!	Para eso, corregimos las velocidades:
	px=px/npd
	py=py/npd
	pz=pz/npd
    
	!Corrijo las velocidades restando el momento lineal total por particula:
    vx = vx-px
    vy = vy-py
    vz = vz-pz
    
	ecinq = etotq-epot

!	Nos aseguramos de que la energía cinética no es negativa:
	if (ecinq <=0.d00) then  
    	write (*,*) "Error: energia cinetica negativa"
        stop
   	end if

!	Calculamos el factor de corrección:
	corr_v=dsqrt(ecinq/ecin)

!	Calculamos las velocidades corregidas:
	vx=corr_v*vx
	vy=corr_v*vy
	vz=corr_v*vz 

!	Comprobamos que el momento sea nulo:
	px=sum(vx)
	py=sum(vy)
	pz=sum(vz)
 
	write(*,*) 'Momento despues de la correccion:',px,py,pz

!	Calculamos la energía cinética:
	ecin=0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))
	write(*,*) 'ecin', ecin
    return

end subroutine corr_ecin
     
