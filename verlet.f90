! Ledicia díaz Lago                                              
! Subrutina VERLET: devuelve posiciones, velocidades y aceleraciones en el tiempo que equilibrarán las nuevas energías.                                                          


	subroutine verlet (np,rx,ry,rz,vx,vy,vz,ax,ay,az,epot,ecin,d_potencial,dd_potencial)

    use definir_precision
    use variables_comunes 
	
	implicit none
    
	!Variables de entrada: np,rx,ry,rz,vx,vy,vz,ax,ay,az
	!Variables de salida:epot,ecin,d_potencial,dd_potencial,rx,ry,rz,vx,vy,vz,ax,ay,az
	real(kind=doblep),intent(in)::np
	real(kind=doblep),dimension(npmax),intent(inout)::rx,ry,rz,vx,vy,vz,ax,ay,az
    real(kind=doblep),intent(out):: d_potencial,dd_potencial,epot,ecin

	
	!Posiciones y velocidades para el primer intervalo de tiempo:
    !write(*,*)'hola'
    rx = rx + vx*dt + ax*dt2 
    ry = ry + vy*dt + ay*dt2
    rz = rz + vz*dt + az*dt2

    vx = vx + ax*dt_medios
    vy = vy + ay*dt_medios
    vz = vz + az*dt_medios
	
    !Para seguir calculando necesitamos los valores de las aceleraciones a través de las nuevas posiciones 
	!Para eso llamamos a la subrutina del Potencial de Lernard-Jones

    call energia_potencial(np, rx, ry, rz, ax, ay, az, epot, d_potencial, dd_potencial)
    
	!Los nuevos valores de las aceleraciones se utilizan para recalcular las velocidades y por tanto la energía cinética:

    vx = vx+ax*dt_medios
    vy = vy+ay*dt_medios
    vz = vz+az*dt_medios
    
    ecin = 0.5d00*(sum(vx*vx)+sum(vy*vy)+sum(vz*vz))

    return
    
	end subroutine verlet

    

    