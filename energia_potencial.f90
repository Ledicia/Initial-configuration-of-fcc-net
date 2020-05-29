!Ledicia Díaz Lago
        subroutine energia_potencial(np, rx, ry, rz, ax, ay, az, epot, d_potencial, dd_potencial)
        !Esta subrutina toma valores de posiciones, velocidades, parámetros da simulación e energía total.
        !Calcula la energía potencial del sistema, y las fuerzas sobre cada partícula (aceleraciones, ya que estamos en coordenadas reducidas)
        !Tambien calcula la suma de todas las derivadas primer y segundo orden de la energía potencial.

	
		use definir_precision
        use variables_comunes
        implicit none
                    
		!Variable de entrada con el número de partículas.
        integer (kind=entero), intent(in):: np
        !Variables de entrada con las posicionesy velocidades de las partículas.
        real (kind=doblep), dimension(npmax), intent(in) :: rx, ry, rz
        !Variables de salida con las aceleraciones de las partículas, y la energía potencial.
        real (kind=doblep), dimension(npmax), intent(out) :: ax, ay, az
        real (kind=doblep), intent(out) :: epot
        !Suma de las derivadas de primeir y segundo orden de la energía potencial.
        real (kind=doblep), intent(out) :: d_potencial, dd_potencial
        !Variables para os bucles.
        integer (kind=entero) :: i, j
        !Variables para guardar la posición da la partícula actual.
        real (kind=doblep) :: xi, yi, zi
        !Variables para as distancias entre partículas. pli es la inversa de pl en coordenadas reducidas.
        real (kind=doblep) :: delta_x, delta_y, delta_z
        !Variables para guardar las distancias, y distintas potencias da inversa da distancia.
        !El número del nombre indica la potencia a la que está elevada la inversa.
        real (kind=doblep) :: rij_cuadrado, rij_inv1, rij_inv2, rij_inv6, rij_inv12
        !Variables para almacenar la derivada del potencial entre dos partículas, y las componentes de la fuerza.
        real (kind=doblep) :: d_potencial_ij, fuerza_ij_x, fuerza_ij_y, fuerza_ij_z
        !Variable auxiliar para la corrección de la energía potencial.
   	    real(kind=doblep):: factor1
        
		!Inicializo a cero la energía potencial y sus derivadas.
        !write(*,*) 'entre', np; pause
        epot = 0.d00
        d_potencial = 0.d00
        dd_potencial = 0.d00
		
		pli = 1.d00 / pl
        !Inicializo a cero todas las aceleraciones (o fuerzas).
        ax=0.d00
    	ay=0.d00
    	az=0.d00
        
        !Calculo la energía potencial total, las derivadas y la fuerza sobre cada partícula.
  		do i=1, np-1
        !Posicion de la particula i
        	xi = rx(i)
            yi = ry(i)
            zi = rz(i)
            do j=i+1, np
            !delta_k = rij: distancia entre las particulas en el eje k=x,y,z.
            delta_x = xi - rx(j)
            delta_y = yi - ry(j)
            delta_z = zi - rz(j)
            
			!Aplico las condiciones de frontera periódicas:
            delta_x = delta_x - pl * dnint(delta_x * pli)   !la funcion DNINT reondea al entero más proximo con doble precisión
            delta_y = delta_y - pl * dnint(delta_y * pli)	!Restamos el lado de la caja multiplicado por el número de cajas: 
            delta_z = delta_z - pl * dnint(delta_z * pli)	!-pl*(rij*pli=rij*1/pl)
                            
            !Calculo la distancia al cuarado y las distintas potencias da inversas necesarias.
            rij_cuadrado = delta_x*delta_x + delta_y*delta_y + delta_z*delta_z
            if (rij_cuadrado < rc2) then
            	!Potencial de Lennard_Jones:  V(r)=4.d00*eps*[(sigma/rij)^12-(sigma/rij)^6]
            	!Se necesitan las inversas de las distancias entre particulas.
            	rij_inv1 = dsqrt(1.d00 / rij_cuadrado)
            	rij_inv2 = 1.d00 / rij_cuadrado
            	rij_inv6 = rij_inv2 * rij_inv2 * rij_inv2
            	rij_inv12 = rij_inv6 * rij_inv6
                                
            	!Calculo el potencial entre las dos particulas i, j, y se lo sumo al potencial total.
            	!Asi al terminar el bucle tendremos las constribuciones de todas las particulas.
            	epot = epot + (rij_inv12 - rij_inv6)
                                
            	!Calculo la derivada del potencial entre las dos partículas, y lo utilizo para calcular las fuerzas en cada una.
            	!No metemos los factores comunes porque los introducimos en los cálculos finales por comodidad.
            	!De la regla de la cadena obtenemos el factor que multiplica *delta_x*rij_inv1
            	d_potencial_ij = (-2.d00 * rij_inv12 + rij_inv6) * rij_inv1
            	fuerza_ij_x = -d_potencial_ij * delta_x * rij_inv1     
            	fuerza_ij_y = -d_potencial_ij * delta_y * rij_inv1
            	fuerza_ij_z = -d_potencial_ij * delta_z * rij_inv1
                                
            	!Sumo esta fuerza para las dos partículas. 
            	!Como estamos trabajando con coordenadas reducidas (masa=1) son equivalentes a las aceleraciones.
            	ax(i) = ax(i) + fuerza_ij_x
            	ay(i) = ay(i) + fuerza_ij_y
            	az(i) = az(i) + fuerza_ij_z
            	!Para la partícula j, la fuerza va en el sentido contrario.
            	ax(j) = ax(j) - fuerza_ij_x
            	ay(j) = ay(j) - fuerza_ij_y
            	az(j) = az(j) - fuerza_ij_z
                                
            	!Sumo la derivada de la energía potencial al total.
            	!Asi al terminar el bucle tendremos las constribuciones de todas las particulas.
            	d_potencial = d_potencial + d_potencial_ij
            	dd_potencial = dd_potencial + (26.d00 * rij_inv12 - 7.d00 * rij_inv6) * rij_inv2
            end if
        end do
     end do
    !Multiplico los valores calculados por factores que antes se eliminaron por simplicidad.
    epot = 4.d00 * epot                               !V(r) lleva un factor 4
    ax = 24.d00 * ax                                  !V´(r) proporcional a 6*4 (en aux1 se ha sacado un factor comun 6).
    ay = 24.d00 * ay
    az = 24.d00 * az
    d_potencial = 24.d00 * d_potencial
    dd_potencial = 24.d00 * dd_potencial

    !Corregimos las energías, ya que, a la hora de calcular el potencial, se ha truncado la funcion para distancias mayores que el radio de corte.
	!Esto conlleva a correcciones en los valores de la energía potencial y de las derivadas, considerando la función de distribución radial igual a 1.

	factor1 = 8.d00*pi*np**2/(vol*rc*rc*rc)             
    correccion_potencial = factor1*(1.d00/(3.d00*rc**6)-1.d00)/3.d00      
    correccion_dpot = 2.d00*factor1*(-2.d00/(3.d00*rc**6)+1.d00)
    correccion_ddpot = 2.d00*factor1*(26.d00/(3.d00*rc**6)-7.d00)
    
	!Valores finales de la energía y de las derivadas:
    epot = epot+correccion_potencial 
    d_potencial = d_potencial+correccion_dpot
    dd_potencial = dd_potencial+correccion_ddpot
    
    d_potencial = d_potencial/(3.d00*vol)
    dd_potencial = (dd_potencial-2.d00*d_potencial)/(9.d00*vol*vol)
    return
    
end subroutine energia_potencial
		
            
    	