!******************************************************************************************************
! Ledicia Díaz Lago                                                                              
! Simulación física de materiales 2018/2019
! Equilibración de la configuración inicial en el tiempo para un diferencial de tiempo dt=1.d-04 con kpasos=5.d05 grabados cada 100.
! La energía total es -550, la densidad ro=0.5 y suponemos 500 partículas.
! Los datos necesarios se leen de un fichero de resultados_simulación.
! Trabajo en unidades reducidas m=1 (masa).                                   
!******************************************************************************************************

   
	program equilibracion

    use variables_comunes
    use definir_precision

    implicit none 


!	Variables necesarias para la equilibración:
!	kpaso: entero que cuenta cada cuanto grabamos el número de iteraciones en el tiempo (intervalo de pasos).
!	ktotal: entero que cuenta el numero de iteraciones en el tiempo (pasos totales).
!	xpaso: kpaso declarado como doble precisión.
!	tiempo: acumulador del tiempo que transcurre.
!	n, ki, kcuenta: enteros que se emplean en los bucles.

!   Variables que definen las 3N posiciones,velocidades y aceleraciones de las partículas.
	real(kind=doblep),dimension(npmax)::rx,ry,rz,vx,vy,vz,ax,ay,az 
!	kpaso: entero que cuenta cada cuanto grabamos el número de iteraciones en el tiempo.
!	ktotal: entero que cuenta el numero de iteraciones en el tiempo.    
    integer(kind=entero):: kpaso, ktotal
    !Numero de particulas, y número de celdas por arista.
	integer (kind=entero):: np!, ka
    !npd: np declarado como doble precisión.
    !real(kind=doblep):: npd
	!xpaso: kpaso declarado en doble precisión.
    real(kind=doblep):: xpaso
    !tiempo: acumulador del tiempo que transcurre.
    real(kind=doblep):: tiempo
    !ki, kcuenta:enteros que se emplean en los bucles.
    integer(kind=entero):: i,ki,kcuenta
    !Energias cinetica y potencial.
	real(kind=doblep)::ecin,epot                                         
    !Energia total del sistema, y energía total que queremos.
	real (kind=doblep) :: etot,etotq                        
    !Derivadas primera y segunda del potencial de Lennard-Jones.                             
    real(kind=doblep):: d_potencial, dd_potencial   
    !Archivos que se emplean para guardar los resultados.                    
	character (len=55) :: fname,gname1,gname2 

	!Recopilamos los datos de la configuración inicial abriendo el fichero resultados_simulacion.dat:
	!Leo las variables desde el archivo de texto ASCII.
                
		write (*,*) 'Fichero de datos simulacion'    !Introducimos el nombre del fichero que queremos leer.
		read (*,9000) fname
9000	format (a25)
        open (10,file=fname)                         !Abre el fichero y lee los datos.
        	read (10,*)                              !Salto la primera línea (no la quiero leer).
        	read (10,*) np,pl,pli,rc,rc2             !Leo la línea 2.
            read (10,*)                              !Salto la tercera línea (no la quiero leer).
    		read (10,*) vol,ro,dt,kpaso,ktotal       !Leo la línea 4.
            read (10,*)                              !Salto la quinta línea (no la quiero leer).
			read (10,*) etot,ecin,epot               !Leo la línea 6.
			read (10,9000) gname1                    !Fichero donde están grabados los resultados de la configuración inicial.
			read (10,9000) gname2                    !Fichero donde están grabadas posiciones, velocidades y aceleraciones de la configurcion inicial.
		close (10)
		open (20,file='resultados_cinematica.dat',form='unformatted')
        	read (20) rx,ry,rz,vx,vy,vz,ax,ay,az
		close (20)
        
	!La energía total que queremos es la que leemos del archivo datos_simulacion.txt:
	etotq = etot
    
	!Defino diferentes variables que necesito: 
	!ka = k                   !Número de celdas por arista.
	np = npmax               !Número de partículas.
	!npd = dble(np)			 !Número de partículas en doble precisión.
    xpaso = dble(kpaso)		 !Kpaso en doble precisión.
    dt_medios = dt/ 2.d00    !Mitad del diferencial del tiempo.
    dt2 = dt_medios*dt       !Mitad del diferencial del tiempo al cuadrado.
	!pli = 1.d00/pl
    !ro = npd/vol
    
	!El sistema va a evolucionar hasta el equilibrio manteniendo constante el valor de la energía total prefijada.
	!El hecho de que las partículas estén perfectamente colocadas (red fcc) hace que la energía al principio de la simulación pege un salto,
    !de forma que, de seguir por este camino, se alcanzaría el equilibrio a una energía distinta de la deseada. 
    !Para corregir este cambio de energía hacemos evolucionar el sistema en un pequeñointervalo temporal (1000 pasos) para corregir la energía cinética 

	!Hacemos evolucionar el sistema 1000 pasos:		
	do i=1,1000
		call verlet (np,rx,ry,rz,vx,vy,vz,ax,ay,az,epot,ecin,d_potencial,dd_potencial)
	end do
	etot = ecin + epot
    write(*,*) etot
	write(*,*) ecin
    write(*,*) epot
    
    call corr_ecin(np,vx,vy,vz,ecin,epot,etotq)
	etot = ecin + epot
	write(*,*) 'Energia total despues de la correccion:',etot


	!A medida que avanza el tiempo se necesita una configuración que equilibre los valores de la energía en el tiempo. 
	!Aparecen modificaciones en posiciones,velocidades y aceleraciones, mediante la subrutina VERLET.

    kcuenta = 0
    do ki = 0,ktotal
      call verlet(np,rx,ry,rz,vx,vy,vz,ax,ay,az,epot,ecin,d_potencial,dd_potencial) 
      etot = epot + ecin

      !Esto podríamos meterlo dentro de la subrutina. 
      if(mod(ki,kpaso)==0) then !graba cada kpasos
        kcuenta = kcuenta + 1
        tiempo = dble((kcuenta-1))*xpaso*dt
        
        !Estoy utilizando el fichero resultados_cinematica.dat
		!Almacenamos los datos para la representación.        
        open(10,file='resultados_energias.dat',access='append')
            write(10,9008) tiempo,etot,ecin,epot
        close(10)   
       end if
    end do

9008 format(1pe19.12,3(2x,e19.12))

    
	!Hecho esto grabamos la configuración final de todos los datos: 
    open(unit=1001, file='resultados_simulacion.txt', action='write')
        write(1001,9001) 'np', 'pl', 'pli', 'rc', 'rc2' 
        write(1001,9002) np,pl,pli,rc,rc2        
        write(1001,9003) 'vol', 'ro', 'dt', 'kpaso', 'ktotal' 
        write(1001,9004) vol,ro,dt,kpaso,ktotal   
        write(1001,9005) 'etot', 'ecin', 'epot'  
        write(1001,9006) etot,ecin,epot        
        write(1001,9000) fname
        write(1001,9000) gname1
        write(1001,9000) gname2
    close(1001)
    
    open(30,file='equilibracion_cinematica.dat',form='unformatted')
        write(30) rx,ry,rz,vx,vy,vz,ax,ay,az
    close(30)

    write(*,*) 'Equilibracion terminada'

!   Formatos utilizados:
9001 format(a3, a14,4(2x, a16))                             !Formato de la línea de texto.
9002 format(i4,2x,1pe19.12,3(2x,e19.12))                    !Formato para los valores numéricos en coordenadas reducidas.
9003 format(a3, 2(2x, a22), 2(2x, a12))  
9004 format(1pe19.12,2x,e19.12,2x,e19.12,2x,i4,2x,i6)   
9005 format(a10,2(2x, a22))  
9006 format(1pe19.12,2x,e19.12,2x,e19.12)
 



	end program 