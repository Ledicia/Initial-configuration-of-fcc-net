!Ledicia Díaz Lago

	        module variables_comunes
            use definir_precision
            implicit none

            !Variables comunes que van ser útiles en todo el programa.
            !k: número de celdas por arista.
            !np: número total de partículas (500).
            integer (kind=entero), parameter :: k=5, npmax=500
            
            !Parámetros para la simulación (se trabaja en coordenadas reducidas masa=1).
            !pl: tamaño de la caja.
            !vol: volume de la caja.
            !rc: radio de máximo alcance del potencial.
            !a: tamaño de cada celda de la red (parametro de corte).
            !etotq: energía total que quiero que tenga el sistema (en unidades reducidas).
            !ro: densidad de partículas (en unidades reducidas).
            real (kind=doblep) :: pl, pli, vol, rc, rc2, ro
            !Correccion_potencial: correción a la energía potencial necesaria por usar un radio de corte finito.
            !Parámetros de corrección de derivadas:
			real (kind=doblep) :: correccion_potencial, correccion_dpot, correccion_ddpot, dt, dt_medios, dt2
            
            !Parámetros del potencial de Lennard-Jones, en Julios y metros.
            !real (kind=doblep), parameter :: eps, sigma

            !Se define pi como un parámetro constante.
            real (kind=doblep), parameter :: pi = 3.14159265359d00

        end module variables_comunes


