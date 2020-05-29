!Ledicia D�az Lago

	        module variables_comunes
            use definir_precision
            implicit none

            !Variables comunes que van ser �tiles en todo el programa.
            !k: n�mero de celdas por arista.
            !np: n�mero total de part�culas (500).
            integer (kind=entero), parameter :: k=5, npmax=500
            
            !Par�metros para la simulaci�n (se trabaja en coordenadas reducidas masa=1).
            !pl: tama�o de la caja.
            !vol: volume de la caja.
            !rc: radio de m�ximo alcance del potencial.
            !a: tama�o de cada celda de la red (parametro de corte).
            !etotq: energ�a total que quiero que tenga el sistema (en unidades reducidas).
            !ro: densidad de part�culas (en unidades reducidas).
            real (kind=doblep) :: pl, pli, vol, rc, rc2, ro
            !Correccion_potencial: correci�n a la energ�a potencial necesaria por usar un radio de corte finito.
            !Par�metros de correcci�n de derivadas:
			real (kind=doblep) :: correccion_potencial, correccion_dpot, correccion_ddpot, dt, dt_medios, dt2
            
            !Par�metros del potencial de Lennard-Jones, en Julios y metros.
            !real (kind=doblep), parameter :: eps, sigma

            !Se define pi como un par�metro constante.
            real (kind=doblep), parameter :: pi = 3.14159265359d00

        end module variables_comunes


