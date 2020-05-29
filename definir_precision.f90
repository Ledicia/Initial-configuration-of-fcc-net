!Ledicia Díaz Lago
	!Definimos este modulo para indicar la precision con la que se trabaja:
	module definir_precision
    
    	implicit none
        integer,parameter:: entero=SELECTED_INT_KIND(9)
        integer, parameter :: doblep=SELECTED_REAL_KIND(15,137)
        
	end module definir_precision