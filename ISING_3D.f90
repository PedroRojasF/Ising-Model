!Programa que implementa el algoritmo de Monte Carlo
!y Metrópoli para resolver un problema de Ising 3D

MODULE DATOS
IMPLICIT NONE
SAVE !El comando "SAVE" conserva los elementos de un subprograma 
      !después de que se ejecutan las instrucciones RETURN o END, 
      !lo que evita que se vuelvan indefinidos.
      
!________Declaración de variables_______________________________________________________________________
!Temperatura inicial/Temperatura final/ Diferencial de temperatura
   REAL,PARAMETER   :: T0=0.5, Tf=10.0, dT=0.005   
   
!Tamaño del 'enrejado' cúbico / Dimensión 
   INTEGER,PARAMETER:: L=20, D=3 
   
!# Espínes vecinos/ #Barrido de'Monte Carlo/ #Barrido.MC-equilibro
   INTEGER,PARAMETER:: Nespv=6, MCb=100, MCbeq=100    
    
!Número de espínes
   REAL,PARAMETER   :: Nesp=L**D, nr1=1.0/(MCb*Nesp*Nespv), nr2=nr1/MCb
   
!Coordenadas de la posición de los espínes/ Contadores de bucle
   INTEGER:: x,y,z,i,j,k                
   INTEGER:: up, down, right, left, above,below
END MODULE
!_______________________________________________________________________________________________________
PROGRAM ISING3D  !____PROGRAMA PRINCIPAL________________________________________________________________
  USE DATOS
IMPLICIT NONE
!-------Declaración de variables------------------------
!Temperatura/ Beta=1/k_{B}T/ Beta^2
  REAL:: T , beta , beta2 
!Variables de acumulación en el bucle
  REAL:: summ0 , summ1 , summ2 , summ3 , summ4 
!Energía/Magnetización/Energía media/ Magnetiz. media/
! calor específico/ Susceptibilidad
  REAL:: E , M, PE , PM , Cespcf, Xu    
  REAL:: ENERGIA , MAGNETIZACION ! external functions
  INTEGER:: barrido, c       !Contadores de bucle
  INTEGER:: spin(L,L,L) = 1 !Configuración inicial:Espines hacia arriba
!_______Guardando resultados____________________________________________________________________________
OPEN(unit=1,file='Resultado.dat')
T = T0 !Inicialización de la temperatura 
 c = 0
10 FORMAT(6f15.6)
WRITE(1,*) " °  Temperatura     Temperatura     Energía     Magnetización      Calor      Susceptibilid  "
WRITE(1,*) '      inversa                        media          media        específico     magnética'
WRITE(1,*) " °    T^(-1)            T              E              M             C_v            Chi   "
WRITE(1,*) "  _______________||____________||_____________||______________||____________||______________"
!_________________________________________________________________________________________________________
DO WHILE (T<=Tf)            !Bucle principal de temperatura
T = T + dT; c = c + 1
summ0 = 0.; summ1 = 0.; summ2 = 0.; summ3 = 0.
 beta = 1./T ; beta2 = beta*beta 
   DO barrido = 1, MCbeq       !Bucle para el equilibrio
       CALL METROPOLIS(spin,beta)
   END DO
   DO barrido = 1, MCb       !Iteraciones de Monte Carlo
       CALL METROPOLIS(spin,beta)
E = ENERGIA(spin)
M = MAGNETIZACION(spin)
      summ0 = summ0 + E      ! E   acumulador
      summ1 = summ1 + E*E    ! E^2 acumulador
      summ2 = summ2 + M      ! M   acumulador
      summ3 = summ3 + M*M    ! M^2 acumulador
   END DO
!_______________________________________________________________________________________________________
!Promedios termodinámicos en cada valor de temperatura
PE = summ0*nr1 !Energía promedio
PM = summ2*nr1 !Magnetización promedio
 Cespcf = (nr1*summ1 - nr2*summ0*summ0)*beta2    !Calor específico
 Xu      = (nr1*summ3 - nr2*summ2*summ2)*beta    !Susceptibilidad
WRITE (1,10)beta, T, PE, PM, Cespcf, Xu
!WRITE (1,10)beta, T, PE, PM, Cespcf, Xu
WRITE (*,*) " >>> Paso de temperatura ", c
END DO
END PROGRAM
!_______________________________________________________________________________________________________

SUBROUTINE METROPOLIS(spin,betaa)!Algoritmo de Metropoli
 USE DATOS
IMPLICIT NONE
REAL, INTENT(IN) :: betaa              !betaa = 1/(k_{B}*T)
REAL             :: R, dE              !Diferencia de energía
REAL             :: x0,y0,z0,aleat     !Números aleatorios                        
INTEGER,INTENT(INOUT):: spin(L,L,L)    !Espines
INTEGER:: sbox                            !Caja para el espín
!________________________________________________________________________________________________________
     DO i = 1, L
        DO j = 1, L
           DO k = 1, L
              CALL random_number(x0)
              CALL random_number(y0)
              CALL random_number(Z0)
                   x = int(1+(L-1)*x0)    !Coordenada en X
                   y = int(1+(L-1)*y0)    !Coordenada en Y
                   z = int(1+(L-1)*z0)    !Coordenada en Z
                   sbox = spin(x,y,z)     !Espín 
!_________________ Condiciones de contorno periódicas___________________________________________________
                  IF (x==1) THEN
                     left  = spin(L  ,y,z)
                     right = spin(2  ,y,z)
                     ELSE IF (x==L) THEN
                     left  = spin(L-1,y,z)
                     right = spin(1  ,y,z)   
                     ELSE
                     left  = spin(x-1,y,z)
                     right = spin(x+1,y,z) 
                  END IF 
               !________________________________________________________________________________________
                  IF (y==1) THEN
                     up   = spin(x, 2 ,z)
                     down = spin(x, L ,z)
                     ELSE IF (y==L) THEN
                     up   = spin(x, 1 ,z)
                     down = spin(x,L-1,z)   
                     ELSE
                     up   = spin(x,y+1,z)
                     down = spin(x,y-1,z) 
                  END IF   
                !________________________________________________________________________________________
                  IF (z==1) THEN
                     above= spin(x,y, 2 )
                     below= spin(x,y, L )
                     ELSE IF (z==L) THEN
                     above = spin(x, y ,1)
                     below = spin(x,y,L-1)   
                  ELSE
                     above = spin(x,y,z+1)
                     below = spin(x,y,z-1) 
                  END IF 
!________________________________________________________________________________________________________     
              R  = up + down + right + left + above + below  
              dE = 2*sbox*R   
CALL random_number(aleat)                  !"Tiramos el dado"   
              IF (dE<0.0) THEN   
                 sbox = -sbox              !Volteamos al espín
                 ELSE IF (aleat < exp(-dE*betaa)) THEN
                 sbox = -sbox              !Nuevamente volteamos al espín
              END IF
              spin(x,y,z) = sbox           !Retornamos el espin ésta vez sin voltearlo    
           END DO
        END DO
     END DO
END SUBROUTINE METROPOLIS
!________________________________________________________________________________________________________
REAL FUNCTION ENERGIA(spin)
  USE DATOS
IMPLICIT NONE
REAL              :: Enrg,R       !Energía inicial
INTEGER,INTENT(IN):: spin(L,L,L)
Enrg=0.
        DO x = 1, L
          DO y = 1, L
            DO z = 1, L
             ! Condiciones periódicas de contorno
                  IF (x==1) THEN
                     left  = spin( L ,y,z)
                     right = spin( 2 ,y,z)
                     ELSE IF (x==L) THEN
                     left  = spin(L-1,y,z)
                     right = spin(1  ,y,z)   
                     ELSE
                     left  = spin(x-1,y,z)
                     right = spin(x+1,y,z) 
                  END IF 
                  !______________________________________________________________________________________
                  IF (y==1) THEN
                     up   = spin(x, 2 ,z)
                     down = spin(x, L ,z)
                     ELSE IF (y==L) THEN
                     up   = spin(x, 1 ,z)
                     down = spin(x,L-1,z)   
                     ELSE
                     up   = spin(x,y+1,z)
                     down = spin(x,y-1,z) 
                  END IF   
                  !______________________________________________________________________________________
                  IF (z==1) THEN
                     above= spin(x,y, 2 )
                     below= spin(x,y, L )
                     ELSE IF (z==L) THEN
                     above = spin(x, y ,1)
                     below = spin(x,y,L-1)   
                  ELSE
                     above = spin(x,y,z+1)
                     below = spin(x,y,z-1) 
                  END IF    
!________________________________________________________________________________________________________
              R  = up + down + right + left + above + below    
              Enrg = Enrg - spin(x,y,z)*R  !Acumulación de energía
            END DO
          END DO
        END DO  
   ENERGIA = Enrg !Energía     
END FUNCTION
!________________________________________________________________________________________________________
REAL FUNCTION MAGNETIZACION(spin)
  USE DATOS  
IMPLICIT NONE
INTEGER,INTENT(IN):: spin(L,L,L)  !Espín
MAGNETIZACION = sum(spin)         ! Magnetización
END FUNCTION
 






