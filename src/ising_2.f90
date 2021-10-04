program issing_1
	use ziggurat
	implicit none
	
	logical :: es
	integer :: seed, i, j, nising, niter,nmatriz,n,f,c,f1,c1,M,cont_no_cambia,cont_dE_menor_que_uni,cont_dE_menor_que_cero,E
	!logical :: parar
	real:: suma,x, dE, beta, KT,p
	real, parameter:: Jota=1, Ho=0, K=1, T=5 
	!real (kind=8), allocatable :: 
	integer, allocatable:: S(:,:),Magnetizacion(:),En(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!! Revisar formulas 
!! Agregar un contador para todos los casos 
!! Agregar el espejo 
!! Agregar el comentario del paper donde explica bien cada paso 
!! Agregar el significado de J 
!!
!! Analizar que pasa con la temperatura en paramagnetico 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!	 
![NO TOCAR] Inicializa generador de número random

	inquire(file='seed.dat',exist=es)
	if(es) then
		open(unit=10,file='seed.dat',status='old')
		read(10,*) seed
		close(10)
		print *,"  * Leyendo semilla de archivo seed.dat"
	else
		seed = 24583490
	end if
	call zigset(seed)
![FIN NO TOCAR]    

!!!!!!!!!! 
!!! abro los archivos 
open(1, file = 'input.dat',status='old')
!open(2, file = 'output.dat',status='old')
!open(3, file = 'distribucion.dat', status='old') 
!open(4, file=  'histograma.dat', status='old')
open(5, file= 'Magnetizacion.dat',status='old')
open(16, file= 'S.dat',status='old') !! Ojo porque uni=6 es pantalla no se puede usar 
open(7, file='S1.dat',status='old')
!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!! inicializo 
nising=20
nmatriz=50!! esto esta por ahora para no repetir la filas y columas 
!niter=1
read(1,*) niter
close(1)
print *,  " " 
print *, "Numero de Iteraciones : ", niter 

KT=K*T

!!!!!
!!! alocate 
allocate(S(nmatriz,nmatriz), Magnetizacion(niter),En(niter))
!!!!!!!

!!!!!!!!!!!!!!!!!
!!! Modelo de Ising 

!! ver como es la mejor forma de ordenar 
do i=1,nmatriz 
  do j=1,nmatriz 
   x=uni()
   if (x<=0.5) then 
   x=1
   else 
   x=-1
   end if
   S(i,j)=x 
   !write(16,*) "S(",i,",",j,"):", x 
   !print *, "i: ", i , "j: ", j, x
end do 
end do

do f=10,nising+9 
   write(16,*) (S(f,c),c=10,nising+9)
   !write(16,fmt="(i5)") (S(f,c), c=1, nising)  
   end do 
close(16)

cont_no_cambia=0
cont_dE_menor_que_cero=0
cont_dE_menor_que_uni=0
!print *, " " 
!print *, "Selecciono las filas y columas"

write(5,*) "#   n // T[K] // f // c // M // E"
do n=1,niter
   !do i=1,nissing
     ! do j=1,nissing 
        !! lo corro para que no de el 1,1 
        !! cuando lo ponga en forma de toroide usar 1       
        f=nint(uni()*nising)+10 !init redondea 
        c=nint(uni()*nising)+10
        !print *," "
        !print *, " Eleccion, fila: ", f, "columna :",c
        !! calculo el delta E 
        !dE=-Jota*S(f,c)*(S(f+1,c)+S(f,c+1)+S(f-1,c)+S(f,c-1))-Ho*S(f,c) 
        dE=2*Jota*S(f,c)*(S(f+1,c)+S(f,c+1)+S(f-1,c)+S(f,c-1))-Ho*S(f,c) !! sacado de matlab
        !print *, "fila: ", f , "columna: ", c 
        !print *, "Delta E : ", dE  
        !! Algoritmo de Markov 
        if (dE<=0)then !! agregar como un and
          !print *, "delta E <0, cambio de estado"
          !! cambio el estado 
          S(f,c)=-S(f,c)
          cont_dE_menor_que_cero=cont_dE_menor_que_cero+1
        
        else 
        !! Boltzmann probability of flipping
        !! prob = exp(-dE / kT);
          p=exp(-dE/KT)
          if (p<uni()) then !! no entiendo bien esta parte por que un estado del sistema
          !print *, "Se p<uni(), cambio de estado"
          cont_dE_menor_que_uni=cont_dE_menor_que_uni+1
          S(f,c)=-S(f,c)
          else 
          !print *, "No se cambia el spin"
          cont_no_cambia=cont_no_cambia+1
          end if
              
        end if  
          M=0
          E=0
          !! esta desde 10 porque esta corrida
          do f1=10,nising+9
            do c1=10,nising+9
              M=M+S(f1,c1)
              E=E+S(f1,c1)*(S(f1+1,c1)+S(f1,c1+1)+S(f1-1,c1)+S(f1,c1-1))-Ho*S(f1,c1)      
              !print *," fila: ",f1, "Columna", c1 
            end do 
         end do 
         Magnetizacion(n)=M
         En(n)=E      
         write(5,fmt="(i7,x,f10.3,x,i5,x,i5,x,x,i5,x,x,i5)") n, T,f,c, Magnetizacion(n), En(n)
  

end do 

close(5)

do f=10,nising+9
   write(7,*) (S(f,c),c=10,nising+9)  
end do 
close(7)


print *, " " 
print *, "Temperatura [K]:", T 
print *, " "
print '(" dE<0  %:" f10.4)', 100*dble(cont_dE_menor_que_cero)/dble(niter)
print *, " "
print '(" p<uni() %:" f10.4)',100*dble(cont_dE_menor_que_uni)/dble(niter)
print *, " "
print  '(" Cambia %:", f10.4)', 100*dble(niter-cont_no_cambia)/dble(niter)
print *, " "
print  '(" No cambia %:" f10.4)', 100*dble(cont_no_cambia)/dble(niter)
print *, " " 
print *, "Fin del programa"
print *, " "  


!! 
!! FIN FIN edicion
!! 
![No TOCAR]
! Escribir la última semilla para continuar con la cadena de numeros aleatorios 

		open(unit=10,file='seed.dat',status='unknown')
		seed = shr3() 
		write(10,*) seed
		close(10)
![FIN no Tocar]        

end program issing_1
