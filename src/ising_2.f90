program issing_1
	use ziggurat
	implicit none
	
	logical :: es
	integer :: seed, i, j, nising, niter,nmatriz,n,f,c,f1,c1,M,cont_no_cambia,cont_dE_menor_que_uni,cont_dE_menor_que_cero,E
	!logical :: parar
	real:: suma,x, dE, beta, KT,p,Mmedia,Emedia,Mcuadrado,Ecuadrado,desvM,desvE
	real, parameter:: Jota=1, Ho=0, K=1, T=4
	!real (kind=8), allocatable :: 
	integer, allocatable:: S(:,:),Mn(:),En(:),aceptado(:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!! Revisar formulas 
!! Agregar el espejo 
!!  Agregar el significado de J !!  
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
allocate(S(nmatriz,nmatriz), Mn(niter),En(niter),aceptado(niter))
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
   end do 
end do

do f=10,nising+9 
   write(16,*) (S(f,c),c=10,nising+9)
   end do 
close(16)

cont_no_cambia=0
cont_dE_menor_que_cero=0
cont_dE_menor_que_uni=0


write(5,*) "#   n // T[K] // f // c // M // E"
do n=1,niter
        !! lo corro para que no de el 1,1 
        !! cuando lo ponga en forma de toroide usar 1       
        f=nint(uni()*nising)+10 !init redondea 
        c=nint(uni()*nising)+10
        !! cambio un spin en la muestra 
        S(f,c)=-S(f,c)
        !! calculo el delta E 
        dE=2*Jota*S(f,c)*(S(f+1,c)+S(f,c+1)+S(f-1,c)+S(f,c-1))-Ho*S(f,c) 
        !! Algoritmo de Markov 
        if (dE<=0)then !! agregar como un and
          !! Acepto el nuevo estado  
          cont_dE_menor_que_cero=cont_dE_menor_que_cero+1        
          aceptado(n)=1
        else 
          !! Boltzmann probability of flipping
          p=exp(-dE/KT) !!Boltzman
          if (uni()<p) then  
              !! acepto el nuevo estado 
              cont_dE_menor_que_uni=cont_dE_menor_que_uni+1
              aceptado(n)=1 
          else 
             !! no acepto el nuevo estado 
             !! vuelvo al anterior
             cont_no_cambia=cont_no_cambia+1
             aceptado(n)=0
             S(f,c)=-S(f,c)
             dE=0
          end if               
        end if  
        
        if (mod(n,1000)==0) then   
        do f=10,nising+9
          write(7,*) (S(f,c),c=10,nising+9)  
        end do    
        end if                
       
          M=0
          E=0
          if (n==1) then 
          !! esta desde 10 porque esta corrida
            do f1=10,nising+9
              do c1=10,nising+9
                !! es el Hamiltoniano
                En(1)=E+-J*S(f1,c1)*(S(f1+1,c1)+S(f1,c1+1)+S(f1-1,c1)+S(f1,c1-1))-Ho*S(f1,c1)      
                Mn(1)=M+S(f1,c1)
              end do 
            end do         
          else 
               En(n)=dE+En(n-1)            
               Mn(n)=2*S(f,c)+Mn(n-1)
          end if             
            
         write(5,fmt="(i7,x,f10.3,x,i5,x,x,i5,x,x,i5,x,x,i5,x,x,i5)") n, T,f,c,aceptado(n), Mn(n), En(n)
end do 
close(5)

        do f=10,nising+9
          write(7,*) (S(f,c),c=10,nising+9)  
        end do 
        close(7)

Mmedia=0
Emedia=0
Ecuadrado=0
Mcuadrado=0
do i=1,niter
    Emedia=Emedia+En(i)        
    Mmedia=Mmedia+Mn(i)      
    Mcuadrado=Mcuadrado+Mn(i)**2
    Ecuadrado=Ecuadrado+En(i)**2 
end do 

Emedia=Emedia/dble(niter*nising**2)
Mmedia=Mmedia/dble(niter*nising**2)
Mcuadrado=Mcuadrado/dble(niter*nising**2)
Ecuadrado=Ecuadrado/dble(niter*nising**2)

desvE=(Ecuadrado-Emedia**2)**0.5
desvM=(Mcuadrado-Mmedia)**0.5

print *, " " 
print *, "Temperatura [K]:", T 
print *, " "
print '(" dE<0  %:" f10.4)', 100*dble(cont_dE_menor_que_cero)/dble(niter)
print *, " "
print '(" uni()<p %:" f10.4)',100*dble(cont_dE_menor_que_uni)/dble(niter)
print *, " --------------------"
print  '(" Cambia %:", f10.4)', 100*dble(niter-cont_no_cambia)/dble(niter)
print *, " "
print  '(" No cambia %:" f10.4)', 100*dble(cont_no_cambia)/dble(niter)
print *, " --------------------"
print *, " " 
print '(" Energia Media por Spin:" f10.5)', Emedia
!print *, " "
!print '(" Energia cuadrado por Spin" f10.4)', Ecuadrado
print *, " " 
print '(" Desviacion de E:" f10.4)',desvE
print *, " "   
print '(" M Media por Spin:" f10.4)', Mmedia
print *," " 
!print '(" M cuadrado por Spin" f10.4)', Mcuadrado
!print *, " " 
print '(" Desviacion de M: " f10.4)', desvM
print *, " --------------------"
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
