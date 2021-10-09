program issing_1
	use ziggurat
	implicit none
	
	logical :: es
	integer :: seed, i, j,L, niter,nmatriz,n,f,c,f1,c1,M,cont_no_cambia,cont_dE_menor_que_uni,cont_dE_menor_que_cero,E
	real (kind=8):: suma,x, dE, beta, KT,p,Mmedia,Emedia,Mcuadrado,Ecuadrado,desvM,desvE,cambia,no_cambia,Ho, T,Jota
	real, parameter::  K=1
	integer, allocatable:: S(:,:),Mn(:),aceptado(:)
	real (kind=8), allocatable :: En(:) !por que multiplico por 0.5 es real 
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
open(11, file = 'input.dat',status='old')
open(12, file= 'output.dat', status='old')
open(15, file= 'Magnetizacion.dat',status='old')

!open(8, file='M_T.dat',status='old')
!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!11
!!! inicializo 
read(11,*) niter
read(11,*) T
read(11,*) L 
read(11,*) Jota
read(11,*) Ho
close(11)
print *,  " " 
print *, "Numero de Iteraciones : ", niter 
print *,  " " 
print  '(" Temperatura:" f10.4)', T
print *,  " " 
print *, "Tamaño de matriz : ", L
print *,  " " 
print  '(" Valor de J:" f10.4)', Jota
print *, " "
print  '(" Campo externo:" f10.4)', Ho



!!!!!
!!! alocate 
allocate(S(0:L+1,0:L+1), Mn(niter),En(niter),aceptado(niter))
!!!!!!!

!!!!!!!!!!!!!!!!!
!!! Modelo de Ising
inquire(file='S.dat',exist=es)
	if(es) then
		open(unit=17,file='S.dat',status='old')
		print *," "
	    print *,"Se leyeron los espines de S.dat"
		do f=1,L
            read(17,*) (S(f,c) , c=1,L) 		    
		end do
	    close(17)
	else
	  print *," " 
	  print *,"Se genera una nueva matriz de spines S " 
      do i=1,L
        do j=1,L 
          !S(i,j)=1
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
	end if


!! espejo de la matriz
S(0,:)=S(L,:)
S(:,0)=S(:,L)
S(:,L+1)=S(:,1)
S(L+1,:)=S(1,:)


KT=K*T


cont_no_cambia=0
cont_dE_menor_que_cero=0
cont_dE_menor_que_uni=0
 
write(2,*)  "#1   n // #2 M // #3 E // T: ", T

do n=1,niter         
        f=nint(uni()*(L-1))+1 !init redondea 
        c=nint(uni()*(L-1))+1
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
             !dM=0
          end if               
        end if  
        
                  
          !! suma de primeros vecinos      
          M=0
          E=0
          if (n==1) then 
              do f1=1,L
                do c1=1,L
                !! es el Hamiltoniano
                E=E-Jota*S(f1,c1)*(S(f1+1,c1)+S(f1,c1+1)+S(f1-1,c1)+S(f1,c1-1))-Ho*S(f1,c1)      
                M=M+S(f1,c1)
                end do 
              end do
             
             
              En(1)=0.5*E !!porque cuento el doble 
              Mn(1)=M
                      
          else 
               En(n)=dE+En(n-1)            
               if (aceptado(n)==1) then 
               Mn(n)=2*S(f,c)+Mn(n-1)
               else 
               Mn(n)=Mn(n-1) !!dM=0
               end if 
          end if
            
          !! espejo de la matriz
          S(0,:)=S(L,:)
          S(:,0)=S(:,L)
          S(:,L+1)=S(:,1)
          S(L+1,:)=S(1,:)           
          
                    
        !if (mod(n,1000)==0) then   
        !   do f=1,L
        !     write(7,*) (S(f,c),c=1,L)            
        !  end do            
        !end if    
                        
         
         
         write(12,fmt="(i7,x,i5,x,x,f10.5,x,x,f10.5)") n,aceptado(n),dble(Mn(n))/dble(L**2),En(n)/dble(L**2)             
             
end do 

!!grabo la matriz final
         if (es) then 
          open(17, file='S.dat',status='old')
          do f=1,L
            write(17,*)  (S(f,c),c=1,L)  
          end do 
          close(17)  
          print *, " "
          print *, "Se graba la nueva matriz en S.dat"
         else           
          open(17, file='S.dat',status='new')
          do f=1,L
            write(17,*)  (S(f,c),c=1,L)    
            write(17,*) " "          
          end do 
          print*," " 
          print*,  "Se crea un nuevo archivo S.dat"
          close(17)  
         end if 
Emedia=0
Emedia=0
Ecuadrado=0
Mcuadrado=0
do i=1,niter
    Emedia=Emedia+En(i)        
    Mmedia=Mmedia+Mn(i)      
    Mcuadrado=Mcuadrado+Mn(i)**2
    Ecuadrado=Ecuadrado+En(i)**2 
end do 


Emedia=Emedia/dble(niter)
Mmedia=Mmedia/dble(niter)
Mcuadrado=Mcuadrado/dble(niter)
Ecuadrado=Ecuadrado/dble(niter)

desvE=(Ecuadrado-Emedia**2)**0.5
desvM=(Mcuadrado-Mmedia**2)**0.5

cont_dE_menor_que_cero=dble(cont_dE_menor_que_cero)/dble(niter)
cont_dE_menor_que_uni=dble(cont_dE_menor_que_uni)/dble(niter)
cambia=dble(niter-cont_no_cambia)/dble(niter)  
no_cambia=dble(cont_no_cambia)/dble(niter)

!print *, " " 
!print *, "Temperatura [K]:", T 
!print *, " "
!print '(" dE<0  %:" f10.4)', 100*dble(cont_dE_menor_que_cero)/dble(niter)
!print *, " "
!print '(" uni()<p %:" f10.4)',100*dble(cont_dE_menor_que_uni)/dble(niter)
!print *, " --------------------"
!print  '(" Cambia %:", f10.4)', 100*dble(niter-cont_no_cambia)/dble(niter)
!print *, " "
!print  '(" No cambia %:" f10.4)', 100*dble(cont_no_cambia)/dble(niter)
!print *, " --------------------"
!print *, " " 
!!print *, " "
!!print '(" Energia cuadrado por Spin" f10.4)', Ecuadrado
!print *, " " 
!print '(" Desviacion de E:" f10.4)',desvE
!print *, " "   
!print '(" M Media por Spin:" f10.4)', Mmedia
!print *," " 
!!print '(" M cuadrado por Spin" f10.4)', Mcuadrado
!!print *, " " 
!print '(" Desviacion de M: " f10.4)', desvM
!print *, " --------------------"

write(15,*) " #1 Emedi // #2 desvE // #3 Mmedia // #4 desM // #5 cambia // #6 no_cambia // T: ", T
write(15,fmt="(6f10.5)") Emedia, desvE, Mmedia, desvM, cambia, no_cambia
close(15)
close(12)    

print *, " "
print *, "Se grabaron los vectores M y E en output.dat"
print *, " "
print *, "Se grabaron las variables de interes en Magnetizacion.dat"
print *, " "
print *," " 
print *,"Fin del programa"    


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
