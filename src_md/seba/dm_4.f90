program dinamica_1
	use ziggurat
    use globals
    use module_dm
	implicit none
	
	logical :: es
	integer :: seed, i,j,cont,parar,k,ncorte,ncorte1,ncorte2,ok
	real (kind=8):: Fij,densidad,a,a1,a2
	real (kind=8):: d,deltat,dt,dt2,tfinal,masa,radio_corte,b
	real, parameter::  sigma=1, ep=1,kb=1
	real (kind=8), allocatable :: fuera(:,:)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! 
!!  
!!    
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


!!! abro los archivos 
open(11, file = 'input.dat',status='old')
 
read(11,*) L 
read(11,*) N !numero de particulas
read(11,*) ncorte
read(11,*) ncorte1
read(11,*) ncorte2
read(11,*) densidad
read(11,*) deltat
read(11,*) Temperatura
close(11)
  
dt=deltat
masa=1
tfinal=10

!!Cambio el L 
L=dble(N)/densidad !!falta agregar sigma 
L=L**0.3333    !! no funciona dble(1/3)


!!!!!
!!! alocate 
allocate(R(3,N))
allocate(Velocidad(3,N))
allocate(F(3,N))
allocate(delta_R(3))
allocate(fuera(2,N))

!!!! Matriz de posiciones iniciales
call init_posiciones(sigma,kb)   
   
!! Potencial de Lennard-Jones


radio_corte=2.5*sigma
call potencial_lennard_jones(sigma,ep,radio_corte,masa,tfinal)


print*, " "
print*, " Termalicacion T=1.5*sigma/kb"
print*, " "

open(unit=51, file='movie_1.vtf', status='old')
write(51,*) "atom 0:200 radius 0.5 name S"
open(unit=52, file='Termalizacion.dat', status='old')


a2=0
dt=deltat
Etotal=0
cont=1
parar=1
do while (parar==1) !!hacer un vector o como el temperatura
F(:,:)=0
fuera(:,:)=0
Ec=0
V=0

a2=0
do i=1,N 
   Fij=0
  do j=i+1,N
          
    delta_R(1)=R(1,i)-R(1,j)
    delta_R(2)=R(2,i)-R(2,j)
    delta_R(3)=R(3,i)-R(3,j)
    
    !! condicones de contorno
    !! Analizar  
    delta_R(1)=delta_R(1)-L*int(2*delta_R(1)/L)
    delta_R(2)=delta_R(2)-L*int(2*delta_R(2)/L)
    delta_R(3)=delta_R(3)-L*int(2*delta_R(3)/L)
        
    !! distancia
    d=(delta_R(1)**2+delta_R(2)**2+delta_R(3)**2)**0.5
    

     
    if (d<radio_corte) then 
      V=V+4.0*ep*((sigma/d)**12-(sigma/d)**6)! -1/60*ep !!ultimo ternmino de radio de corte
              
          
          do k=1,3  
            fij=(24*ep/d**2)*(-(sigma/d)**6+2*(sigma/d)**12)
            F(1,i)=F(1,i)+fij*delta_R(k) !!fuerza en x,y,z
            F(1,j)=F(1,j)-fij*delta_R(k)
         end do 
         
     end if  
 
  end do     
end do 


do i=1,N 
      do k=1,3 
       R(k,i)=R(k,i)+Velocidad(k,i)*dt+(F(k,i)*0.5*dt**2)/masa !! /m   pero es uno 
       Velocidad(k,i)=Velocidad(k,i)+(F(k,i)*dt)/masa              
        
       if (R(k,i)>L) then 
           R(k,i)=R(k,i)-L
           !R(k,i)=R(k,i)-int(R(k,i)/L)*L
           else
         if (R(k,i)<0) then  
           R(k,i)=R(k,i)+L
           !R(k,i)=-R(k,i)
           !R(k,i)=R(k,i)-int(R(k,i)/L)*L
         end if 
        end if  
      end do       
      

end do   
  !Ec=Ec/(N-1) !! energia por par 
  V=V/(N-1) !! energia por par 
  
  !Etotal=Ec+V !!energia por par 

 !!No me funciona 
 write(51,*) "atom 0:200 radius 0.5 name S"
 write(51,*) "timestep"
 write(51,*) 
    do i=1,N
       write (51,*) R(:,i)
    end do 
 
  write(52,fmt='(i12,e12.3,x,f15.6,x,e12.10,x,e12.10)') cont,dt,V
  print  '(i12,e12.3,x,f15.6,x,e12.5,x,e12.5)', cont, dt, V
    
  dt=dt+deltat
  if ((dt>=tfinal).or.(cont==ncorte1))then 
    parar=-1
  end if 
  cont=cont+1
  
end do ! Termalizacion   
close(51)  
close(52)  


!print*, " "
!print*, " Corrida No funciona VERLET"
!print*, " "

open(unit=55, file='movie_2.vtf', status='old')
write(55,*) "atom 0:200 radius 0.5 name S"
open(unit=56, file='Corrida.dat', status='old')


a2=0
dt=deltat
Etotal=0
cont=1
parar=1
do while (parar==1) !!hacer un vector o como el temperatura
F(:,:)=0
fuera(:,:)=0
Ec=0
V=0

do i=1,N 
   Fij=0
  do j=i+1,N
    
    do k=1,3
        delta_R(k)=R(k,i)-R(k,j)
    end do 
        
    !! condicones de contorno
    do k=1,3
       delta_R(k)=delta_R(k)-L*int(2*delta_R(1)/L)
    end do 

        
    !! distancia
    d=(delta_R(1)**2+delta_R(2)**2+delta_R(3)**2)**0.5
    
        
    if (d<radio_corte) then 
      V=V+4.0*ep*((sigma/d)**12-(sigma/d)**6)-1/60*ep !!ultimo ternmino de radio de corte          
    end if  
 
    do k=1,3
       R(k,i)=R(k,i)+Velocidad(k,i)*dt+(F(k,i)*0.5*dt**2)/masa !! /m   pero es uno 
       Velocidad(k,i)=Velocidad(k,i)+0.5*(F(k,i)*dt)/masa !! Velocidad(t+dt/2)             
       
         
       fij=(24*ep/d**2)*(-(sigma/d)**6+2*(sigma/d)**12)
       F(1,i)=F(1,i)+fij*delta_R(k) !!fuerza en x,y,z
       F(1,j)=F(1,j)-fij*delta_R(k)
       
       
       Velocidad(k,i)=Velocidad(k,i)+0.5*(F(k,i)*dt)/masa
       
                   
       Ec=Ec+0.5*masa*Velocidad(k,i)**2 
       
       if (R(k,i)>L) then 
           R(k,i)=R(k,i)-L
           !R(k,i)=R(k,i)-int(R(k,j))*L
        end if 
        
       if (R(k,i)<0) then    
           R(k,i)=R(k,i)+L
           !R(k,i)=R(k,i)+int(R(k,j))*L
       end if 
        
     end do     
     
  end do     
end do 


  Ec=Ec/(N-1) !! energia por par 
  V=V/(N-1) !! energia por par 
  
  Etotal=Ec+V !!energia por par 

 !!No me funciona 

 write(55,*) "timestep"
 write(55,*) 
    do i=1,N
       write (55,*) R(:,i)
    end do 
 
  !write(56,*) cont,dt,V,Ec,Etotal 
  write(52,fmt='(i12,e12.3,x,e12.10,x,e12.10,x,e12.10,x,e12.10,x,e12.10)')  & 
                cont,dt,V,Ec,Etotal
  print  '(i12,e12.3,x,f15.6,x,e12.5,x,e12.5)', cont, dt, V, Ec, Etotal
    
  dt=dt+deltat
  dt2=dt/2
  if ((dt>=tfinal).or.(cont==ncorte2))then 
    parar=-1
  end if 
  cont=cont+1
  
end do ! Termalizacion   
close(55)  
close(56)  


!open(42, file='fuera.dat',status='old')
! do i=1,N
!   if (fuera(1,i)==1) then 
!      a=0
!      !write(42,*) fuera(1,i), fuera(2,i)
!   end if 
!end do 

print *,  " "
print *,  "----------------------------"
print *,  " " 
print *, "Numero de Iteraciones : ", cont-1
print *,  " " 
print  '(" Tamaño de la caja recalculado :" f10.4)', L
print *,  " " 
print *, "Numero de particulas", N
print *, " "
print  '(" Densidad :" f10.4)', densidad
print *, " " 
print  '(" Delta t  :" e12.3)', deltat !!pasarlo  a exponencial
print  *," " 
print  '(" Temperatura  :" f10.4)', Temperatura
print *, " "
print *, " Chequeo de distacias ", a,a1,a2
print *, " ---------------- " 
print *, "Fin del programa"    
print *, " ---------------- " 
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

 
       
end program dinamica_1
