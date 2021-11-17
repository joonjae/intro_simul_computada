MODULE module_dm
    use ziggurat
    use globals
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: init_posiciones, potencial_lennard_jones, write_conf

CONTAINS

! Subrutina para inicio de las posiciones random
! (IN)  NA
! (OUT) NA
! (GLOBAL)  L: Longitud de la caja
!           r: posicion de las particulas
subroutine init_posiciones(sigma,kb)
    implicit none
    real (kind=4), intent(in) :: sigma
    real (kind=4), intent(in) :: kb
    integer i

    do  i=1,N        
        R(1,i)=L*uni()
        R(2,i)=L*uni() 
        R(3,i)=L*uni()   
        
          
        Velocidad(1,i)=rnor()*sqrt(Temperatura*sigma/kb) 
        Velocidad(2,i)=rnor()*sqrt(Temperatura*sigma/kb)
        Velocidad(3,i)=rnor()*sqrt(Temperatura*sigma/kb)                  
    end do            
    return  
end subroutine init_posiciones


! Subrutina para calculo de energia potencial y fuerza de las particulas
! (IN)  sigma: medida del diametro efectivo de las particulas
! (IN)  esp: epsilon, medida del pozo potencial (intensidad de la interaccion)
! (OUT) 
! (GLOBAL)  L: Longitud de la caja
!           r: posicion de las particulas
!           f: fuerza total sobre cada particula
subroutine potencial_lennard_jones(sigma,ep,radio_corte,masa,tfinal)
    implicit none
    real (kind=4), intent(in) :: sigma
    real (kind=4), intent(in) :: ep
    real (kind=8), intent(in) :: radio_corte
    real (kind=8), intent(in) :: masa

    integer :: i,j,N,cont,k,ncorte
    integer :: parar
	real (kind=8):: L,fij,densidad,a,a1,a2,Ec,Etotal
	real (kind=8):: d,dt,deltat,tfinal

    open(unit=32, file='Minimizacion.dat', status='old')
    !open(unit=40, file='Distancia.dat', status='old')
    open(unit=50, file='movie.vtf', status='old')
    write(50,*) "atom 0:200 radius 0.5 name S"
    open(unit=100,file='se_va.dat',status='old')
    print*, " "
    print*, " Minimizacion"
    print*, " "


    a=0
    parar=1
    cont=1
    do while (parar==1) 
    
        V=0 !potencial total 
        f = 0.
        
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
                    V=V+4.0*ep*((sigma/d)**12-(sigma/d)**6) -1/60*ep !!ultimo ternmino de radio de corte
                      
                    do k=1,3  
                        fij=(24*ep/d**2)*(-(sigma/d)**6+2*(sigma/d)**12)
                        F(1,i)=F(1,i)+fij*delta_R(k) !!fuerza en x,y,z
                        F(1,j)=F(1,j)-fij*delta_R(k)
                    end do 
                   
                end if  
                      
            end do     
        end do 
           
        a1=0 
        !!!minimizo el potencial inicial 
        do i=1,N 
           !! x=1/2*a*t^2 !!no hay velocidad unicial 
           do k=1,3 
               R(k,i)=R(k,i)+(F(k,i)*0.5*dt**2)/masa  !! /m   pero es uno 
               
               if (R(k,i)>L) then 
                   R(k,i)=R(k,i)-L
                   !R(k,i)=R(k,i)-int(R(k,i)/L)*L  
                   write(100,*) R(i,k) , L           
               end if 
               if (R(k,i)<0) then  
                   R(k,i)=R(k,i)+L
                   !R(k,i)=R(k,i)+(int(-R(k,i)/L)+1)*L
                   write(100,*) R(i,k)           
               end if  
           end do 
        end do   
         
        write(50,*) "timestep"
        write(50,*) 
        do i=1,N
            write (50,*) R(:,i)
        end do 
        
        write(32,*) cont,dt,V
        print '(i12,x,e12.3,x,e12.5)', cont, dt, V 
          
        if ((dt>=tfinal).or.(cont==ncorte))then 
            parar=-1
        end if 
        cont=cont+1
        dt=dt+deltat
    
    end do ! de tiempo Minimizacion   

    close(100)
    close(50)
    close(32)

    return
end subroutine potencial_lennard_jones

subroutine write_conf(mode)
    implicit none
    integer, intent (in)::mode
    integer             ::i
    
    select case(mode)
        case(0)
            open(unit=20,file="movie.vtf",status="unknown")
            write(20,*) "atom 0:99 radius 1 name Ar"
        case(1)
            write(20,*) "timestep"
            write(20,*)
            do i=1,n
                write(20,*) r(:,i)
            end do
        case(2)
            close(20)
    end select

end subroutine write_conf

END MODULE module_dm
