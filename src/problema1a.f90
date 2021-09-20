program simple 
    use ziggurat
    implicit none
    logical :: es
    integer :: seed,i ,j,k
    real :: x , y

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

! Ej: Número random en [0,1]: uni()


!!        do i = 1, 500
!!            print *,uni()
!!        end do
        do i = 1, 100
                x = uni()
                y = uni()
                if( x < y )then
                        print *,i,": X<Y ", "(", x ,"<", y ,")"
                end if
        end do

!! 
!! EDITAR AQUI 
!! 


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


end program simple
