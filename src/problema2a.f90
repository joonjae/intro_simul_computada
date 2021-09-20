program simple 
    use ziggurat
    implicit none
    logical :: es
    integer :: seed, i, imatch
    integer :: Ntest=10000
    real :: pi = 3.14159

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

        imatch = 0
        do i = 1, Ntest
                if(uni()<(pi/4))then
                        imatch = imatch + 1
                end if
        end do

        print *,"Pi=Ncirc*4/Ntot"
        print *,"Pi=",real(imatch*4)/real(Ntest)

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
