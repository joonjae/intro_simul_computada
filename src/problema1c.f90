program simple 
    use ziggurat
    implicit none
    logical :: es
    integer :: seed,N,i ,j,k
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


    inquire(file='input.dat',exist=es)
    if(es)then
            open(unit=100,file='input.dat',status='old')
            read(100,*) N
            close(100)
            print *,"N=",N

            open(unit=101,file='output.dat',status='unknown')   ! abro o creo archivo output.dat
            j = 0
            do i = 1, N
                    x = uni()
                    y = uni()
                    if( x < y )then
                            print *,i,": X<Y ", "(", x ,"<", y ,")"
                            write(101,*) x,y
                            j = j+1
                    end if
            end do
            close(101)                                          ! cierro archivo output.dat

            print *,"Porcentaje de salida: ",real(j)/real(N)
    else
            print *,"Finaliza el programa porque no existe input.dat"
    end if

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
