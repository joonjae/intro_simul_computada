program simple 
    use ziggurat
    implicit none
    logical :: es
    integer :: seed,N,i ,j,k
    real :: x , y
    real , allocatable :: vector_x(:), vector_y(:), fdx(:)

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

            allocate(vector_x(N))
            allocate(vector_y(N))
            allocate(fdx(N))
            

            open(unit=101,file='output.dat',status='unknown')   ! abro o creo archivo output.dat
            j = 0
            do i = 1, N
                    x = uni()*2.
                    y = uni()

                    vector_x(i) = x
                    vector_y(i) = y
                    fdx(i) = (x - 1.)**2

                    !! Rutina del ejercicio punto C 
                    if( (y <= fdx(i)) )then
                            write(101,*) x,y
                            j = j+1
                    end if
                    !! Fin punto C del ejercicio
            end do
            close(101)                                          ! cierro archivo output.dat

            print *,"Porcentaje de salida: ",real(j)/real(N)
            !!call execute_command_line('gnuplot output.dat')
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
