program simple 
    use ziggurat
    use mymodule

    implicit none
    logical :: es
    integer :: seed , N , i , uni_hits , rexp_hits ,rnor_hits , k
    integer :: file_rexp , file_rnor
    real :: x , y , randx , randy , time_begin , time_end , time_proc_uni , time_proc_rexp , time_proc_rnor

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

            !! Simulacion con rutina "uni()"
            open(unit=101,file='output.dat',status='unknown')   ! abro o creo archivo output.dat
            uni_hits = 0
            call CPU_TIME(time_begin)
            do i = 1, N
                    x = uni()*2.
                    y = uni()
                    
                    if( funcion_x( x , y )==1 ) then
                        uni_hits = uni_hits + 1
                        write(101 ,*) x , y
                    end if
            end do
            call CPU_TIME(time_end)
            time_proc_uni = time_end - time_begin
            close(101)                                          ! cierro archivo output.dat

            !! Creacion de funcion para el plot de histograma
            !! El eje X va de 0 a 2, entonces el delta de x es
!            delta_x = 0.1
!            open(unit=102,file='histogram.dat',status='unknown')   ! abro o creo archivo output.dat
            
!            close(102)

            !! Simulacion con rutina "rexp()"
            file_rexp = 102
            open(unit=file_rexp,file='output_rexp.dat',status='unknown')   ! abro o creo archivo output.dat
            rexp_hits = 0
            call CPU_TIME(time_begin)
            i = 0
            do while ( i < N ) 
!                write(file_rexp,*) -rexp()+2,rexp()
!                i = i+1
                randx = rexp()
                if( randx>=0 .and. randx<=2)then
                    if(mod(i,2)==0)then
                        x = randx+0.177
                        !x = randx-1
                    else
                        x = -randx+1.813
                        !x = -randx+3
                    end if
                    if( x>0 .and. x<2)then
                        y = rexp()
                        if(y<1)then
                            if( funcion_x( x , y ) == 1 )then
                                    rexp_hits = rexp_hits + 1
                                    write(file_rexp ,*) x , y
                            end if
                            i = i+1
                        end if
                    end if
                end if
            end do
            call CPU_TIME(time_end)
            time_proc_rexp = time_end - time_begin
            close(file_rexp)                                          ! cierro archivo output.dat

            !! Simulacion con rutina "rnor()"
            file_rnor = 103
            open(unit=file_rnor,file='output_rnor.dat',status='unknown')   ! abro o creo archivo output.dat
            rnor_hits = 0
            call CPU_TIME(time_begin)
            i = 0
            do while ( i < N )
                randy = rnor()
                randx = rnor()
!                if(randy>0)then 
!                        if(randx>0)then
!                                write(file_rnor,*) randx-1,randy
!                        else
!                                write(file_rnor,*) randx+3,randy
!                        end if
!                        i = i+1
!                end if
                if(randy>=0)then
                    y = randy
                    if( randx > 0 ) then
                        !x = randx - 1.45 
                        !x = randx - 1.42
                        x = randx - 1.41
                        !x = randx - 1.35
                        if( x > 0 .and. x < 2 )then
                                if( funcion_x( x , y )==1 )then
                                        rnor_hits = rnor_hits + 1
                                        write(file_rnor ,*) x , y
                                end if
                                i = i+1
                        end if
                    else
                        !x = randx + 3.45 
                        !x = randx + 3.42 
                        x = randx + 3.41
                        !x = randx + 3.35
                        if( x < 2 .and. x > 0 )then
                                if( funcion_x( x , y )==1 )then
                                        rnor_hits = rnor_hits + 1
                                        write(file_rnor ,*) x , y
                                end if
                                i = i+1
                        end if
                    end if
                end if
            end do
            call CPU_TIME(time_end)
            time_proc_rnor = time_end - time_begin
            close(file_rnor)                                          ! cierro archivo output.dat


            print *,"Tiempo de computo por uni es ", time_proc_uni
            print *,"Porcentaje de salida uni: ",real(uni_hits)/real(N)
            print *,"Tiempo de computo por rexp es ", time_proc_rexp
            print *,"Porcentaje de salida rexp: ",real(rexp_hits)/real(N)
            print *,"Tiempo de computo por rnor es ", time_proc_rnor
            print *,"Porcentaje de salida rnor: ",real(rnor_hits)/real(N)


            !call execute_command_line("gnuplot plot_output.sh")
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
