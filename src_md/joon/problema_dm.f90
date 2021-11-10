program problema_dm
    ! tareas y funciones de otro archivo
    use ziggurat
    use mymodule
    use globals
    implicit none

    ! Variables
    logical :: es
    integer :: seed, i , j
    integer :: N_MC

    ! Inicio generador nro random
    inquire(file='seed.dat',exist=es)
    if(es) then
        open(unit=10,file='seed.dat',status='old')
        read(10,*) seed
        close(10)
    !    print *,"  * Leyendo semilla de archivo seed.dat"
    else
        seed = 24583490
    end if
    call zigset(seed)
    !----------------------------------------------------

    !----------------------------------------------------
    ! Inicio lectura de input para el ising
    inquire(file='input.dat',exist=es)
    if(es) then
        open(unit=1,file='input.dat',status='old')
        read(1,*)
        read(1,*) N_MC, T, n, c, L
        close(1)
!        print *,"  * leyendo datos de input_ising.dat"
    else
        N_MC=1000
        T=1
        n=4
        c=3
        L=10
    end if
    close(1)
    print *,'Datos input:'
    print *,'N_MC:',N_MC
    print *,'T:',T
    print *,'n:',n
    print *,'c:',c
    print *,''
    !----------------------------------------------------

    allocate(r(c,n), v(c,n), f(c,n))

    print *,'Calculo las posiciones iniciales:'
    call verlet_posiciones()
    
    call fuerza(1.,1.)
    print *,'Fuerzas:'
    do i=1,n
        print *,f(:,i)
    end do

    ! Escribir la última semilla para continuar con la cadena de numeros aleatorios 
    open(unit=10,file='seed.dat',status='unknown')
    seed = shr3() 
    write(10,*) seed
    close(10)
end program problema_dm

