program problema_ising
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
        read(1,*) N_MC, T, N, D, L
        close(1)
!        print *,"  * leyendo datos de input_ising.dat"
    else
        N_MC=1000
        T=1
        N=4
        D=3
        L=10
    end if
    close(1)
    print *,N_MC, T, N, D
    !----------------------------------------------------

    allocate(R(D,N), V(D,N), F(D,N))

    call init_posiciones(D,N)
    print *, R(1,2)
!    call init_mapa(mapa,L)          ! inicio el mapa 
!
!    area = L*L
!    densidad = 1./real(area)    ! densidad de area
!    pasos_mc = 1./real(N_MC)       ! densidad o pasos de monte carlo
!
!    allocate(ip(L))             ! reservo memoria para suma vecinos
!    allocate(im(L))             ! reservo memoria para restar vecinos

    

    ! Escribir la última semilla para continuar con la cadena de numeros aleatorios 
    open(unit=10,file='seed.dat',status='unknown')
    seed = shr3() 
    write(10,*) seed
    close(10)
end program problema_ising

