program problema_ising
    ! tareas y funciones de otro archivo
    use ziggurat
    use mymodule
    implicit none

    ! Parametros
    !parameter(L=6) ! Parametro para la dimension de la matriz

    ! Variables
    logical :: es
    integer :: seed, i , j
    integer , allocatable :: mapa(:,:)
    integer :: L = 6, dim_mapa
    integer , allocatable :: init_pos(:)

    ! Inicio generador d nro random
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
    !----------------------------------------------------


    allocate(mapa(L,L))
    dim_mapa = size(shape(mapa))
    allocate(init_pos(dim_mapa))

    call init_mapa(mapa,L)
    do j=1,L
        write(*,*) (mapa(i,j), i=1,L)
    end do

    call rpos_mapa(mapa,L,init_pos,dim_mapa)
    print *,init_pos

    ! Escribir la última semilla para continuar con la cadena de numeros aleatorios 
    open(unit=10,file='seed.dat',status='unknown')
    seed = shr3() 
    write(10,*) seed
    close(10)
end program problema_ising

