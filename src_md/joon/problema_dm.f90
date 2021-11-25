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
    integer :: i_mc

    integer :: file_input = 1
    integer :: file_output = 2

    real :: v_sist
    
    real :: tic, toc

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

!    !----------------------------------------------------
!    ! Inicio lectura de input
!    inquire(file='input.dat',exist=es)
!    if(es) then
!        open(unit=file_input,file='input.dat',status='old')
!        read(file_input,*)
!        read(file_input,*) N_MC, T, n, c, L
!        close(file_input)
!!        print *,"  * leyendo datos de input_ising.dat"
!    else
!        N_MC=1000
!        T=1
!        n=4
!        c=3
!        L=10
!    end if
!    close(file_input)
!    print *,'Datos input:'
!    print *,'N_MC:',N_MC
!    print *,'T:',T
!    print *,'n:',n
!    print *,'c:',c
!    print *,''
!    !----------------------------------------------------
!    open(unit=file_output,file='output.dat',status='unknown')
!
!    allocate(r(c,n), v(c,n), f(c,n))
!    call write_conf(0)
!
!    print *,'Posiciones iniciales:'
!    call init_posiciones()
!    
!    do i_mc=1, N_MC
!        call write_conf(1)
!        
!        call fuerza(1.,1.,2.5,v_sist)
!        !do i=1,n
!        !    print *,'Fuerzas:',f(:,i)
!        !end do
!
!        !print *,'V sistema:',v_sist
!        write(file_output,*) v_sist
!    end do
!    call write_conf(2)
!
!    close(file_output)


    
    call cpu_time(tic)
    call setup_parameters
    call initial_condition
    call equilibrate
    call evolve
    
    call cpu_time(toc)
    write(*,*) "Total time taken by code = ",toc-tic
    
    ! Escribir la última semilla para continuar con la cadena de numeros aleatorios 
    open(unit=10,file='seed.dat',status='unknown')
    seed = shr3() 
    write(10,*) seed
    close(10)
end program problema_dm

