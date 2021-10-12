program problema_ising
    ! tareas y funciones de otro archivo
    use ziggurat
    use mymodule
    implicit none

    ! Parametros
    !parameter(L=6) ! Parametro para la dimension de la matriz

    ! Variables
    logical :: es
    integer :: seed, i , j, icolor
    integer :: L = 20 ! tamanio de la matriz. ej. L x L
    integer , allocatable :: mapa(:,:)
    integer , allocatable :: init_pos(:)
    integer :: dim_mapa

    integer :: ier,pgbeg ! para referenciar al inicio de pgplot

    integer, allocatable :: ip(:),im(:) ! para ir mapeando cuales son los vecinos
    real :: N0, U, T0, T, H
    integer :: N_MC, N_INTER
    integer :: imc, i_in
    integer :: i_area, area
    integer :: isum
    real :: energia_antigua, energia_nueva, delta_energ_beta, ising
    real :: sumep, summp, xetot, xe2tot, xmtot, xm2tot

    real :: ht, ut
    real :: r, prob_acept
    real :: densidad, pasos
    
    integer :: n_changes

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

    ! Preparo el pgplot
    ier = pgbeg(0,'/xserve',1,1)
    if(ier .ne. 1) STOP
    call pgenv(0,real(L+1),0,real(L+1),1,0)
    !----------------------------------------------------

    ! Inicio lectura de input para el ising
    ! N0: nro de spin
    ! U: Energia de interaccion
    ! T: Temperatura
    ! H: Campo magnetico
    ! N_MC: nro de iteraciones de Monte Carlo
    inquire(file='input_ising.dat',exist=es)
    if(es) then
        open(unit=1,file='input_ising.dat',status='old')
        read(1,*)
        read(1,*) N0 ,U ,T0, H
        read(1,*)
        read(1,*) N_MC, N_INTER
        close(1)
!        print *,"  * leyendo datos de input_ising.dat"
    else
        N0=0
        U=1
        T0=1
        H=0
        N_MC=1000
    end if
    close(1)
!    print *,"N0:",N0," U:",U ,"T:",T ,"H",H, "N_MC:",N_MC

    !----------------------------------------------------



    allocate(mapa(L,L))             ! reservo memoria para el mapa
    dim_mapa = size(shape(mapa))    ! verifico la dimension
    allocate(init_pos(dim_mapa))    ! reservo memoria para la posicion inicial dentro del mapa

    !print *,"dim:",dim_mapa," size:",shape(mapa)

    call init_mapa(mapa,L)          ! inicio el mapa 
    !do j=1,L                            
    !    write(*,*) (mapa(i,j), i=1,L)   ! muestro por pantalla el mapa
    !end do

    call plot_mapa(ier,mapa,L)

    call rand_pos(L,dim_mapa,init_pos)      ! pido una posicion aleatoria en el "init_pos"
    !print *,init_pos                    ! muestro por pantalla "init_pos"

    area = L*L
    densidad = 1./real(area)
    !pasos = 1./real(N_MC)

    !
    allocate(ip(L))             ! reservo memoria para suma vecinos
    allocate(im(L))             ! reservo memoria para restar vecinos

    do i=1, L
        ip(i) = i+1
        im(i) = i-1
    end do
    ip(L) = 1
    im(1) = L

    open(3, file='ising.dat',status='unknown')
    open(4, file='cambios.dat',status='unknown')
    T = T0
    do while(T < T0+3.5)
        ut = U/T
        ht = H/T

        xetot = 0.
        xmtot = 0.
        xe2tot = 0.
        xm2tot = 0.

        n_changes = 0
        do imc=1, N_MC
            do i_in=1, N_INTER
                do i_area=1, area
                    i = nint(L*uni()) + 1
                    j = nint(L*uni()) + 1
                    if(i>L)then
                        i=L
                    end if
                    if(j>L)then
                        j=L
                    end if
                    isum = mapa(im(i),j)+mapa(ip(i),j)+mapa(i,im(j))+mapa(i,ip(j))
                    !print *,"posicion",i,j
                    !print *,"suma",isum
                    energia_antigua = -ut*mapa(i,j)*isum - ht*mapa(i,j)
                    energia_nueva   =  ut*mapa(i,j)*isum + ht*mapa(i,j)
                    delta_energ_beta   = energia_nueva - energia_antigua
                    if(delta_energ_beta<0)then
                        mapa(i,j) = -1.*mapa(i,j) ! invertimos
                        n_changes = n_changes + 1 ! cuenta el cambio
!                        print *,"invierte",i,j

                        !para plotear modificando color
                        if( mapa(i,j) == 1 )then
                            call pgsci(2) ! rojo
                        else
                            call pgsci(1) ! blanco
                        end if
                        call pgcirc(real(i),real(j),0.25)
                    else
                        prob_acept = exp(-delta_energ_beta)
                        r = uni()
                        if( r < prob_acept )then

                            mapa(i,j) = -1.*mapa(i,j) ! invertimos
                            n_changes = n_changes + 1 ! cuenta el cambio
                            !print *,"invierte",i,j

                            !para plotear modificando color
                            if( mapa(i,j) == 1 )then
                                call pgsci(2) ! rojo
                            else
                                call pgsci(1) ! blanco
                            end if
                            call pgcirc(real(i),real(j),0.25)
                        end if
                    end if
                    !call plot_mapa(ier,mapa,L)
                end do
            end do
            
            sumep = 0.
            summp = 0.

            do j=1,L
                do i=1,L
                    isum = mapa(im(i),j)+mapa(ip(i),j)+mapa(i,im(j))+mapa(i,ip(j))
                    sumep = sumep - U*mapa(i,j)*isum - H*mapa(i,j)
                    summp = summp + mapa(i,j)
                end do
            end do

            sumep = 0.5*sumep*densidad
            xetot = xetot + sumep
            xe2tot = xe2tot + sumep*sumep

            summp = summp*densidad
            xmtot = xmtot + abs(summp)
            xm2tot = xm2tot + summp*summp

            if(mod(imc,1)==0)then
                print *, imc, sumep, xetot/real(imc), summp, xmtot/real(imc)
            end if

        end do
        !open(2, file='conf',status='unknown')
        !do i=1,L
        !    write(2,*) (mapa(i,j),j=1,L)
        !end do
        !close(2)
        write(3,*) T, sumep, xetot/real(imc), summp, xmtot/real(imc)
        write(4,*) T, real(n_changes)/real(imc)
        T = T + 0.2
    end do
    close(3)
    close(4)

    ! Escribir la última semilla para continuar con la cadena de numeros aleatorios 
    open(unit=10,file='seed.dat',status='unknown')
    seed = shr3() 
    write(10,*) seed
    close(10)
end program problema_ising

