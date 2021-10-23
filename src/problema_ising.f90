program problema_ising
    ! tareas y funciones de otro archivo
    use ziggurat
    use mymodule
    implicit none

    ! Variables
    logical :: es
    integer :: seed, i , j, icolor
    integer :: L! tamanio de la matriz. ej. L x L
    integer , allocatable :: mapa(:,:)
    integer , allocatable :: init_pos(:)
    integer :: dim_mapa

    integer :: ier,pgbeg ! para referenciar al inicio de pgplot

    integer, allocatable :: ip(:),im(:) ! para ir mapeando cuales son los vecinos
    real :: N0, Jta, T0, T, H
    integer :: N_MC, N_INTER
    integer :: imc, i_in
    integer :: i_area, area
    integer :: isum
    real :: energia_antigua, energia_nueva, delta_energ_beta, ising
    real :: sumep, summp, xetot, xe2tot, xmtot, xm2tot
    real :: xmagnet, xenergy, xmag2, xener2
    real :: var_e, var_m, suscept, cv

    real :: ht, beta
    real :: r, prob_acept
    real :: densidad, pasos_mc
    
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

    !----------------------------------------------------

    ! Inicio lectura de input para el ising
    ! N0: nro de spin
    ! Jta: Energia de interaccion
    ! T: Temperatura
    ! H: Campo magnetico
    ! N_MC: nro de iteraciones de Monte Carlo
    inquire(file='input_ising.dat',exist=es)
    if(es) then
        open(unit=1,file='input_ising.dat',status='old')
        read(1,*)
        read(1,*) N0 ,Jta ,T0, H
        read(1,*)
        read(1,*) N_MC, N_INTER, L
        close(1)
!        print *,"  * leyendo datos de input_ising.dat"
    else
        N0=0
        Jta=1
        T0=1
        H=0
        N_MC=1000
    end if
    close(1)

    !----------------------------------------------------



    allocate(mapa(L,L))             ! reservo memoria para el mapa
    dim_mapa = size(shape(mapa))    ! verifico la dimension
    allocate(init_pos(dim_mapa))    ! reservo memoria para la posicion inicial dentro del mapa


    call init_mapa(mapa,L)          ! inicio el mapa 
    !do j=1,L                            
    !    write(*,*) (mapa(i,j), i=1,L)   ! muestro por pantalla el mapa
    !end do
    !do j=1,L
    !    do i=1,L
    !        mapa(i,j)=1 ! inicio el mapa en 1
    !    end do
    !end do


    !call rand_pos(L,dim_mapa,init_pos)      ! pido una posicion aleatoria en el "init_pos"
    !print *,init_pos                    ! muestro por pantalla "init_pos"

    area = L*L
    densidad = 1./real(area)    ! densidad de area
    pasos_mc = 1./real(N_MC)       ! densidad o pasos de monte carlo

    allocate(ip(L))             ! reservo memoria para suma vecinos
    allocate(im(L))             ! reservo memoria para restar vecinos

    ! me va a ayudar para el calculo de condicion de contorno
    do i=1, L
        ip(i) = i+1
        im(i) = i-1
    end do
    ip(L) = 1
    im(1) = L

    open(3, file='ising.dat',status='unknown')
    open(4, file='fraccion_aceptados.dat',status='unknown')
    open(5, file='caloresp_suscept.dat',status='unknown')
    T = T0
    do while(T < T0+4)
        beta = Jta/T    ! 1/K.T
        ht = H/T

        xenergy = 0
        xener2 = 0
        xmagnet = 0
        xmag2 = 0

        n_changes = 0

        do i_in=1, N_INTER
            xetot = 0.
            xmtot = 0.
            xe2tot = 0.
            xm2tot = 0.

            do imc=1, N_MC
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
                    energia_antigua = -beta*mapa(i,j)*isum - ht*mapa(i,j)
                    energia_nueva   =  beta*mapa(i,j)*isum + ht*mapa(i,j)
                    delta_energ_beta   = energia_nueva - energia_antigua
                    if(delta_energ_beta<0)then
                        mapa(i,j) = -1.*mapa(i,j) ! invertimos
                        n_changes = n_changes + 1 ! cuenta el cambio
!                        print *,"invierte",i,j
                    else
                        prob_acept = exp(-delta_energ_beta)
                        r = uni()
                        if( r < prob_acept )then

                            mapa(i,j) = -1.*mapa(i,j) ! invertimos
                            n_changes = n_changes + 1 ! cuenta el cambio
                            !print *,"invierte",i,j

                        end if
                    end if
                end do ! area

                sumep = 0.
                summp = 0.

                ! Calculo de energia y magnetizacion
                do j=1,L
                    do i=1,L
                        isum = mapa(im(i),j)+mapa(ip(i),j)+mapa(i,im(j))+mapa(i,ip(j))
                        sumep = sumep - Jta*mapa(i,j)*isum - H*mapa(i,j)
                        summp = summp + mapa(i,j)
                    end do
                end do

                sumep = sumep*densidad          ! Energia interna
                xetot = xetot + sumep           ! Para calcular la media de la energia
                xe2tot = xe2tot + sumep**2      ! Para calcular la media de la energia cuadrada

                summp = summp*densidad          ! Magnetizacion media por espin
                xmtot = xmtot + abs(summp)      ! Para calcular la media de la magnetizacion
                xm2tot = xm2tot + summp**2      ! Para calcular la media de la magnetizacion cuadrada

                !if(mod(imc,1)==0)then
                !    print *, imc, sumep, xetot/real(imc), summp, xmtot/real(imc)
                !end if

            end do !N_MC

            if(i_in>1)then
                xenergy = xenergy + (xetot / real(imc)) ! la media de la energia por spin
                xener2 = xener2 + (xe2tot / real(imc)) ! la media de la energia cuadrada
        
                var_e = var_e + (xe2tot - xetot**2)/real(imc)

                xmagnet = xmagnet + (xmtot / real(imc)) ! la media del magnetismo por spin
                xmag2 =  xmag2 + (xm2tot / real(imc)) ! la media del magnetismo al cuadrado

                var_m = var_m + (xm2tot - xmtot**2)/real(imc)

            end if
        end do ! N_INTER

        !cv = var(E)/Kb.N.T^2
        var_e = (xener2/real(N_INTER-1)-(xenergy/real(N_INTER-1))**2)
        cv = ((var_e * beta) / (T))  ! Calor Especifico 

        !suscrp = N.var(M)/Kb.T
        var_m = (xmag2/real(N_INTER-1)-(xmagnet/real(N_INTER-1))**2)
        suscept = (var_m * real(imc) * beta)    ! Susceptibilidad magnetica

        write(3,*) T, xenergy/real(N_INTER-1), xmagnet/real(N_INTER-1), xener2/real(N_INTER-1), xmag2/real(N_INTER-1)
        write(4,*) T, real(n_changes)/real(imc*(N_INTER-1)*area)
        write(5,*) T, var_e , var_m, cv, suscept
        T = T + 0.01
    
    end do ! temperatura

    close(3)
    close(4)
    close(5)

    ! Escribir la última semilla para continuar con la cadena de numeros aleatorios 
    open(unit=10,file='seed.dat',status='unknown')
    seed = shr3() 
    write(10,*) seed
    close(10)
end program problema_ising

