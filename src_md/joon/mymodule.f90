MODULE mymodule
    use ziggurat
    use globals
    IMPLICIT NONE

    !------------------------------------------------------------------------------
    integer npart
    real, allocatable :: pos(:,:),vel(:,:),acc(:,:),mass(:),mass_inv(:),pos_t0(:,:)
    real pot, energy
    real dt,total_time,curr_time, reset_temp_time, eq_time, final_time
    integer nsteps, nstep_reset, nstep_eq
    real epsil,sigma,sigma6,eps4,eps48, pi, cons
    real density,temperature,box_len, temp_in
    real E_avg,Esq_avg,E_std_dev
    real V_C, sigma6_rc6, rc_2, msd
    integer flag_write, flag_vmd
    !------------------------------------------------------------------------------

    PRIVATE

    PUBLIC :: init_posiciones, fuerza, write_conf
    PUBLIC :: setup_parameters, initial_condition, equilibrate, evolve
    PUBLIC :: velocity_verlet

CONTAINS

!------------------------------------------------------------------------------
subroutine distance(i,j,rijsq,vec)
    !! input: i and j, that represent atoms i and j
    !! Output: vec(3) = pos(i,:)-pos(j,:) (vector from atom j to i)
    !! Output: rijsq: square of the distance between atoms i and j
    implicit none
    integer,intent(in) :: i,j
    real,intent(out) :: rijsq,vec(3)

    integer k

    do k=1,3
        vec(k)=pos(i,k)-pos(j,k)

        if(vec(k) > 0.5*box_len) vec(k)= vec(k)-box_len
        if(vec(k) < -0.5*box_len) vec(k)= vec(k)+box_len
    end do

    rijsq = (sum(vec*vec))

end subroutine distance

subroutine LJ_potential(xsq,V_LJ,dVLJ_dx)
    !! Takes xsq=x_squares as the input
    !! Calculates V_LJ=LJ potential and dVLJ_dx=(1/x).derivative of potential
    implicit none
    real,intent(in)::xsq
    real,intent(out)::V_LJ,dVLJ_dx
    real sig6_x6

    sig6_x6=sigma6/xsq**3

    V_LJ=(eps4*sig6_x6*(sig6_x6-1.)) - V_C
    dVLJ_dx=eps48/xsq*sig6_x6*(0.5-sig6_x6)

end subroutine LJ_potential


subroutine compute_pot
        implicit none
        integer i,j
        real rijsq,vec(3),V_LJ,dVLJ_drij
        real dpot_dx(3)

        pot=0.
        acc=0.
        do i=1,npart-1
            do j=i+1,npart
                call distance(i,j,rijsq,vec)
                if(rijsq<rc_2) then
                    call LJ_potential(rijsq,V_LJ,dVLJ_drij)
                    pot=pot+V_LJ
                    dpot_dx=dVLJ_drij*vec
                    acc(i,:)=acc(i,:)-mass_inv(i) * (dpot_dx)
                    acc(j,:)=acc(j,:)+mass_inv(j) * (dpot_dx)
                end if
            enddo
        enddo

end subroutine compute_pot

subroutine velocity_verlet
    implicit none
    
    pos = pos + vel*dt + 0.5*acc*(dt**2)
    vel = vel + 0.5*acc*dt
    call compute_pot
    vel = vel + 0.5*acc*dt

end subroutine velocity_verlet



subroutine compute_energy
    implicit none
    integer i

    energy=pot
    do i=1,npart
      energy=energy+0.5*mass(i)*sum(vel(i,:)*vel(i,:))
    enddo

    !E_avg=E_avg+energy
    !Esq_avg=Esq_avg+energy**2

end subroutine compute_energy


subroutine get_inside_the_box
    implicit none
    integer i,j

    do i=1,npart                              !! Loop over all particles
       do j=1,3                               !! Loop along x,y,z direction

          if (pos(i,j)>box_len) then
            pos(i,j) = pos(i,j)-box_len
          else if (pos(i,j) <0) then
            pos(i,j) = pos(i,j)+box_len
          end if

       end do
    end do

end subroutine get_inside_the_box


subroutine thermal_velocities(temperature)
     implicit none
     !double precision,intent(in) :: temperature
     real ,intent(in) :: temperature
     !double precision rnd,kb
     real rnd,kb
     integer i,j

     kb=1.

     do i=1,npart
        do j=1,3
            !call gaussian_random_number(rnd)            !! Random number with gaussian distribution with mean=0, sigma=1
            rnd = rnor()
            vel(i,j)=rnd*sqrt(kb*temperature/mass(i))   !! Random number with gaussian distribution with mean=0, sigma=sqrt(kb*temperature/mass(i))
        enddo
     enddo

end subroutine thermal_velocities



subroutine write_output_msd
    implicit none
    write(10,*) curr_time,msd
end subroutine write_output_msd

subroutine write_energy
    implicit none
    write(11,*) curr_time,energy
end subroutine write_energy

subroutine write_temp
    implicit none
    write(18,*) curr_time, temp_in
end subroutine write_temp

subroutine write_vmd
    implicit none
    integer i

    write(13,*) npart
    write(13,*)
    do i=1,npart
      write(13,*) "C",pos(i,:)
    enddo

end subroutine write_vmd

subroutine equilibrate
    implicit none
    integer i

    nstep_eq=nint(eq_time/dt)
    nstep_reset=nint(reset_temp_time/dt)

    open(18,file="temp_time.out")
    open(11,file="energy.out")
    open(13,file="vmd.vtf")

    do i=1,nstep_eq
        call compute_energy
        if(mod(i,nstep_reset)==1) then
            call thermal_velocities(temperature)
        end if
        temp_in = cons*(energy-pot)
        if(flag_write==1 .and. mod(i,50)==1) then
            call write_temp
        end if
        if(flag_write==1 .and. mod(i,50)==1) then
            call write_energy
        end if
        if(flag_vmd==1 .and. mod(i,100)==1) then
            call write_vmd
        end if
        call velocity_verlet
        call get_inside_the_box
        curr_time=curr_time+dt
    enddo

    close(18); close(11); close(13)

    pos_t0 = pos

end subroutine equilibrate


subroutine setup_parameters
    implicit none
    real mass_atom

    open(10,file="md_LJ_fluid.inp")
    read(10,*) npart
    read(10,*) dt
    read(10,*) reset_temp_time
    read(10,*) eq_time
    read(10,*) final_time
    read(10,*) mass_atom
    read(10,*) sigma
    read(10,*) epsil
    read(10,*) density
    read(10,*) temperature
    read(10,*) flag_write
    read(10,*) flag_vmd
    close(10)

    allocate(pos(npart,3),vel(npart,3),acc(npart,3),mass(npart),mass_inv(npart), pos_t0(npart,3))
    mass=mass_atom
    mass_inv=1./mass

    sigma6=sigma**6
    eps4=epsil*4
    eps48=epsil*48
    rc_2= 6.25*(sigma**2)
    sigma6_rc6 = sigma6/rc_2**3
    V_C = eps4*sigma6_rc6*(sigma6_rc6 - 1.)
    pi = 3.141592653
    cons = 2./(3.*npart)

end subroutine setup_parameters


subroutine initial_condition
    implicit none
    integer i,j,k
    integer npart_side,part_num
    real lattice_len

    npart_side=nint(1.*npart**(1/3.))   !! no. of particles along each direction
    box_len=(npart/density)**(1/3.)       !! length of simulation box
    lattice_len = box_len/real(npart_side)  !! length of the cubic unit cell

    write(*,*) "Initial conditions: placing particles in simple cubic lattice:"
    write(*,*) "total number of particles = ",npart
    write(*,*) "particles along each side = ",npart_side
    write(*,*) "box length = ",box_len
    write(*,*) "lattice length = ",lattice_len

    pos=0.
    !vel=0.
    call thermal_velocities(temperature)     !! initial velocity

    !! putting particles on a simple cubic lattice
    part_num=1
    do i=1,npart_side
        do j=1,npart_side
            do k=1,npart_side
                pos(part_num,1)=(i-1)*lattice_len
                pos(part_num,2)=(j-1)*lattice_len
                pos(part_num,3)=(k-1)*lattice_len
                part_num=part_num+1
            enddo
        enddo
    enddo

    !! Writing initial conditions in a vmd compatible xyz file
    open(13,file="init_cond.xyz")
    call write_vmd
    close(13)

    call compute_pot
    call compute_energy
    curr_time=0.

    E_avg=0.
    Esq_avg=0.

end subroutine initial_condition

subroutine compute_msd
    implicit none

    integer i, k
    real*8 r_squared, vect(3)

    msd=0.
    do i=1,npart
        vect=pos(i,:)-pos_t0(i,:)
        do k=1,3
            if(vect(k) > 0.5*box_len) vect(k)= vect(k)-box_len
            if(vect(k) < -0.5*box_len) vect(k)= vect(k)+box_len
        end do
        r_squared=sum(vect*vect)
        msd=msd+r_squared
    end do
    msd=msd/npart

end subroutine compute_msd

subroutine evolve
    implicit none
    integer i

    nsteps=nint(final_time/dt)
    curr_time = 0.
    open(10,file="msd_vs_time.out")

    do i=1,nsteps
        call compute_energy
        if(flag_write==1 .and. mod(i,5)==1) call write_output_msd
        call compute_msd
        call velocity_verlet
        call get_inside_the_box
        curr_time=curr_time+dt
    enddo
    !E_avg=E_avg/real(nsteps)
    !Esq_avg=Esq_avg/real(nsteps)
    !E_std_dev = dsqrt(Esq_avg-E_avg**2)
    !write(*,*) "dt (a.u.) standard deviation of E (a.u)",dt, E_std_dev
    close(10)

end subroutine evolve



! Subrutina para inicio de las posiciones random
! (IN)  NA
! (OUT) NA
! (GLOBAL)  L: Longitud de la caja
!           r: posicion de las particulas
subroutine init_posiciones()
    implicit none
    integer j,i

    do i=1,n
        do j=1,c
            r(j,i) = uni()*L
            print *, j,i," - posicion: ",r(j,i)
        end do
        print *,''
    end do
    return  
end subroutine init_posiciones


! Subrutina para calculo de energia potencial y fuerza de las particulas
! (IN)  sigma: medida del diametro efectivo de las particulas
! (IN)  espsl: epsilon, medida del pozo potencial (intensidad de la interaccion)
! (OUT) 
! (GLOBAL)  L: Longitud de la caja
!           r: posicion de las particulas
!           f: fuerza total sobre cada particula
subroutine fuerza(sigma,epsl,r_corte,v_sist)
    implicit none
    real, intent(in) :: sigma
    real, intent(in) :: epsl
    real, intent(in) :: r_corte
    real, intent(out):: v_sist

    integer :: j,i
    real    :: d_sq, sr2, sr6, sr12
    real    :: v_r
    real    :: caja_sq, d_corte_caja_sq
    real, dimension(c) :: d
    real, dimension(c) :: fij


    !d_corte_caja_sq = (r_corte / L)**2 ! esto me va a ser util para comparar con la condicion pediodica de contorno
    d_corte_caja_sq = (2*sigma)**2
    caja_sq = L ** 2

    f = 0.
    do i=1,n - 1 ! particula i
        do j=i+1,n ! particula j

            d(:) = r(:,i) - r(:,j) ! distancia entre par de particulas
            d(:) = d(:) - L * INT( 2*d(:)/L ) ! Condiciones periodicas de contorno
            d_sq = SUM ( d**2 )

            !if( d_sq < d_corte_caja_sq )then
            if( d_sq < d_corte_caja_sq )then
                !d_sq = d_sq * caja_sq
                !d(:) = d(:) * L
                sr2 = (sigma**2) / sqrt(d_sq)
                sr6 = sr2**3
                sr12= sr6**2

                v_r = (sr12 - sr6)

                v_sist = v_sist + v_r

                ! f(r) = -DELTA(V_r) = -(d/dr) V_r * r_vector = 4*epsilon*(12*sigma^12/r^13 - 6*sigma^6/r^7) * r_vector
                ! f(r) = 4*epsilon*6*(2*sigma^12/r^12 - sigma^6/r^2) * r_vector/r
                ! f(r) = 24*epsilon * (2*sigma^12/r^12 - sigma^6/r^2) * r_vector/r
                fij = (sr12 + v_r) * d/sqrt(d_sq)
                f(:,i) = f(:,i) + fij
                f(:,j) = f(:,j) - fij
            end if
        end do
    end do

    v_sist = 4*epsl * v_sist
    f = 24*epsl * f

    return
end subroutine fuerza

subroutine write_conf(mode)
    implicit none
    integer, intent (in)::mode
    integer             ::i
    
    select case(mode)
        case(0)
            open(unit=20,file="movie.vtf",status="unknown")
            write(20,*) "atom 0:99 radius 1 name Ar"
        case(1)
            write(20,*) "timestep"
            write(20,*)
            do i=1,n
                write(20,*) r(:,i)
            end do
        case(2)
            close(20)
    end select

end subroutine write_conf

END MODULE mymodule
