MODULE mymodule
    use ziggurat
    use globals
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: verlet_posiciones, fuerza, write_conf

CONTAINS

! Subrutina para inicio de las posiciones random
! (IN)  NA
! (OUT) NA
! (GLOBAL)  L: Longitud de la caja
!           r: posicion de las particulas
subroutine verlet_posiciones()
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
end subroutine verlet_posiciones


! Subrutina para calculo de energia potencial y fuerza de las particulas
! (IN)  sigma: medida del diametro efectivo de las particulas
! (IN)  espsl: epsilon, medida del pozo potencial (intensidad de la interaccion)
! (OUT) 
! (GLOBAL)  L: Longitud de la caja
!           r: posicion de las particulas
!           f: fuerza total sobre cada particula
subroutine fuerza(sigma,epsl)
    implicit none
    real, intent(in) :: sigma
    real, intent(in) :: epsl

    integer :: j,i
    real    :: d_sq, sr2, sr6, sr12
    real    :: v_r, v_sist
    real, dimension(c) :: d

    f = 0.
    do i=1,n - 1 ! particula i
        do j=i+1,n ! particula j

            d(:) = r(:,i) - r(:,j) ! distancia entre par de particulas
            d_sq = SUM ( d**2 )
            sr2 = (sigma**2) / d_sq

            sr6 = sr2**3
            sr12= sr6**2

            v_r = 4*epsl * (sr12 - sr6)

            v_sist = v_sist + v_r

            f(:,i) = f(:,i) + d * v_r
            f(:,j) = f(:,j) - d * v_r

        end do
    end do

    return
end subroutine fuerza

subroutine write_conf(mode)
    implicit none
    integer, intent (in)::mode
    integer             ::i
    
    select case(mode)
        case(0)
            open(unit=20,file="movie.vtf",status="unknown")
            write(20,*) "atom 0:99 radius 0.5 name Ar"
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
