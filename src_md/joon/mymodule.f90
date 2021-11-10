MODULE mymodule
    use ziggurat
    use globals
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: verlet_posiciones, fuerza

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
            R(j,i) = uni()*L
            print *, j,i," - posicion: ",R(j,i)
        end do
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
    real, dimension(3) :: d

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


END MODULE mymodule
