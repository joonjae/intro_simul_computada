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
end subroutine verlet_positiones


! Subrutina para calculo de energia potencial y fuerza de las particulas
! (IN)  sigma: medida del diametro efectivo de las particulas
! (IN)  espsl: epsilon, medida del pozo potencial (intensidad de la interaccion)
! (OUT) 
! (GLOBAL)  L: Longitud de la caja
!           R: posicion de las particulas
!           V: energia potencial de las particulas
!           F: fuerza total sobre cada particula
subroutine fuerza(sigma,epsl)
    implicit none
    real, intent(in) :: sigma
    real, intent(in) :: epsl

    integer :: j,i
    real    :: sr2, sr6, sr12, pot_corte
    real, dimension(3) :: rij, fij

    r_corte_caja = r_corte / caja

    F = 0.
    do i=1,n - 1
        do j=i+1,n

            d(:) = r(:,i) - r(:,j) ! distancia entre par de particulas
            d_sq = SUM ( d**2 )
            sr2 = (sigma**2) / d_sq

            sr6 = sr2**3
            sr12= sr6**2

            v_r = 4*epsil * (sr12 - sr6)

            v_sist = v_sist + v_r


        end do
    end do

    return
end subroutine fuerza


END MODULE mymodule
