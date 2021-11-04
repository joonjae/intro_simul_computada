MODULE mymodule
    use ziggurat
    use globals
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: init_posiciones

CONTAINS

! Subrutina para crear la matriz 
! (IN)  n_dimensiones: 
! (IN)  n_particulas: 
! (OUT) mapa: 
!subroutine init_posiciones(mapa,dimensiones,n_particulas)
subroutine init_posiciones(n_dimensiones,n_particulas)
    implicit none
    integer j,i
    integer, intent(in)  :: n_particulas
    integer, intent(in)  :: n_dimensiones
!    integer, intent(out) :: mapa(n,n)

    do j=1,n_dimensiones
        do i=1,n_particulas
            R(j,i) = uni()*L
        end do
    end do
    return  
end subroutine init_posiciones



END MODULE mymodule
