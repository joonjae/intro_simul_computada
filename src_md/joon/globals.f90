module globals
    integer :: n ! Numero de atomos
    integer :: c ! coordenadas
    real(kind=8) :: L, Vtot, T
    real(kind=8), allocatable :: r(:,:) ! Posicion (3,n)
    real(kind=8), allocatable :: v(:,:) ! velocidad (3,n)
    real(kind=8), allocatable :: f(:,:) ! fuerza (3,n)
end module globals
