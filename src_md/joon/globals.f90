module globals
    integer :: n ! Numero de atomos
    integer :: c ! coordenadas
    real(kind=4) :: L, Vtot, T
    real(kind=4), allocatable :: r(:,:) ! Posicion (3,n)
    real(kind=4), allocatable :: v(:,:) ! velocidad (3,n)
    real(kind=4), allocatable :: f(:,:) ! fuerza (3,n)
end module globals
