module globals
    real(kind=8), allocatable :: r(:,:) ! Posicion (3,n)
    real(kind=8), allocatable :: velocidad(:,:) ! velocidad (3,n)
    real(kind=8), allocatable :: f(:,:) ! fuerza (3,n)
    real(kind=8)              :: V ! Potencial
	real(kind=8)              :: Temperatura
	real(kind=8)              :: L,Ec,Etotal
    real(kind=8), allocatable :: delta_R(:)
    integer                   :: N
end module globals
