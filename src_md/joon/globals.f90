module globals
    real(kind=8), allocatable :: R(:,:), V(:,:), F(:,:) ! para representar posiciones, velocidades y fuerzas de N particulas
    real(kind=8) :: L, Vtot, T
    integer :: N, D
end module globals
