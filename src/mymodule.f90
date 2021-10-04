MODULE mymodule
    use ziggurat
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: funcion_x, init_mapa, rpos_mapa

CONTAINS

! Subrutina para crear la matriz 
subroutine init_mapa(mapa,n)
    implicit none
    integer j,i
    integer, intent(in)  :: n
    integer, intent(out) :: mapa(n,n)

    do j=1,n
        do i=1,n
            if(uni()<0.5)then
                mapa(i,j)=  1
            else
                mapa(i,j)= -1
            end if
        end do
    end do
    return  
end subroutine init_mapa

! Selecciona una posicion random de la matriz de mapa
subroutine rpos_mapa(mapa,n,pos,npos)
    implicit none
    integer i
    real rvar
    integer, intent(in) :: n, npos  ! n: tamanio del mapa n x n
                                    ! npos: dimension de la matriz (n x n => 2)
    integer, intent(in) :: mapa(n,n)
    integer, intent(out):: pos(0:npos-1)

    i = 0
    do while (i<npos)
        rvar = abs(shr3())
        !if( 0<=rvar .and. rvar<=n )then
        if( rvar<=n )then
            pos(i) = rvar
            i = i + 1
        end if
    end do
    return
end subroutine rpos_mapa



FUNCTION funcion_x( x , y ) RESULT( hits )
      real :: x , y , fdx
      integer :: hits
      
      fdx = (x - 1.)**2
      
      hits = 0
      if( y <= fdx ) then
              hits = 1
      end if

      RETURN
END FUNCTION funcion_x


END MODULE mymodule
