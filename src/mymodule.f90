MODULE mymodule
    use ziggurat
      IMPLICIT NONE

      PRIVATE

      !PUBLIC :: funcion_x, init_mapa, rpos_mapa
      PUBLIC :: funcion_x, init_mapa, rand_pos, plot_mapa

CONTAINS

! Subrutina para mostrar el mapa
! Para llamar esta subrutina se debe haber iniciado el pgplot
! (IN)  ier: referencia de pgbeg()
! (IN)  mapa: matriz unos y ceros que representa el mapa
! (IN)  nsize: tamanio de la matriz (mapa). ej. nsize=3 => 3 x 3
subroutine plot_mapa(ier,mapa,nsize)
    implicit none
    integer j,i
    integer, intent(in) :: ier
    integer, intent(in) :: nsize
    integer, intent(in) :: mapa(nsize,nsize)

    do j=1,nsize
        do i=1,nsize
            if(mapa(i,j)==1) then
                !call pgsfs(1) ! circulo relleno
                !call pgsci(3) ! verde
                call pgsci(2) ! rojo
            else
                !call pgsfs(2) ! circulo con contorno
                call pgsci(1) ! blanco
            end if
            call pgcirc(real(i),real(j),0.25)
        end do
    end do 
    return
end subroutine 

! Subrutina para crear la matriz 
! (IN)  mapa: escribe aleatoriamente en esta matriz unos y ceros
! (OUT) n: es el tamanio del mapa
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
                !mapa(i,j)= 0
                mapa(i,j)= -1
            end if
        end do
    end do
    return  
end subroutine init_mapa

! Subrutina que elige una posicion random dentro de la matriz del mapa
! (IN)  nsize: tamanio de la matriz (mapa). ej. nsize=3 => 3 x 3
! (IN)  ndim:  la dimension de la matriz (mapa). ej. ndim=3 => 3 x 3 x 3
! (OUT) rpos:  retorna una posicion de la matriz elegida aleatoriamente
subroutine rand_pos(nsize,ndim,rpos)
    implicit none
    integer i
    integer rvar
    integer, intent(in) :: nsize, ndim
    integer, intent(out) :: rpos(ndim)

    do i=1,ndim
        rpos(i) = nint(uni()*nsize)+1
        if(rpos(i)>nsize)then
            rpos(i) = nsize
        endif
    end do
    return

end subroutine rand_pos
!subroutine rand_pos(nsize,ndim,rpos)
!    implicit none
!    integer i
!    integer rvar
!    integer, intent(in) :: nsize, ndim
!    integer, intent(out) :: rpos(0:ndim-1)
!
!    i = 0
!    do while (i<ndim)
!        rvar = abs(shr3())
!        !if( 0<rvar .and. rvar<=nsize )then
!        if( rvar<=nsize )then
!            rpos(i) = rvar
!            i = i + 1
!        end if
!    end do
!    return
!
!end subroutine rand_pos
!subroutine rpos_mapa(mapa,n,npos,pos)
!    implicit none
!    integer i
!    real rvar
!    integer, intent(in) :: n, npos  ! n: tamanio del mapa n x n
!                                    ! npos: dimension de la matriz (n x n => 2)
!    integer, intent(in) :: mapa(n,n)
!    integer, intent(out):: pos(0:npos-1)
!
!
!    i = 0
!    do while (i<npos)
!        rvar = abs(shr3())
!        !if( 0<=rvar .and. rvar<=n )then
!        if( rvar<=n )then
!            pos(i) = rvar
!            i = i + 1
!        end if
!    end do
!    return
!end subroutine rpos_mapa


! Ejercicio
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
