MODULE mymodule
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: funcion_x

CONTAINS

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
