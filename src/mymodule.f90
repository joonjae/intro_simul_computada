MODULE mymodule
      IMPLICIT NONE

      PRIVATE

      PUBLIC :: funcion_x

CONTAINS

FUNCTION funcion_x( x , y , filenumb ) RESULT( hits )
      real :: x , y , fdx
      integer :: hits, filenumb
      
      fdx = (x - 1.)**2
      
      hits = 0
      if( y <= fdx ) then
              write(filenumb ,*) x , y
              hits = 1
      end if

      print *,"hola mundo"
      RETURN
END FUNCTION funcion_x

END MODULE mymodule