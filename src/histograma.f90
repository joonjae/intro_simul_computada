SUBROUTINE histograma(x,n,minimo,maximo,nbins,hist)

 implicit none 

   integer :: i,j,cont  
   integer, intent(in) :: nbins,n
   real, intent(in):: x,minimo,maximo
   real :: b, centro, delta, suma    
   real (kind=8), allocatable :: hist(:),limites(:)
   

 
open(900, file=  'histograma.dat', status='old')

delta=(maximo-minimo)/nbins


allocate(hist(nbins-1),limites(nbins))

limites(1)=minimo  

   do i=2,nbins
     limites(i)=limites(i-1)+delta
     !print *,limites(i)
   end do 

hist(:)=0
cont=0

       do i=1,n 
         do j=1,nbins-1 
             if (x(i)>limites(j).and.x(i)<=limites(j+1)) then  
                hist(j)=hist(j)+1 
                cont=cont+1
               end if 
          end do 
       end do        
      
hist(:)=hist(:)/dble(cont)      

b=limites(2)-limites(1)
b=b/2.0
     
suma=0
!print *," "
do i=1,nbins-1
     suma=suma+hist(i)
     centro=limites(i)+b 
!     print *,"Bin : ", i, " Xo: ", limites(i), "X1: ", limites(i+1)   ,"Centro : ", centro,  "Valor  : " ,hist(i)
     write(900,*) centro, hist(i)
end do
!print *, "  " 
!print *, " Suma : " , suma  
close(900)

END SUBROUTINE 
