program gausiana_1
	use ziggurat
	implicit none
	logical :: es
	integer :: seed, i, j, n, cont_g ,cont_ac,nbins
	real:: x,y,a,b,suma,centro, delta,k, maximo, minimo
	real, parameter:: pi = 3.1415927 
	real (kind=8), allocatable :: x1(:), y1(:),hist(:),hist_n(:),limites(:)
	 
![NO TOCAR] Inicializa generador de número random

	inquire(file='seed.dat',exist=es)
	if(es) then
		open(unit=10,file='seed.dat',status='old')
		read(10,*) seed
		close(10)
		print *,"  * Leyendo semilla de archivo seed.dat"
	else
		seed = 24583490
	end if

	call zigset(seed)
![FIN NO TOCAR]    


open(1, file = 'input.dat',status='old')
open(2, file = 'output.dat',status='old')
open(3, file = 'distribucion.dat', status='old') 
open(4, file=  'histograma.dat', status='old')
! Ej: Número random en [0,1]: uni()

read(1,*) n
close(1)

cont_g=0
cont_ac=0


allocate(x1(n),y1(n))
k=1.0
	   do i = 1, n
	       x=uni() 
	       y=uni()
	       if (x<0.5) then 
		   x=rnor() 
		   else 
		   x=rnor()+2
		   end if 
		     
		   	   
		   if (x>=0.and.x<=2) then 
		   cont_g=cont_g+1
		   write(3,*) i, x, y  
		   
		   b=(x-1-0)**2.0
		   b=k*b
		   b=b/exp(-x**2.0/2.0)
		   b=b*(2*Pi)**0.5
		   !b=b*2
		   if (b>=y) then 
			  cont_ac=cont_ac+1
      		  !print *,i,x,y,cont               
			  x1(i)=x
			  y1(i)=y
			  write(2,*) i, x, y, cont_g,cont_ac
			  !write(2,fmt="(1i5,x,2f16.5,x,1i5,x,a)") i, b, x, y, cont  
		   end if   
		   end if 
		end do

!!!!HISTOGRANA HACER COMO FUNCION
!!! que calcule los bines
!!! modificar para que haga con nbins no con delta 
 
delta=0.1
maximo=2.0
minimo=0
  

nbins=nint((maximo-minimo)/delta)+1 ! lo tengo que carga

!print *, nbins   
   

allocate(hist(nbins-1),hist_n(nbins-1),limites(nbins))

limites(1)=minimo  

   
 
  do i=2,nbins
     limites(i)=limites(i-1)+delta
     !print *,limites(i)
   end do 

hist(:)=0
       do i=1,n 
         do j=1,nbins-1 
             if (x1(i)>limites(j).and.x1(i)<=limites(j+1)) then  
                hist(j)=hist(j)+1 
               end if 
          end do 
       end do        
      
hist(:)=hist(:)/dble(cont_ac)      

b=limites(2)-limites(1)
b=b/2.0
     
suma=0
print *," "
do i=1,nbins-1
     suma=suma+hist(i)
     centro=limites(i)+b 
     print *,"Bin : ", i, " Xo: ", limites(i), "X1: ", limites(i+1)   ,"Centro : ", centro,  "Valor  : " ,hist(i)
     write(4,*) centro, hist(i)
end do
print *, "  " 
print *, " Suma : " , suma  

a=dble(cont_g)/dble(n)
b=dble(cont_ac)/dble(cont_g)

print *, "    "   
print *, "El numero N es:", n
print *, "El numero de Generados sampliando una gausiana:", cont_g
print *, "La proporcion Generados/N:", a
print *, "Numero de aceptados:", cont_ac
print *, "La proporcion Aceptados/generasdos:", b
print *, "La integral es:", b/2.0
print *, "La integral teorica es 1/3"

close(2)
close(3)
!! Llamo a multiplico
!a=100.0 
!b=500.0
!call multiplico(a,b,z) 
!entrada x,y integer
!devuelve z
!print *,z

!! 
!! EDITAR AQUI 
!! 


!! 
!! FIN FIN edicion
!! 
![No TOCAR]
! Escribir la última semilla para continuar con la cadena de numeros aleatorios 

		open(unit=10,file='seed.dat',status='unknown')
		seed = shr3() 
		 write(10,*) seed
		close(10)
![FIN no Tocar]        

end program gausiana_1
