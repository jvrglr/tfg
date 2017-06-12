        program DESORDEN_ESPACIAL
          implicit none
          integer*4 k,cruz,i,a,j,m !k y cruz guardan los numeros aleatorios/ a,j,i=indices
          !m=controla que solo mida en cada MC
!/////////////////////////////////////////////////////
          integer*4 lt,n,n_a,n_b,tmax,l,celdas,r,x,y,especie
          integer*4 x_vec,y_vec,coexistencia
          parameter (l=5,celdas=25)
          !n=sitios activos en total
          !n_a=sitios activos de especie a estado(i,j,o)=1
          !n_b=sitios activos de especie b estado(i,j,o)=2
          !especie guarda la especie del nodo que crea descendencia
          !tmax=tiempo maximo de ejecucion
          !lt="tiempo logaritmico"
          !l= tamaño de la red 1-D
          !celdas=lxl= nº total de nodos
          !r=numero de realizaciones
          !x,y,x_vec e y_vec serviran para guardar coordenadas de los nodos
          !especie= auxiliar que uso para guardar estado (x,y,0) (tipo de nodo)
          !coexistencia=veces que la simulación acaba con coexistencia (ni n_a ni n_b=0)
!/////////////////////////////////////////////////////////////////////7
          real*8 p,t
          !p=probabilidad de que una modificación sea una creacion
          !t=tiempo
          real*8 realran !Auxiliares
          real ran_u
          integer*4 i_ran
          integer*4 lista(celdas)
          !lista = (x1,x2,...,xk,...,xn,0,...,0)
          !Coordenadas del nodo k-ésimo ocupado: y=resto(xk/l)+1 , x=xk/l
          real*8 mean_n_a(500),mean_n_b(500),mean_p(500)
          integer*4 eventos(500),eventos_total(500)
          integer*4 estado(l,l,0:5)
          !estado (x,y,0)= estado del nodo (x,y)
          !estado (x,y,5)= preferencia del nodo (x,y) DESORDEN ESPACIAL
          !estado (x,y,k) (k!=0)= coordenada del vecino k-ésimo del nodo (x,y)
          real*8 mean_t_debil,mean_t_fuerte
          !mean_t_debil=tiempo medio que sobrevive la especie que se extingue
          !mean_t_fuerte= tiempo medio que sobrevive la especie que no se extingue (satura a tmax)
          integer aux !auxiliar que controla si hay o no coexistencia
          !aux =1 while hay coexistencia
          !PARAMETROS DE DESORDEN ESPACIAL
          real*8 epsilon, p0
          !epsilon=incremento de probabilidad debido a preferencia espacial
          !p0 y control uxiliares
          integer*4 control

          !Seed de numeros aleatorios, abrimos fichero
          !call ran_ini(time())
          call ran_ini(104723)
      open(unit=1, file="datos.dat", status="unknown",position="append")
          !CARACTERISTICAS DEL SISTEMA
          r=10000
          tmax=50000
          p0=0.7
          !INICIALIZACIÓN DE VARIABLES
          mean_n_a=0.0
          mean_n_b=0.0
          mean_p=0.0
          mean_t_fuerte=0.0
          mean_t_debil=0.0
          eventos=0
          eventos_total=0
          epsilon=0.01
          coexistencia=0


          !/-------------------_/
          !Relleno vecinos en tensor estado (i,j,k)
          do i=1,l
            do j=1,l
              estado(i,j,1) =i+j*l !vecino en la fila superior
              estado(i,j,2) =i-2+j*l !vecino en la fila inferior
              estado(i,j,3) =i-1+(j+1)*l !vecino en la columna de la derecha
              estado(i,j,4) =i-1+(j-1)*l !vecino en la columna de la izq.
            enddo
          enddo


          !Condiciones de contorno periódicas
          do i=1,l
            estado(1,i,2)=l-1+i*l !la  primera fila une por abajo con la última
            estado(l,i,1)=i*l !la última fila une por arriba con la primera
            estado(i,1,4)=i-1+l*l !la primera columna unida por la izquierda con la última
            estado(i,l,3)=i-1+l !la última columna unida por la derecha con la primera
          !/--------------------/
          enddo

            do j=1,r
              !CONDICIONES INICIALES homogéneas
              aux=1 !hay coexistencia
              t=1.0
              i=0
              n_a=0
              n_b=0
              do x=1,l
                do y=1,l
                  i=i+1
                  especie=i_ran(2) !especie que puebla el nodo (j,i) en t=0
                  estado(x,y,0)=especie
                  estado(x,y,5)=i_ran(2) !preferencia del nodo (j,i) ¡no dinámico!
                  if (especie.eq.1) then
                    n_a=n_a+1
                  else
                    n_b=n_b+1
                  endif
                  lista(i)=y-1+l*x
                enddo
              enddo

              n=celdas !numero de nodos ocupados en t=0, es una función del tiempo n=n(t)

            do while ((t.lt.dble(tmax)).and.(n.gt.0))
                t=t+1.0/dble(n) !dble converts integer to double precision
                k=i_ran(n) !escojo nodo ocupado aleatoriamente
                realran=ran_u() !numero aleatorio uniformemente distribuido en [0,1]
                x=lista(k)/l !división entera
                y=mod(lista(k),l)+1 !resto entero de la división +1
                especie=estado(y,x,0)
! VECINO EN EL QUE CREA?---------------------------------
                cruz=i_ran(4)
                a=estado(y,x,cruz) !vecino en el que crea descendencia el nodo k
                x_vec=a/l !columna en la que se encuentra el vecino
                y_vec=mod(a,l)+1 !fila en la que se encuentra el vecino
                control=estado(y_vec,x_vec,5)
                !¿es privilegiado?---------------------
                if (especie.eq.control) then
                  p=p0+epsilon
                else
                  p=p0
                endif
                !------------------------
                if(realran.le.p) then !CREA
                  if (estado(y_vec,x_vec,0).eq.0) then !solo ocupa si (x,y) está vacio
                    n=n+1
                    lista(n)=a !creo en lista() de forma ordenada
                    estado(y_vec,x_vec,0)=especie
                    if (especie.eq.1) then
                      n_a=n_a+1
                    else
                      n_b=n_b+1
                    endif
                  endif

                else !DESTRUYE
                  estado(y,x,0)=0
                  lista(k)=lista(n)!lista()=(x1,x2,...,xn,...,xn,0,...,0)
                  lista(n)=0 !lista()=(x1,x2,...,xn,...,xn-1,0,0,...,0)
                  n=n-1
                  if (especie.eq.1) then
                    n_a=n_a-1
                  else
                    n_b=n_b-1
                  endif
                endif
!DE MOMENTO NO MIDO VARIABLES DINÁMICAS-------------
                !if (int(t).eq.m) then
                  !lt=int(50.0*log10(real(m)))+1
                  !lt=int(50.0*log10(real(t)))+1
                  !mean_n_a(lt)=mean_n_a(lt)+n_a
                  !mean_n_b(lt)=mean_n_b(lt)+n_b
                  !eventos(lt)=eventos(lt)+1
                  !m=m+1 ! volvemos a medir en el proximo paso de MC
                !endif
!---------------------------------------------------------------------
            if(((n_a.eq.0).or.(n_b.eq.0)).and.(aux.eq.1)) then
              aux=0
              mean_t_debil=mean_t_debil+t
              exit !SOLO ME INTERESA MEAN_T_DEBIL!!!!!
            endif

              enddo !fin de una realización
              mean_t_fuerte=mean_t_fuerte+t
              if (aux.eq.1) then
! La simulación ha terminado con coexistencia entre las dos especies
                mean_t_debil=mean_t_debil+t
                coexistencia=coexistencia+1
              endif

            enddo !fin de bucle con índice j sobre realizaciones
            write(*,*) "epsilon=",epsilon
            write(*,*) "t_debil=",mean_t_debil/dble(r)
            write(*,*) "t_fuerte=",mean_t_fuerte/dble(r)
            write(*,*) "COEXISTENCIA=",coexistencia
            write(*,*) "tmax=",tmax

            t=mean_t_debil/dble(r)
            realran=mean_t_fuerte/dble(r)
          write(1,*) epsilon,t,realran,celdas,p0,coexistencia,tmax

            close(1)
          end
