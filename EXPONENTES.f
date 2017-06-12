        program diagrama_fase
          implicit none
          integer*4 k,cruz,i,a,j,m !k y cruz guardan los numeros aleatorios/ a,j,m=indices
!/////////////////////////////////////////////////////
          integer*4 lt,n,tmax,l,r
          parameter (l=10000)
          !n=sitios activos
          !tmax=tiempo maximo de ejecucion
          !lt="tiempo logaritmico"
          !l= tamaño de la red 1-D
          !r=numero de realizaciones
!/////////////////////////////////////////////////////////////////////7
          real*8 p,t
          !p=probabilidad de que una modificación sea una creacion
          !t=tiempo
          real*8 realran,sumaR
          !Auxiliares
          real ran_u
          integer*4 i_ran
          integer*4 lista(l),state(l)!CAMBIAR
          real*8 mean_n(500),mean_p(500),mean_R(500)
          integer*4 eventos(500),eventos_total(500)
          integer*4 vecino(l,2) !l filas y tantas columnas como vecinos

          !Seed de numeros aleatorios, abrimos fichero
          !call ran_ini(time())
          call ran_ini(104723)
          open(unit=1, file="n(t).dat", status="unknown")
          open(unit=2, file="eventos.dat", status="unknown")
          !CARACTERISTICAS DEL SISTEMA
          r=10000
          tmax=17000
          p=0.76732 !p=pc
          mean_n=0.0
          mean_p=0.0
          mean_R=0.0
          eventos=0
          eventos_total=0

          !/-------------------_/
          !Relleno matriz vecino
          !donde vecino(i,j)=coordenada de uno de los vecinos del nodo i-ésimo
          !/--------------------/
          do i=1,l
            vecino(i,1)=i-1
            vecino(i,2)=i+1
          enddo

!Para evitar correlaciones, solo mido en cada paso de Montecarlo
            do i=1,int(tmax)
              t=dble(i)
              lt=int(50.0*log10(t))+1
          eventos_total(lt)=eventos_total(lt)+1
            enddo

            do lt=1,212
        write(2,*) lt,10.0**(((dble(lt))-1)/50.0),eventos_total(lt)
            enddo

            do j=1,r
              !CONDICIONES INICIALES spreading
              m=2
              sumaR=0.0
              t=1.0
              state=0
              lista=0
              n=1 !numero de nodos ocupados, es una función del tiempo n=n(t)
              lista(1)=l/2 ! vector ordenado (x1,x2,x3,---,xn,o,o,---,o)
              state(l/2)=1 ! state(i)=i si i es un nodo ocupado
                           ! state(i)=0 si i es un nodo vacío



            do while ((t.lt.dble(tmax)).and.(n.gt.0))
                t=t+1.0/dble(n) !dble converts integer to double precision
                k=i_ran(n) !escojo nodo ocupado aleatoriamente
                realran=ran_u() !numero aleatorio uniformemente distribuido en [0,1]
                if(realran.le.p) then !CREA
                  cruz=i_ran(2)
                  a=vecino(lista(k),cruz) !vecino en el que crea descendencia el nodo k
                  if (state(a).eq.0) then !si state(a) ya esta ocupado no se ocupa (unica utilidad de state)
                    n=n+1
                    lista(n)=a !creo en lista() de forma ordenada
                    state(a)=1
                    sumaR=sumaR+(a-l/2)**2 !si se crea, se actualiza sumaR
                    !CONDICIONES DE CONTORNO, interrumpo bucle si toco frontera
                    if ((a.eq.l).or.(a.eq.1)) then !CONDICIONES DE CONTORNO
                      exit !si toca frontera interrumpo loop
                    endif
                  endif

                else !DESTRUYE
                  sumaR=sumaR-(lista(k)-l/2)**2 !si se destruye, se actualiza sumaR
                  state (lista(k))=0 !si uso lista(k) dentro de state, !he de eliminar primero en state!
        !he de borrar en lista() de forma ordenada lista()=(x1,x2,...,xk,...,xn,o,...,0) he de eliminar xk
                  lista(k)=lista(n)!lista()=(x1,x2,...,xn,...,xn,0,...,0)
                  lista(n)=0 !lista()=(x1,x2,...,xn,...,xn-1,0,0,...,0)
                  n=n-1
                endif
                if (int(t).eq.m) then
                lt=int(50.0*log10(real(m)))+1
                mean_n(lt)=mean_n(lt)+n
                mean_R(lt)=mean_R(lt)+sumaR
                eventos(lt)=eventos(lt)+1
                m=m+1 ! volvemos a medir en el proximo paso de MC
                endif


              enddo !fin de una realización


        enddo !fin de bucle con índice i sobre realizaciones

          do j=1,212 !CAMBIAR!
            if (eventos(j).ne.0) then
            mean_p(j)=dble(eventos(j))/dble(eventos_total(j))
            realran=mean_p(j)*mean_n(j)/(dble(eventos(j))*dble(r))
            sumaR=mean_p(j)*mean_n(j)/(dble(r))
            t=10.0**(((dble(j))-1)/50.0)


            !write(1,*) t,realran,mean_p(j)/dble(r),mean_R(j)/mean_n(j)
        write(1,*) t,realran,mean_p(j)/dble(r),mean_R(j)/mean_n(j)

            write(2,*) j,t,eventos(j)
            endif
          enddo
          close(1)
          close(2)
        end
