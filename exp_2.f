        program contact_process
          integer k,cruz,n,nbin,tmax,l,i,a,r,j,sumR,m
          real p,t,realran,h,pbin,Rbin
          integer lista(10000),sigma(10000)
          real mean_n(0:30),mean_p(0:30),mean_R(0:30)

          !Seed de numeros aleatorios, abrimos fichero
          call ran_ini(time())
          open(unit=1, file="data_p.dat", status="unknown")
          open(unit=2, file="data_n.dat", status="unknown")
          open(unit=3, file="data_R.dat", status="unknown")
          !CARACTERISTICAS DEL SISTEMA
          p=0.76732 !p=pc
          l=10000 !Ha de ser par
          r=10000
          tmax=10000
          mean_n=0.0 !media de n sobre las realizaciones para un tiempo t
          mean_p=0.0
          mean_R=0.0
          do j=1,r
            !CONDICIONES INICIALES
            sumR=0 !sumR= sumatorio (distancia a L/2)^2
            i=0 !indice que lleva el temaño de cada bin
            h=1 !primer bin
            t=0.0
            lista=0
            sigma=0
            n=1 !numero de nodos ocupados, es una función del tiempo n=n(t)
            nbin=0 !media temporal de n en un bin
            Rbin=0
            pbin=0
            lista(1)=l/2
            sigma(l/2)=1
            do while ((t.lt.real(tmax)).and.(n.ge.1))
              t=t+(1.0)/real(n)
              k=i_ran(n) !escojo nodo ocupado aleatoriamente
              realran=ran_u() !numero aleatorio uniformemente distribuido en [0,1]
              if(realran.le.p) then !CREA
                cruz=i_ran(2) !numero aleatorio entero entre 1 y 2, no puedo-
                if (cruz.eq.1) then !-guardarlo en realran porque realran es real, podria haber problemas
                  a=lista(k)+1!crea a la derch
                else
                  a=lista(k)-1!crea a la izq
                endif
                if (sigma(a).eq.0) then !si sigma(a) ya esta ocupado no se ocupa (unica utilidad de sigma)
                  n=n+1
                  lista(n)=a
                  sigma(a)=1
                  sumR=sumR+(a-l/2)**2 !si se crea, se actualiza sumR
                  !CONDICIONES DE CONTORNO: el programa finaliza al tocar la forntera
                endif
                if ((a.eq.l).or.(a.eq.1)) then
                  goto 111
                endif

              else !DESTRUYE
                sumR=sumR-(lista(k)-l/2)**2 !si se destruye, se actualiza sumR
                sigma (lista(k))=0 !si uso lista(k) dentro de sigma, !he de eliminar primero en sigma!
                lista(k)=lista(n)
                lista(n)=0 !he de borrar en lista() de forma ordenada
                n=n-1
              endif
              if (n.ge.1) then !veo si el sistema esta o no en la fase absorvente, actualizo media sobre los bin
                pbin=pbin+1/real(n)
                nbin=nbin+1 !nbin se incrementa siempre en 1 (n(t)*Dt=1), si n(t) diferente de 0
                Rbin=Rbin+real(sumR)/real(n)
              endif
              if (t.ge.h) then !GUARDO INFORMACIÓN DEL BIN
                mean_p(i)=mean_p(i)+pbin/h
                mean_n(i)=mean_n(i)+real(nbin)/h
                mean_R(i)=mean_R(i)+Rbin/h
                i=i+1
                h=(1.5**i) !h=hmin*base^a
                nbin=0
              endif
            enddo

111         if ((t.lt.real(tmax)).and.((t.ne.h))) then !si t<tmax se ha activado la CC, si t/=h no se actualizaron  mean_
                mean_p(i)=mean_p(i)+pbin/h
                mean_n(i)=mean_n(i)+real(nbin)/h
                mean_R(i)=mean_R(i)+Rbin/h
            endif

          enddo

          do i=0,30
            h=(1.5**i)
            write(1,*) h, mean_p(i)/real(r) !CAMBIAR INDICES
            write(2,*) h, mean_n(i)/real(r)
            write(3,*) h, mean_R(i)/real(r)
          enddo
          close(1)
          close(2)
          close(3)
        end
