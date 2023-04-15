      implicit real*8(a-h,o-z)
       parameter (nterms=100,nseries=201)
       dimension em0(0:nseries),ar0(0:nseries)
       dimension f(0:nseries,0:nterms)
              open(10,file='sup-centers')
              open(11,file='sup-series')
              open(12,file='table-sup')
           
      GAMMA=1.4
      df=2.0/(gamma-1.0)
      delta=(df+1.d0)*0.5
      c=(df/(df+1.0d0))**delta
      beta=1.d0/df

              do k1=0,nseries-1
           read(10,*)em0(k1),ar0(k1)
           read(11,*)(f(k1,k),k=0,nterms-1)
           enddo
           tn=0
           error=0.0
           arv=1.0
10           arv=arv+0.01         
!              k1=(arv-1.001)*200.0/(50.0-1.001)+0.5
               k1=(arv-1.001)*(nseries-1)/(25.0-1.001)+0.5
              if(k1.lt.1) k1=1
             tn=tn+1 
      emc=0.0d0
      do k=0,nterms-1
      emc=emc+f(k1,k)*(arv-ar0(k1))**k
      enddo
      error=error+(arv-eps(emc,beta,delta,c))**2
      write(*,*)error,k1
      write(12,21)arv,emc
!,eps(emc,beta,delta,c)
      if(arv.lt.2) go to 10
      arv=2.0-0.1
20    arv=arv+0.1
              k1=(arv-1.001)*(nseries-1)/(25.0-1.001)+0.5
              if(k1.lt.1) k1=1
             tn=tn+1 
      emc=0.0d0
      do k=0,nterms-1
      emc=emc+f(k1,k)*(arv-ar0(k1))**k
      enddo
      error=error+(arv-eps(emc,beta,delta,c))**2
      write(12,21)arv,emc
!,eps(emc,beta,delta,c)
      if(arv.lt.25) go to 20

21   format(1x,f6.2,1x,'&',1x,E20.12,1x,'\','\')
      
      write(*,*)error/tn,error,tn
              stop
              end
              
              function eps(x,b,d,c)
                    implicit real*8(a-h,o-z)
              eps=c*(1.0d0+b*x*x)**d/x
              return
              end
