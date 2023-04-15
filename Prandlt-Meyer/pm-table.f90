      implicit real*16(a-h,o-z)
       parameter (nterms=10,nseries=50000)
       dimension em0(0:nseries),pm0(0:nseries)
       dimension f(0:nseries,0:nterms)
             pval(x,y)=atan(x*y)/y-atan(x)
     open(10,file='pm-centers')
      open(11,file='pm-series')
      
              open(12,file='table-PM')
           pi=4.*atan(1.0)
      gamma=7.0d0/5.0d0
      alph2=(gamma-1)/(gamma+1)
      alpha=sqrt(alph2)
              do k1=0,nseries-1
           read(10,*)em0(k1),pm0(k1)
           read(11,*)(f(k1,k),k=0,nterms-1)
           enddo
           tn=0
           error=0.0
           pmv=0.0
10           pmv=pmv+0.01         
               k1=int(pmv*(nseries-1)/1.48+0.5)
              if(k1.lt.1) k1=1
              
        
             tn=tn+1 
      emc=f(k1,0)
      do k=1,nterms-1
     emc=emc+f(k1,k)*(pmv-pm0(k1))**k
      enddo
      emn=sqrt(1+emc*emc)
               
                 
                       write(*,*)pmv,emn,pval(emc,alpha)
      error=error+(pmv-pval(emc,alpha))**2
  
      write(12,21)pmv,emn
      if(pmv.lt. 1.47) go to 10

21   format(1x,f8.4,1x,'&',1x,E20.12,1x,'\','\')
      
      write(*,*)error/tn,error,tn
              stop
              end
              

