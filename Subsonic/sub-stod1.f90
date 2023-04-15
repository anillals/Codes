      implicit real*8(a-h,o-z)
      parameter(n=50)
       dimension f(0:n)
      eps(x,b,d,c)=c*(1.0d0+b*x*x)**d/x
     open(10,file='sub-centers')
      GAMMA=1.4
      df=2.0/(gamma-1.0)
      delta=(df+1.d0)*0.5
      c=(df/(df+1.0d0))**delta
      beta=1.d0/df
     aratio=10.0
     nseries=200
     k1=(aratio-1.001)*(nseries-1)/(50.0-1.001)+0.5
     write(*,*)k1
     do k=1,k1
     read(10,*)em0,ar0
     enddo 
     write(*,*)em0,ar0
     call stodcoefs(em0,ar0,beta,f,n)
     em=f(0)
     do i=1,50
     em=em+f(i)*(aratio-ar0)**i
     enddo
     write(*,*)aratio,em,eps(em,beta,delta,c)
     stop
     end
     
     
              

      subroutine stodcoefs(em0,eps0,beta,f,n)
      implicit real*8(a-h,o-z)
      dimension f(0:n),p(0:n),q(0:n),u(0:n),v(0:n)   
      f(0)=em0
      den=eps0*(em0*em0-1.0)
      do k=0,n-2
      p(k) =0.0d0
      q(k) =0.0d0
      u(k)=0.0d0 
      v(k)=0.0    
      
      do l=0,k
      p(k)=p(k)+f(l)*f(k-l)
      enddo
 
      do l=0,k
      q(k)=q(k)+p(l)*f(k-l)
      enddo   
        
      do l=0,k
      u(k)=u(k)+l*f(l)*p(k-l)
      enddo   
 
      do l=0,k-1
      if(k.gt.0) v(k)=v(k)+(l+1)*f(l+1)*p(k-l)
      enddo
     f(k+1)=((k+1)*f(k)+beta*q(k)-u(k)-eps0*v(k))/den/(k+1.d0)
      enddo
      return
      end
