      implicit real*8(a-h,o-z)
      parameter(n=20)
       dimension f(0:n)
             pval(x,y)=atan(x*y)/y-atan(x)
     open(10,file='pm-centers')
     pi=4.*atan(1.0)
      gamma=7.0d0/5.0d0
      alph2=(gamma-1)/(gamma+1)
      alpha=sqrt(alph2)
     pmangle=1.181
     nseries=50000
     k1=int(pmangle*(nseries-1)/1.48+0.5)
     write(*,*)k1
     do k=1,k1+30
     read(10,*)em0,pm0
     enddo 
     x0=sqrt(em0*em0-1.0)
     write(*,*)x0,pm0
     call pmcoefs(gamma,x0,f,n)
     x=f(0)
     do i=1,20
     x=x+f(i)*(pmangle-pm0)**i
     enddo
     write(*,*)pmangle,sqrt(1.0+x*x),pval(x,alpha)
     stop
     end
     
     
              

      subroutine pmcoefs(gam,x0,x,n)
      implicit real*8(a-h,o-z)
      dimension x(0:n)
      dimension u(0:n),v(0:n),t(0:n)
      alph2=(gam-1)/(gam+1)
      alpha=sqrt(alph2)
      opal2=1.0d0+alph2
      omal2=1.0d0-alph2       
      x(0)=x0
      u(0)=x0*x0
      v(0)=u(0)*u(0)
      x(1)=(1.d0+opal2*u(0)+alph2*v(0))/u(0)/omal2
      do k=1,n-1
      u(k)=0.0
      v(k)=0.0
      t(k)=0.0
      do l=0,k
      u(k)=u(k)+x(l)*x(k-l)
      enddo
      do l=0,k
      v(k)=v(k)+u(l)*u(k-l)
      enddo
      do l=0,k-1
      t(k)=t(k)+(l+1)*x(l+1)*u(k-l)
      enddo
     x(k+1)=(opal2*u(k)+alph2*v(k)-t(k)*omal2)/omal2/u(0)/(k+1)
      enddo
      return
      end 
