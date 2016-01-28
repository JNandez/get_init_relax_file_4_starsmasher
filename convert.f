      Program Converter 

      integer ndim 
      parameter (ndim=9)
      double precision M,R,P,rho,T,kap,ATG,TGD
      double precision H,He,C,N,O,Ne,Mg,L,q,a
      double precision eth,enuc,enu,S,U,NA,K
      double precision newS,metals,betab,mu,mue,me,mui
      double precision xx(ndim),Aii(ndim),planck,planckb
      double precision arad,nn(ndim),nz(ndim)
      parameter(NA=6.0221417d23)
      parameter(k=1.380648d-16,
     $  me=9.10938215d-28,pi=3.141592654,
     $  planck=6.6260755d-27,
     $  planckb=planck/(2.d0*pi),
     $  arad=7.5646d-15)

      integer Nmesh,Nvar,Model,mod,i,j,ii
      double precision age,SM,SR,Dover,sum1,roche
      character*50 name1

      write(*,*) 'Write the name of the file *.mdl:'
      read(*,*) name1
     
      write(*,*) 'Write the number of the model:'
      read(*,*) mod 
      
      write(*,*) 'Write the mass ratio M1/M2:'
      read(*,*) q

      write(*,*) 'Write the initial separation:'
      read(*,*) a

      SM=1.98844d33
      SR=6.9598d10
  
      open(2,file=name1,status='unknown')
      open(1,file='eg.last1.muse_s2mm')
      open(3,file='initial.dat')
      open(4,file='entropy.dat')

      read(2,*)Nmesh,Nvar,Dover
      write(*,*)Nmesh,Nvar,Dover

      nn(1)=0
      nn(2)=2
      nn(3)=6
      nn(4)=7
      nn(5)=8
      nn(6)=10
      nn(7)=12
      nn(8)=14
      nn(9)=30

      nz(1)=1
      nz(2)=2
      nz(3)=6
      nz(4)=7
      nz(5)=8
      nz(6)=10
      nz(7)=12
      nz(8)=14
      nz(9)=26


      Aii(1)=1
      Aii(2)=4
      Aii(3)=12
      Aii(4)=14
      Aii(5)=16
      Aii(6)=20
      Aii(7)=24
      Aii(8)=28
      Aii(9)=56

      do j=1,10000

        read(2,*)Model,age 
c        write(*,*)Model,age
        if(Model.eq.mod) then
          write(*,*)'Model ',Model, "age ",age 
          exit
        endif 
        if(Model.gt.mod) then 
          write(*,*)'It is impossible to find that model'
          stop
          exit
        endif
        do i=1,Nmesh,1
         read(2,*)M,R,P,rho,T,kap,ATG,TGD,H,He,C,N,O,Ne,Mg,L,eth,enuc,
     $   enu,S,U
c         write(*,*)M,R,P,rho,T,kap,ATG,TGD,H,He,C,N,O,Ne,Mg,L,eth,enuc,
c     $   enu,S,U
        enddo
      enddo

      do i=1,Nmesh,1

        read(2,*)M,R,P,rho,T,kap,ATG,TGD,H,He,C,N,O,Ne,Mg,L,eth,enuc,
     $   enu,S,U
        metals=0.d0
        mui=0.d0
        mue=0.d0
        sum1=1-(H+He+C+N+O+Ne+Mg)
        xx(1)=H
        xx(2)=He
        xx(3)=C
        xx(4)=N
        xx(5)=O
        xx(6)=Ne
        xx(7)=Mg
        xx(8)=sum1*0.5d0
        xx(9)=sum1*0.5d0
c        mu=28.d0/(40.d0*H+5.d0*He+16)
c        mue=2.d0/(H+1.d0)
        do ii=1,ndim,1
          mue=mue+dble(nz(ii))*xx(ii)/dble(Aii(ii))
          mui=mui+dble(xx(ii))/dble(Aii(ii))
          if (xx(ii).ne.0) then
            metals=metals+(dble(xx(ii))/dble(Aii(ii)))*
     $           log(dble(Aii(ii))**2.5d0/dble(xx(ii)))
          endif 
        end do
        mue=1.d0/mue
        mui=1.d0/mui
        mu=1.d0/(1.d0/mue+1.d0/mui)
        betab=(1.d0/mu)*(2.5d0
     $  +log((k/(NA**(5.d0/3.d0)*2.d0*pi*planckb**2))**1.5d0))
     $  +(1.d0/mue)*log(2.d0*(me*NA)**1.5d0*mue)
     $  +metals
        newS=(1.d0/mu)*log(T**1.5d0/rho)
     $  +(4.d0/3.d0)*(arad/(k*NA))*T**3.d0/rho+betab
c        write(*,*)sum1*0.5d0
        write(1,100)M*SM,R*SR,P,rho,H,He,C,N,O,Ne,Mg,sum1*0.5d0,
     &       sum1*0.5d0
c        write(*,*)H+He+C+N+O+Ne+Mg+sum1*9/19+sum1*10/19 
        write(4,*)R,M,S/(K*NA),U,eth,enuc,enu,T,rho,P,newS,betab,
     $  H,He,C,N,O,Ne,Mg
      enddo
 100  format(e24.16,2x,e21.16,2x,e21.16,2x,e21.16,2x,e21.16,
     &  2x,e21.16,2x,e21.16,2x,e21.16,2x,e21.16,2x,
     &   e21.16,2x,e21.16,2x,e21.16,2x,e21.16)

 150  format(e24.15,2x,e21.15,2x,e21.15,2x,e21.15,2x,e21.15,
     &  2x,e21.15,2x,e21.15,2x,e21.15,2x,
     &   e21.15,2x,e21.15,2x,e21.15,2x,e21.15)
    
      write(3,*)'Mass ratio: ',q
      write(3,*)'Initial separation: ',a
      write(3,*)'Roche radius (Rrlof): ',roche(q)*a
      write(3,*)'Roche radius (r-rlof): ',roche(q)
      write(3,*)'Star radius: ', R
      write(3,*)'Age: ',age,' years'
      write(3,*)'Model: ', mod
      write(3,*)'Total mass:',M

      end

      double precision function roche(q)
      double precision alpha,beta,q
      alpha=2.d0/3.d0
      beta=1.d0/3.d0
      roche=0.49*(q**alpha)/(0.6*(q**alpha)+log(1+(q**beta)))
      return
      end 
