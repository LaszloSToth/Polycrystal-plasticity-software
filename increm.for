c     
c
c  deformation texture simulation based on rate sensitive slip
c                               Laszlo S. Toth, 2002.febr. Metz 
c  simplified program for strain rate sensitive crystal plasticity
c   - no twinning  
c   - no random orientation generation
c   - can be used as a subroutine for a given polycrystal for one increment
c   - input: 
c          - orientations with stress state and crss values
c          - available slip systems
c          - imposed velocity gradient
c          - time increment
c          - relaxation matrix in case relaxed constrains conditions are applied
c          - and so on, see below


      
      subroutine increm
     &(ngrain,nsys,mc,nab,h,gam0,nhard,irel,vel,fi1,fi,fi2,f,tauc0,igu,dt)
     	include 'POLYCR.DIM'
     
c input parameters:
c ngrain: total number of grains
c nsys= total number of slip systems
c nab(48,48) hardening matrix
c h= strain rate sensitivity index: h=1/m
c nhard: type of hardening theory used
c      nhard=1 Broknkhorst type with latent hardening
c      nhard=2 dislocation cell model 
c irel(i,j): the relaxed constraints condition in a matrix form: 
c     if a component is relaxed, it is =1
c     if a component is not relaxed, it is =0
c vel(i,j): prescribed velocity gradient    
c fi1(ngrain),fi(ngrain),fi2(ngrain),Euler angles of each grain
c f(ngrain): volume fraction of grains  
c st(ngrain,3,3): stress state in each grain (during the previous deformation step)
c roc(ngrain),row(ngrain): dislocation densities in cell interior and cell wall within each grain if cell hardening is used 
c     to use in common block with cell hardening subroutine "estrin"
c fw(ngrain): volume fraction of cell walls in each grain if cell hardening is used 
c     to use in common block with cell hardening subroutine "estrin"                                    
c ressum(ngrain): accummulated resolved shear strain in each grain, for the use in cell hardening
c     to use in common block with cell hardening subroutine "estrin"
c tauc0(ngrain,48): the reference stresses in each slip system (maximum 48 systems)
c igu: 0 or 1, if it is zero, it means that there is no initial stress state (no initial guess)
c              if it is 1, there is a previous stress state to use in the present increment

c output parameters:
c fi1,fi,fi2,row,roc,fw,ressum,tauc: with renewed values
c avst(i,j): average stress state on the polycrystal
                
      implicit real*8 (a-h,o-z)

      integer*2 nab   
c nab(48,48) hardening matrix      

      real*8 mc,ms,m,k0

c grain orientations (three Euler angles and volume of each grain):
      dimension fi1(NGRMX),fi(NGRMX),fi2(NGRMX),f(NGRMX),sc(3,3)

c reference strengths of slip systems, tauc0.  
      dimension tauc0(NGRMX,NSYSMX),tauc(NSYSMX)
                  
      dimension vel(3,3),irel(3,3),isrel(5),isnc(5)
      dimension mc(3,3,NSYSMX),ms(3,3,NSYSMX),tab(3,3),sns(NSYSMX,3)
      dimension st1(3,3)
      dimension nab(NSYSMX,NSYSMX),gamr(NSYSMX)
      dimension e(5),s(5),m(NSYSMX,5),gam(NSYSMX),stmem(NGRMX,5) 
      common/cell1/ res,ressum(NGRMX),row(NGRMX),roc(NGRMX),f1,f2,
     &res_bar,alfa,G,b,hm,hn,holt,alfs,bets,k0,res0,kvol
      common/cell2/tres,f_inf,res1,res2,fw(NGRMX)
c /cell1/,/cell2/: hardening parameters for dislocation-cell type hardening
      common/bronk/ q(4),h0,tsat,a
c /bronk/: hardening parameters for Bronkhorst type of hardening
c /slip1/: common block for slip distribution
      common/slip1/ slip(NGRMX,NSYSMX) 
	common/velcr/ velg(NGRMX,3,3),fgrad(NGRMX,3,3)			   ! velg: grain velocity gradient made by slip, fgrad: finite grain strain gradient
	common/strs/ st(NGRMX,3,3),avst(3,3),defg(NGRMX,3,3),avdef(3,3),avvelg(3,3)   
	common/trace/ sn(NSYSMX,3),beta(NGRMX,NSYSMX)  ! slip plane normals in crystal reference system, not normalized, and the orientation of the slip line with respect to the x-axis (the surface plane must be x-y)

c constants:
      pi=dacos(-1.d0)
      deg=180.d0/pi  
      sq2=sqrt(2.d0)
      sq3=sqrt(3.d0)
    
c fixed parameters:
c tol: required precision in the stress state     
      tol=1.d-6        

c nil initial values:
      avst=0.d0		! average stress state
      avdef=0.d0		! average obtained strain rate tensor
	avvelg=0.d0 	! average obtained velocity gradient tensor 

c checks if the relaxed conditions are allright:
      call relcheck(irel)

c total volume of grains:
      fsum=0.d0
      do i=1,ngrain
      fsum=fsum+f(i)
      end do                                         
c the prescribed strain rate in vector form:
      call evector(vel,irel,e,nc,is,isrel,isnc)
        
c inputs:                                          
c vel(3,3): velocity gradient
c irel(3,3): relaxation matrix
c outputs:
c isrel(5)= relaxation vector for the stress state
c isnc(5)= relaxation vector to use in the construction of the stress vector 
c is= also used in the construction of the stress vector 
c e(5): the prescribed strain rate vector
c nc=number of prescribed strain rate components

c*********************************************************************
c big loop for the number of grain orientations

      do 5000 ig=1,ngrain

c      write(*,'(a,i5)') '  grain  ',ig

c conversion of Euler angles from degrees to radians and construction of the 
c transformation matrix sc(i,j): transformation from the sample to the crystal system:
      afi1=fi1(ig)/deg
      afi=fi(ig)/deg
      afi2=fi2(ig)/deg
      call euler(sc,afi1,afi,afi2)

c normalise tauc's so that the maximum reference strength is 1.d0. (Original values remain in tauc0)
      taumax=0.d0
      do i=1,nsys
      if(tauc0(ig,i).gt.taumax) taumax=tauc0(ig,i)
      end do
      do i=1,nsys
      tauc(i)=tauc0(ig,i)/taumax
      end do

c expression of the slip systems normal vectors in the sample reference system (for slip line identification)
      do 160 L=1,nsys
      do 160 i=1,3
	sns(L,i)=0.d0
      do 160 j=1,3
160   sns(L,i)=sns(L,i)+sc(j,i)*sn(L,j)
c orientation of the slip line in the x-y reference system, with respect to x, positive angle
      do 161 L=1,nsys
	if((sns(L,1).eq.0).and.(sns(L,2).eq.0)) then
	beta(ig,L)=720.	 ! particular case, when the slip plane is perpendicluar to the sample surface
      else
	beta(ig,L)=datan2D(sns(L,2),(-sns(L,1)))
	end if
	if(beta(ig,L).lt.0.) beta(ig,L)=beta(ig,L)+180.
161   continue

c transformation of the slip systems into the sample reference system
      do 17 L=1,nsys
      do 16 i=1,3
      do 16 j=1,3
      tab(i,j)=0.d0
      do 16 k=1,3
16    tab(i,j)=tab(i,j)+mc(i,k,L)*sc(k,j)
      do 17 i=1,3
      do 17 j=1,3
      ms(i,j,L)=0.d0
      do 17 k=1,3
17    ms(i,j,L)=ms(i,j,L)+sc(k,i)*tab(k,j)

c assign the slip system M vectors
      do i=1,nsys
      call mvectors(ms,irel,i,m,sq2,sq3)
      end do
c input of mvectors:
c ms(3,3,slip system index) : Schmid orientation tensor in the sample system
c i: slip system index                            
c irel(3,3): matrix for the relaxed constraints condition
c output:
c m(slips system index,5) : Schmid orientation tensor in five component vector form      
c                           in the sample system

c use previous stress state as guess in the present step: 
      do i=1,5
      s(i)=stmem(ig,i)
      end do
           
c call the rate sensitive subroutine to solve the single crystal problem:

      call rtsens(gam0,h,e,s,nc,m,tauc,nsys,igu,gam,tol,info)  
                  

c*******************************
c results obtained from subroutine rtsens:
c gam(i)= strain rate in slip system i
c s(5)= stress vector for the grain in the sample reference system
c slip distribution:
      do i=1,nsys
      slip(ig,i)=gam(i)
	end do

      gs=0.d0                    
c gs= sum of increments of absolute value shears      
      gsum=0.d0
c gsum= sum of slip rate absolute values      
      res=0.d0           
c res= sum of absolute values of resolved shear strain rates 
      resmax=0.d0      ! to find maximum resolved shear strain rate

      do 664 i=1,nsys 
      res=res+dabs(gam(i))
      gamr(i)=gam(i)   
c gamr(i): shear rate in system i      
      if(dabs(gamr(i)).gt.resmax) resmax=dabs(gamr(i))
      gam(i)=gam(i)*dt
c gam(i): shear increment in system i      
      gs=gs+dabs(gam(i))
664   gsum=gsum+dabs(gamr(i))

c acummulated total resolved shear strain in grain igr 
      ressum(ig)=ressum(ig)+gs

c recalculates the imposed strain rate (for checking):
      do 60 i=1,3
      do 60 j=1,3
      defg(ig,i,j)=0.d0
      do 60 L=1,nsys
60    defg(ig,i,j)=defg(ig,i,j)+(ms(i,j,L)+ms(j,i,L))*gamr(L)/2.d0
c def(ig,i,j): obtained strain rate

c average strain rate:
      do i=1,3
      do j=1,3
      avdef(i,j)=avdef(i,j)+defg(ig,i,j)*f(ig)/fsum
      end do
      end do   
      
c Stresses:
c Conversion from the nc component stress vector into 3x3 tensor:

      call sconv(isrel,is,isnc,s,nc,st1,sq2,sq3)
      
c memorize s(i) to use in the next increment:
      do i=1,5
      stmem(ig,i)=s(i)      
      end do
      
c cell-hardening is employed:
      if(nhard.eq.2) then
        call hardestrin(ig,dt)
        do i=1,nsys
        tauc0(ig,i)=tres
        end do
        taumax=tres
      end if

c rescale the stress state according to taumax:
       do 987 i=1,3
       do 987 j=1,3
987    st(ig,i,j)=st1(i,j)*taumax

c average stress:
      do i=1,3
      do j=1,3
      avst(i,j)=avst(i,j)+st(ig,i,j)*f(ig)/fsum
      end do
      end do

c orientation change:
      if(dt.eq.0.d0) go to 529
      call rot(vel,ig,dt,irel,ms,gamr,nsys,afi1,afi,afi2)
529   continue
c conversion of new Euler angles from radians to degrees 
      fi1(ig)=afi1*deg
      fi(ig)=afi*deg
      fi2(ig)=afi2*deg
             

c average obtained velocity gradient:
      do i=1,3
      do j=1,3
      avvelg(i,j)=avvelg(i,j)+velg(ig,i,j)*f(ig)/fsum
      end do
      end do

c update crss in case of Bronkhorst type of hardening
c ig: index of grain
      if(nhard.eq.1) then
       do 8633 i=1,nsys
       taudot=0.d0
       do 8632 j=1,nsys
       k=nab(i,j)
8632   taudot=taudot+dabs(gamr(j))*q(k)*h0*(dabs(1.d0-tauc0(ig,i)/tsat))**a
       tauc0(ig,i)=tauc0(ig,i)+taudot*dt       
       if(tauc0(ig,i).gt.tsat) tauc0(ig,i)=tsat
8633  continue
      end if

5000  continue

      return
      end


c***********************************************
c subroutine to prescribe the strain rate vector:

      subroutine evector(vel,irel,e,nc,is,isrel,isnc)
      implicit real*8 (a-h,o-z)
      dimension vel(3,3),irel(3,3),e(5),isrel(5),isnc(5)

c inputs:                                          
c vel(3,3): velocity gradient
c irel(3,3): relaxation matrix

c outputs:
c isrel(5)= relaxation vector for the stress state
c isnc(5)= relaxation vector to use in the construction of the stress vector 
c is= also used in the construction of the stress vector 
c e(5): the prescribed strain rate vector  

      do i=1,5
      isrel(i)=1    ! 1 means: relaxed, 0 means: prescribed
      end do
c for isrel(1)=-1 and isrel(2)=-1 means that two diagonal components of the
c strain rate vector are relaxed: special case, handled by the permutation
c factor: 'is' (=1,2, or 3)
      is=0
c isnc(i), i=1,nc: contains the order of the non-zero stress components

      sq2=dsqrt(2.d0)
      sq3=dsqrt(3.d0)
      nc=0      ! number of prescribed components of the strain rate vector
c for the diagonal components:
      if(irel(1,1).eq.0.and.irel(2,2).eq.0.and.irel(3,3).eq.0) then
        nc=nc+1
        isrel(1)=0   ! 0 means: this deformation component is prescribed
        isnc(nc)=1
        e(nc)=(vel(2,2)-vel(1,1))/sq2
        nc=nc+1
        isnc(nc)=2
        isrel(2)=0
        e(nc)=vel(3,3)*sq3/sq2
      end if
      if(irel(1,1).eq.1.and.irel(2,2).eq.1.and.irel(3,3).eq.1) go to 300

c case when two diagonal components of the strain rate tensor are relaxed:
      if(irel(1,1).eq.1.and.irel(2,2).eq.1.and.irel(3,3).eq.0) then
        isrel(1)=-1
        isrel(2)=-1
        isnc(1)=-1
        nc=nc+1
        is=3
        e(nc)=vel(3,3)*sq3/sq2
      end if
      if(irel(1,1).eq.1.and.irel(2,2).eq.0.and.irel(3,3).eq.1) then
        isrel(1)=-1
        isrel(2)=-1
        isnc(1)=-1
        nc=nc+1
        is=2
        e(nc)=vel(2,2)*sq3/sq2
      end if
      if(irel(1,1).eq.0.and.irel(2,2).eq.1.and.irel(3,3).eq.1) then
        isrel(1)=-1
        isrel(2)=-1
        isnc(1)=-1
        nc=nc+1
        is=1
        e(nc)=vel(1,1)*sq3/sq2
      end if

c shear components:
300   if(irel(2,3).eq.0.and.irel(3,2).eq.0) then
        nc=nc+1
        isrel(3)=0
        isnc(nc)=3
        e(nc)=(vel(2,3)+vel(3,2))*sq2/2.d0
      end if
      if(irel(1,3).eq.0.and.irel(3,1).eq.0) then
        nc=nc+1
        isrel(4)=0
        isnc(nc)=4
        e(nc)=(vel(1,3)+vel(3,1))*sq2/2.d0
      end if
      if(irel(1,2).eq.0.and.irel(2,1).eq.0) then
        nc=nc+1
        isrel(5)=0
        isnc(nc)=5
        e(nc)=(vel(1,2)+vel(2,1))*sq2/2.d0
      end if

      return
      end


c**************************************************
c subroutine to prescribe the slip system orientation vectors (m-vectors)

      subroutine mvectors(ms,irel,i,m,sq2,sq3)
	include 'POLYCR.DIM'
      implicit real*8 (a-h,o-z)
      real*8 ms,m
      dimension ms(3,3,NSYSMX),irel(3,3),m(NSYSMX,5)
      
c input:
c ms(3,3,slip system index) : Schmid orientation tensor in the sample system
c i: slip system index                            
c irel(3,3): matrix for the relaxed constraints condition

c output:
c m(slips system index,5) : Schmid orientation tensor in five component vector form      
c                           in the sample system


      nc=0      ! number of prescribed components of the strain rate vector
      if(irel(1,1).eq.0.and.irel(2,2).eq.0.and.irel(3,3).eq.0) then
        nc=nc+1
        m(i,nc)=(ms(2,2,i)-ms(1,1,i))/sq2
        nc=nc+1
        m(i,nc)=ms(3,3,i)*sq3/sq2
      end if
      if(irel(1,1).eq.1.and.irel(2,2).eq.1.and.irel(3,3).eq.1) go to 300
      if(irel(1,1).eq.1.and.irel(2,2).eq.1.and.irel(3,3).eq.0) then
        nc=nc+1
        m(i,nc)=ms(3,3,i)*sq3/sq2
      end if
      if(irel(1,1).eq.1.and.irel(2,2).eq.0.and.irel(3,3).eq.1) then
        nc=nc+1
        m(i,nc)=ms(2,2,i)*sq3/sq2
      end if
      if(irel(1,1).eq.0.and.irel(2,2).eq.1.and.irel(3,3).eq.1) then
        nc=nc+1
        m(i,nc)=ms(1,1,i)*sq3/sq2
      end if
300   if(irel(2,3).eq.0.and.irel(3,2).eq.0) then
        nc=nc+1
        m(i,nc)=(ms(2,3,i)+ms(3,2,i))/sq2
      end if
      if(irel(1,3).eq.0.and.irel(3,1).eq.0) then
        nc=nc+1
        m(i,nc)=(ms(1,3,i)+ms(3,1,i))/sq2
      end if
      if(irel(1,2).eq.0.and.irel(2,1).eq.0) then
        nc=nc+1
        m(i,nc)=(ms(1,2,i)+ms(2,1,i))/sq2
      end if

      return
      end


      SUBROUTINE EULER(SC,FI1,FI,FI2)
      IMPLICIT REAL*8 (A-H,O-Z)
      DIMENSION SC(3,3)
      CFI= DCOS(FI)
      SFI= DSIN(FI)
      CFI1=DCOS(FI1)
      SFI1=DSIN(FI1)
      CFI2=DCOS(FI2)
      SFI2=DSIN(FI2)
      SC(1,1)=CFI1*CFI2-SFI1*SFI2*CFI
      SC(1,2)=SFI1*CFI2+CFI1*SFI2*CFI
      SC(1,3)=SFI2*SFI
      SC(2,1)=-CFI1*SFI2-SFI1*CFI2*CFI
      SC(2,2)=-SFI1*SFI2+CFI1*CFI2*CFI
      SC(2,3)=CFI2*SFI
      SC(3,1)=SFI1*SFI
      SC(3,2)=-CFI1*SFI
      SC(3,3)=CFI
      RETURN
      END

c**********************************************************************
c subroutine to check if the relaxed constraint conditions are allright:

      subroutine relcheck(irel)
      dimension irel(3,3)

c input:
c irel(i,j): relaxation matrix for the relaxed constraints conditions

c output:
c there is no output; if the relaxation matrix is incorrect, the calculation stops
c with a message

      n=0
      do 1 i=1,3
      do 1 j=1,3
1     if(irel(i,j).eq.1) n=n+1
      if(n.eq.9) then
         write(*,100)
         write(3,100)
         stop
      end if
100   format(/'   You cannot relax all velocity components!')
      n=0
      do 2 i=1,3
2     if(irel(i,i).eq.1) n=n+1
      if(n.eq.1) then
         write(*,200)
         write(3,200)
         stop
      end if
200   format(/'   Only one normal component in the velocities cannot be
     !relaxed; must be at least two!')
      return
      end

c**********************************************************************
c subroutine to convert the nc component stress vector into 3x3 tensor:

      subroutine sconv(isrel,is,isnc,s,nc,st,sq2,sq3)
      implicit real*8 (a-h,o-z)
      dimension isrel(5),isnc(5),s(5),snew(5),st(3,3)

      do 5 i=1,5
5     snew(i)=0.d0

      if(nc.eq.5) then
          do 20 i=1,5
20        snew(i)=s(i)
          go to 500
      end if

c two normal strains relaxed:
      if(isrel(1).eq.-1) go to 300

      do 10 i=1,nc
10    snew(isnc(i))=s(i)

500   st(1,1)=(-snew(1)-snew(2)/sq3)/sq2
      st(2,2)=(snew(1)-snew(2)/sq3)/sq2
      st(3,3)=snew(2)*sq2/sq3
      st(2,3)=snew(3)/sq2
      st(3,1)=snew(4)/sq2
      st(1,2)=snew(5)/sq2
      st(3,2)=snew(3)/sq2
      st(1,3)=snew(4)/sq2
      st(2,1)=snew(5)/sq2
      return

c two normal strains relaxed:
300   do 301 i=2,nc
301   snew(isnc(i))=s(i)

      do 302 i=1,3
302   st(i,i)=0.d0
      st(is,is)=s(1)*sq2/sq3

      st(2,3)=snew(3)/sq2
      st(3,1)=snew(4)/sq2
      st(1,2)=snew(5)/sq2
      st(3,2)=snew(3)/sq2
      st(1,3)=snew(4)/sq2
      st(2,1)=snew(5)/sq2
      return
      end

c ***********************************************
c subroutine to calculate the orientation change:

      subroutine rot(vel,ig,dt,irel,ms,gamr,nsys,fi1,fi,fi2)
	include 'POLYCR.DIM'
      implicit real*8 (a-h,o-z)
      real*8 ms
      dimension vel(3,3),irel(3,3),ms(3,3,NSYSMX),gamr(NSYSMX),om(3,3)
      dimension ilt(3,3),i1(3),i2(3),temp1(3,3),temp2(3,3)
  	  common/velcr/ velg(NGRMX,3,3),fgrad(NGRMX,3,3)			   ! velg: grain velocity gradient made by slip, fgrad: finite grain strain gradient

c irel: relaxation matrix
c vel: velocity gradient
c dt: time increment
c ms: slip system matrices in the sample system
c gam: slips
c fi1,fi,fi2: Initial (and new) Euler angles
c om: lattice rotation matrix
c om1,om2,om3: the 3 components of the lattice spin vector
c p1,p2,p3: the 3 components of the plastic spin vector
c r1,r2,r3: the 3 components of the rigid body rotation vector
c velg(ig,3,3)=velocity gradient of the grain	made by slip
c fgrad: finite grain strain gradient


      pi2=2.d0*dacos(-1.d0)
      do 1 i=1,3
      do 1 j=1,3
      ilt(i,j)=0
1     om(i,j)=0.d0

c rigid body rotation:
      r1=(vel(2,3)-vel(3,2))/2.d0
      r2=(vel(3,1)-vel(1,3))/2.d0
      r3=(vel(1,2)-vel(2,1))/2.d0

      n=0
      do 10 i=1,3
      do 10 j=1,3
      if(i.ne.j.and.irel(i,j).eq.0.and.ilt(i,j).ne.1.and.n.le.3) then
          ilt(i,j)=1    ! this component is used up
          ilt(j,i)=1
          n=n+1
          i1(n)=i
          i2(n)=j
          om(i,j)=0.d0
          do 20 L=1,nsys
20        om(i,j)=om(i,j)+ms(i,j,L)*gamr(L)
          om(i,j)=vel(i,j)-om(i,j)
      end if
10    continue

C om is antisymmetric:
      do 30 i=1,n
30    om(i2(i),i1(i))=-om(i1(i),i2(i))

c	do i=1,3
c      Write(3,'(3f12.3)') (om(i,j),j=1,3)
c	end do

		do i=1,3
		do j=1,3
		temp1(i,j)=0.d0
          do L=1,nsys
		temp1(i,j)=temp1(i,j)+ms(i,j,L)*gamr(L)
		end do
        velg(ig,i,j)=temp1(i,j)+om(i,j)
		end do
		end do


		do i=1,3
		do j=1,3
		temp2(i,j)=0.d0
		do k=1,3
		temp2(i,j)=temp2(i,j)+velg(ig,i,k)*fgrad(ig,k,j)
		end do
		end do
		end do


		do i=1,3
		do j=1,3
		fgrad(ig,i,j)=fgrad(ig,i,j)+dt*temp2(i,j)
		end do
		end do	

      om1=om(2,3)
      om2=om(3,1)
      om3=om(1,2)

c plastic spin-rotation:
      p1=r1-om1
      p2=r2-om2
      p3=r3-om3
c      if(inf.eq.1) write(3,791) p1,p2,p3
c791   format(' Components of plastic spin: ',3e14.5/)

      if(dabs(dsin(fi)).lt.1.d-3) then
          fi2dot=0.d0
          go to 50
      end if
      fi2dot=-om(2,3)*dsin(fi1)/dsin(fi)+om(3,1)*dcos(fi1)/dsin(fi)
50    fidot=-om(2,3)*dcos(fi1)-om(3,1)*dsin(fi1)
      fi1dot=-om(1,2)-fi2dot*dcos(fi)

      fi1=dt*fi1dot+fi1
      fi=dt*fidot+fi
      fi2=dt*fi2dot+fi2
      if(dabs(fi1).gt.pi2) fi1=fi1-dsign(1.d0,fi1)*pi2
      if(dabs(fi).gt.pi2) fi=fi-dsign(1.d0,fi)*pi2
      if(dabs(fi2).gt.pi2) fi2=fi2-dsign(1.d0,fi2)*pi2
      return
      end

c
c    multiphase cell model for stage II-III-IV-V hardening

      subroutine hardestrin(ig,dt)
	include 'POLYCR.DIM'
      implicit real*8(a-h,o-z)
      real*8 k0
      common/cell1/ res,ressum(NGRMX),row(NGRMX),roc(NGRMX),f1,f2,
     &res_bar,alfa,G,b,hm,hn,holt,alfs,bets,k0,res0,kvol
      common/cell2/ tres,f_inf,res1,res2,fw(NGRMX)

c variable names:
c res=resolved shear strain rate (same in walls and in cell)
c res0=reference equivalent resolved shear strain rate
c fw= volume fraction of cell walls
c f1= volume fraction at shear=0
c f2= volume fraction at a given shear
c res1, res2= resolved shear strains at f1 and f2
c f_inf=volume fraction of cell walls at infinite strain
c res_bar= constant in the exponential law
c tauw=shear stress in walls
c tauc= shear stress in cells
c tres= microscopic total resolved shear stress in a cell
c alfa=constant
c dt=time increment
c   G=shear modulus
c   b=burgers vector length
c   hm=strain rate sensitivity index
c   hn=second strain rate sensitivity index
c   row=dislocation density in walls
c   roc=dislocation density in cell interior
c
c There are 2 rate sensitive type relations for the shear stresses
c in the cells and walls. Dislocation densities are initialized, then
c incremented from the differential equations of Yuri.
c There is no need to solve the problem iteratively.

c estimated maximum dislocation density (only for time step)
      romax=1.d+16
      dtime=0.
      iend=0

c total dislocation density:   
      rotot=fw(ig)*row(ig)+(1.d0-fw(ig))*roc(ig)

c resolved shear strain rates in walls and cell interiors:
      resc=res
      resw=res

100   dt1=dt*rotot/romax
      dtime=dtime+dt1
      if(dtime.gt.dt) then
         dt1=dt-dtime
         iend=1
      end if

c total dislocation density:
      rotot=fw(ig)*row(ig)+(1.d0-fw(ig))*roc(ig)

c cell size:
      d=holt/dsqrt(rotot)

c increments in dislocation densities (3D):
      droc=(1.d0/dsqrt(3.d0))*(alfs/b)*(dsqrt(row(ig)))*resw-6.d0*bets
     1*resc/(b*d)/((1.d0-fw(ig))**(1.d0/3.d0))
     2-k0*(resc/res0)**(-1.d0/hn)*resc*roc(ig)
      drow= 6.*bets*resc*((1.d0-fw(ig))**(2./3.))/(b*fw(ig)*d)
     1+dsqrt(3.d0)*bets*resc*(1.d0-fw(ig))*dsqrt(row(ig))/(b*fw(ig))
     2-k0*(resw/res0)**(-1.d0/hn)*resw*row(ig)

      roc(ig)=roc(ig)+droc*dt1
      row(ig)=row(ig)+drow*dt1

c new resolved shear stress in walls:
      tauw=(alfa*G*b*dsqrt(row(ig))*(resw/res0)**(1./hm))/1.e+6

c new resolved shear stress in cells:
      tauc=(alfa*G*b*dsqrt(roc(ig))*(resc/res0)**(1./hm))/1.e+6

c new microscopic resolved shear stress:
      tres=fw(ig)*tauw+(1.-fw(ig))*tauc

c new volume fraction for walls:
      if(kvol.eq.1) fw(ig)=f_inf+(f1-f_inf)*dexp(-ressum(ig)/res_bar)
      if(kvol.eq.2) fw(ig)=f1+(f2-f1)*ressum(ig)/(res2-res1)
      if(kvol.eq.3) then
        res_bar=(res2-res1)*(f_inf-f2)/(f2-f1)
        if(ressum(ig).le.res2) fw(ig)=f1+(f2-f1)*ressum(ig)/(res2-res1)
        if(ressum(ig).gt.res2)
     &    fw(ig)=f_inf+(f2-f_inf)*dexp(-(ressum(ig)-res2)/res_bar)
      end if

      if(iend.eq.0) go to 100

      return
      end


c
c general subroutine for the rate sensitive solution
c
c               Laszlo S. Toth, 1991. May 1, Leuven

      subroutine rtsens(gam0,h,e,s,nc,m,tauc,nsys,igu,gam,tol,info)
	include 'POLYCR.DIM'
      implicit real*8 (a-h,o-z)
      real*8 m
      dimension indx(5)
      dimension dsi(5),f(5),fij(5,5),b(5)
      dimension e(5),s(5),m(NSYSMX,5),tauc(NSYSMX),gam(NSYSMX)

c h=1/m the value of rate sensitivity
c tol=required precision in the stress state
c tol2=precision in the strain rates (tol/100)
c max=maximum number of iterations
c succes=1 if there is convergence, succes=0 in case of no conv. (integer)
c e=prescribed strain rate vector (will be normalized here)
c s=stress vector (corresponds to unit strain rate vector)
c nc=number of prescribed components in strain rate vector
c m=slip system-orientation vectors
c iguess=0 rtsens called the first time so that previous stress solution
c          cannot be used as initial guess
c iguess=1 the previous solution can be used as initial guess
c gam=slip rates in the slip systems
c nsys=physical dimension for the maximum number of slip systems
c ncmax=physical dimension for the maximum number of unknowns (5)
c info=1: write information on Newton procedure on screen

      h0=h
      h00=h
      dh=50.d0  
      max=100
      tol2=tol/1.d+2
      jsum=0

      ntry=1

      if(igu.eq.1) go to 3000

233   call guess(e,s,m,tauc,nc,nsys)
      if(h.gt.10.d0.and.ntry.eq.1) then
         h=10.d0
         dh=10.d0
      end if

c Solve for the stress state (s) using the Newton-Raphson method

3000  jnew=0
50    continue

      call funct(gam0,m,e,s,nc,h,tauc,gam,f,fij,nsys)

      abf=0.d0
      do 20 i=1,nc
      b(i)=-f(i)
      abf=abf+dabs(f(i))
20    continue
      if(abf.lt.tol2.and.h.eq.h0) then
      jsum=jsum+jnew
      if(info.eq.1) write(*,3001) h,jsum
      go to 1000
      end if

c Solves Fij*DSi=-Fi (with LU decomposition and back substitution)
    
      call ludcmp(fij,nc,indx,d)
      call lubksb(fij,nc,indx,b)

      crit=0.d0
      do 40 i=1,nc
      dsi(i)=b(i)
      ti=dabs(dsi(i))
      if(ti.gt.0.25d0) dsi(i)=0.25d0*dsi(i)/ti
      crit=crit+ti
      s(i)=s(i)+dsi(i)
40    continue

      jnew=jnew+1

      if(jnew.gt.max) go to 100
      if(crit.gt.tol) go to 50
100   continue

      jsum=jsum+jnew
      if(info.eq.1) write(*,3001) h,jsum
3001  format('    h=',f5.1,4x,' total number of iterations=',i3/)

       if(jnew.gt.max.and.ntry.eq.1) then
       dh=50.d0 
       h=50.d0
       ntry=2
       go to 233
       end if

       if(jnew.gt.max.and.ntry.eq.2) then
       dh=20.d0  
       h=20.d0
       ntry=3
       go to 233
       end if

       if(jnew.gt.max.and.ntry.eq.3) then
       dh=10.d0
       h=10.d0
       ntry=4
       go to  233
       end if

       if(jnew.gt.max.and.ntry.eq.4) then
       dh=2.d0
       h=10.d0
       ntry=5    
       go to  233
       end if

       if(jnew.gt.max.and.ntry.eq.5) then
       dh=1.d0
       h=1.d0
       ntry=6
       go to 233
       end if

       if(jnew.gt.max.and.ntry.eq.7) then
       write(*,456)
       h=h0
       pause
       stop
       end if
456    format(/'  Newton-method does not work! Stop in RTSENS.')

       if(h.lt.h0) then
           h=h+dh
           go to 3000
       end if

       if(h.gt.h0) then
         h=h-1.d0
         go to 3000
       end if

       if(dabs(h-h0).lt.2.d0.and.dabs(h-h0).gt.1.d-6) then
         h=h0
         go to 3000
       end if

1000   h=h0
                      
c       write(*,*) jnew               
       return
       end


      subroutine funct(gam0,m,e,s,nc,h,tauc,gam,f,fij,nsys)
	include 'POLYCR.DIM'
      implicit real*8 (a-h,o-z)
      real*8 m
      dimension m(NSYSMX,5),e(5),s(5),tauc(NSYSMX),gam(NSYSMX)
      dimension f(5),fij(5,5),ee(5),dg(NSYSMX,5)

c     Calculates Fi and Fij (derivatives)

c     calculation of the shear rates in the slip systems
       do 4 i=1,nsys
       taus=0.d0
       do 5 j=1,nc
5      taus=taus+m(i,j)*s(j)
c       if(taus.lt.0.d0.and.sh(i).gt.1.d-6) taus=0.d0     !!! this is for twinning systems!!!
       gam(i)=gam0*dsign(1.d0,taus)*dabs((taus/tauc(i)))**h
       do 44 j=1,nc
       dg(i,j)=gam0*h*m(i,j)*(dabs(taus/tauc(i)))**(h-1.d0)/tauc(i)
44     continue
4      continue

c calculation of strain rate vector
       do 6 i=1,nc
       ee(i)=0.d0
       do 7 j=1,nsys
7      ee(i)=ee(i)+m(j,i)*gam(j)
6      f(i)=ee(i)-e(i)

       do 10 i=1,nc
       do 10 j=1,nc
       fij(i,j)=0.d0
       do 8 k=1,nsys
       fij(i,j)=fij(i,j)+m(k,i)*dg(k,j)
8        continue
10       continue
       return
       end

                           

      subroutine guess(e,s,m,tauc,nc,nsys)
	include 'POLYCR.DIM'

c **** provides an initial guess for the stress state ****

c e: the strain rate vector in the sample system
c s: the stress vector corresponding to h=1 (initial guess)
c nc= the number of prescribed components in the strain rate vector
c m: the orientation factors of the slip systems
c tauc: the reference ('critical') shear stresses of the slip systems
c (Reference slip rate is taken equal to 1.)
c nsys=physical dimension for the number of slip systems
c ncmax=physical dimension for the max. number of unknowns (5)

      implicit real*8 (a-h,o-z)
      real*8 m
      dimension e(5),s(5),m(NSYSMX,5),tauc(NSYSMX),a(5,5)
      dimension b(5),indx(5)

c set up the linear equation system

      do 1 i=1,nc
      do 1 j=1,nc
      a(i,j)=0.d0
      do 1 L=1,nsys  
1     a(i,j)=a(i,j)+m(L,i)*m(L,j)/tauc(L)

      do 2 i=1,nc
2     b(i)=e(i)

c call the linear equation solver (LU decomp. and back substitution)

        call ludcmp(a,nc,indx,d)          
        call lubksb(a,nc,indx,b)
 
      big=dabs(b(1))
      do 90 i=2,nc
      if(dabs(b(i)).gt.big) then
      big=dabs(b(i))
      end if
90    continue
      do 3 i=1,nc
3     s(i)=b(i)/big
 
      return
      end



c *******************************************

      SUBROUTINE LUDCMP(A,N,INDX,D)
      implicit real*8(a-h,o-z)
      PARAMETER (NMAX=5,TINY=1.0d-20)
      DIMENSION A(5,5),INDX(N),VV(NMAX)
      D=1.d0
      DO 12 I=1,N
      AAMAX=0.d0
      DO 11 J=1,N
        IF (DABS(A(I,J)).GT.AAMAX) AAMAX=DABS(A(I,J))
11      CONTINUE
      IF (AAMAX.EQ.0.d0) PAUSE 'Singular matrix.'
      VV(I)=1.d0/AAMAX
12    CONTINUE
      DO 19 J=1,N
      IF (J.GT.1) THEN
        DO 14 I=1,J-1
          SUM=A(I,J)
          IF (I.GT.1)THEN
            DO 13 K=1,I-1
       SUM=SUM-A(I,K)*A(K,J)
13            CONTINUE
            A(I,J)=SUM
          ENDIF
14        CONTINUE
      ENDIF
      AAMAX=0.d0
      DO 16 I=J,N
        SUM=A(I,J)
        IF (J.GT.1)THEN
          DO 15 K=1,J-1
            SUM=SUM-A(I,K)*A(K,J)
15          CONTINUE
          A(I,J)=SUM
        ENDIF
        DUM=VV(I)*DABS(SUM)
        IF (DUM.GE.AAMAX) THEN
          IMAX=I
          AAMAX=DUM
        ENDIF
16      CONTINUE
      IF (J.NE.IMAX)THEN
        DO 17 K=1,N
          DUM=A(IMAX,K)
          A(IMAX,K)=A(J,K)
          A(J,K)=DUM
17        CONTINUE
        D=-D
        VV(IMAX)=VV(J)
      ENDIF
      INDX(J)=IMAX
      IF(J.NE.N)THEN
        IF(A(J,J).EQ.0.d0) A(J,J)=TINY
        DUM=1.d0/A(J,J)
        DO 18 I=J+1,N
          A(I,J)=A(I,J)*DUM
18        CONTINUE
      ENDIF
19    CONTINUE
      IF(A(N,N).EQ.0.d0) A(N,N)=TINY
      RETURN
      END


      subroutine LUBKSB(A,N,INDX,B)
      implicit real*8(a-h,o-z)
      DIMENSION A(5,5),INDX(N),B(N)
      II=0
      DO 12 I=1,N
      LL=INDX(I)
      SUM=B(LL)
      B(LL)=B(I)
      IF (II.NE.0) THEN
        DO 11 J=II,I-1
          SUM=SUM-A(I,J)*B(J)
11        CONTINUE
      ELSE IF (SUM.NE.0.d0) THEN
        II=I
      ENDIF
      B(I)=SUM
12    CONTINUE
      DO 14 I=N,1,-1
      SUM=B(I)
      IF(I.LT.N)THEN
        DO 13 J=I+1,N
          SUM=SUM-A(I,J)*B(J)
13        CONTINUE
      ENDIF
      B(I)=SUM/A(I,I)
14    CONTINUE
      return
      END

                           