
c
c   deformation texture simulation based on rate sensitive slip
c   governing program for polycrystal deformation

	program POLYCR  
      include 'POLYCR.DIM'

      implicit real*8 (a-h,o-z)

      integer*2 nab   
c nab(48,48) hardening matrix      

      character*20 file1,file2
      character*40 title

      real*8 mc,k0

c grain orientations (three Euler angles and volume of each grain):
      dimension fi1(NGRMX),fi(NGRMX),fi2(NGRMX),f(NGRMX)
	dimension fi10(NGRMX),fi0(NGRMX),fi20(NGRMX)


c reference strengths of slip systems, tauc0.  
      dimension tauc0(NGRMX,NSYSMX),tau0(NSYSMX)
                  
      dimension vel(3,3),irel(3,3),avfgrad(3,3)
      dimension mc(3,3,NSYSMX)
      dimension nab(NSYSMX,NSYSMX),omega(1000) 
      common/cell1/ res,ressum(NGRMX),row(NGRMX),roc(NGRMX),f1,
     &f2,res_bar,alfa,G,b,hm,hn,holt,alfs,bets,k0,res0,kvol
      common/cell2/tres,f_inf,res1,res2,fw(NGRMX)
c /cell1/,/cell2/: hardening parameters for dislocation-cell type hardening
      common/bronk/ q(4),h0,tsat
c /bronk/: hardening parameters for Bronkhorst type of hardening
c /slip1/: common block for slip distribution
      common/slip1/ slip(NGRMX,NSYSMX)
	common/velcr/ velg(NGRMX,3,3),fgrad(NGRMX,3,3)			   ! velg: grain velocity gradient made by slip, fgrad: finite grain strain gradient
	common/strs/ st(NGRMX,3,3),avst(3,3),defg(NGRMX,3,3),avdef(3,3),avvelg(3,3)   
	common/trace/ sn(NSYSMX,3),beta(NGRMX,NSYSMX)  ! slip plane normals in crystal reference system, not normalized, and the orientation of the slip line with respect to the x-axis (the surface plane must be x-y)
	  
c output listing file:
      open(3,file='polycr.lst',status='unknown')
      open(9,file='grains.lst',status='unknown')

c read input parameter file: polycr.ctl
      open(2,file='polycr.ctl',status='old')
      read(2,1) file1      ! name of input file containing discrete distr.
1     format(a)
      write(3,'(a)') ' input distribution file: ',file1
      write(9,'(a)') ' input distribution file: ',file1
      read(2,1) file2      ! name of slip systems file
      write(3,'(a)') ' name of slip system file: ',file2
      write(9,'(a)') ' name of slip system file: ',file2
      if(file2.eq.'cubicsys.dat') go to 844
      if(file2.eq.'hcpsyst.dat') go to 844
      write(*,'(a)') ' input slip system file must be:  cubicsys.dat
     & or hcpsyst.dat!'
      stop
844   continue
      read(2,1) title      ! title of run
      write(3,'(a)') ' title of this run: ',title
      write(9,'(a)') ' title of this run: ',title
      read(2,*) h    ! value of rate sensitivity h=1/m
      write(3,814) h
      write(9,814) h
814   format(' value of rate sensitivity (1/m): ',f10.5)
      read(2,*) Lankford      ! Lankford parameter (R value) simulation) (sheet normal is axis 3)
	if(Lankford.eq.0) read(2,*) 
	if(Lankford.eq.0) read(2,*)
	if(Lankford.eq.1) read(2,*) NLank 
	if(Lankford.eq.1) read(2,*) (omega(i),i=1,NLank) 
      read(2,*) nhard
c nhard=2 isotropic cell hardening is activated
c nhard=1 Bronkhorst latent  hardening
                           ! hardening file name: hard.dat
c nhard=0 no hardening
      write(3,'(a,i1)') ' hardening parameter: ',nhard
      write(9,'(a,i1)') ' hardening parameter: ',nhard
      read(2,*) tol
      write(3,816) tol       ! required precision in the stress state
      write(9,816) tol       ! required precision in the stress state
816   format(' required precision in the stress state: ',f16.14)
      read(2,*) nrun       ! number of runs
      write(3,'(a,i5)') ' number of increments: ',nrun
      write(9,'(a,i5)') ' number of increments: ',nrun
      read(2,*) deps      ! von Mises strain increment in one step
      write(3,'(a,e14.5)') ' von Mises strain increment in one step: ',deps
      write(9,'(a,e14.5)') ' von Mises strain increment in one step: ',deps
	read(2,*) vel(1,1),vel(1,2),vel(1,3)
	read(2,*) vel(2,1),vel(2,2),vel(2,3)
	read(2,*) vel(3,1),vel(3,2),vel(3,3)
	if(Lankford.eq.1) then
c	vel=0.
	vel(1,1)=1.d0	   ! tensile strain in direction 1
	end if
c read the relaxed condition in a matrix form: if 1:relaxed
      read(2,4) ((irel(i,j),j=1,3),i=1,3)
4     format(3i5/3i5/3i5)
      write(3,'(a)') ' relaxation matrix:'
      write(3,4) ((irel(i,j),j=1,3),i=1,3)
      write(9,'(a)') ' relaxation matrix:'
      write(9,4) ((irel(i,j),j=1,3),i=1,3)     
	close (2)

	open(9,file='grains.lst',status='unknown')
	open(7,file='Lankford.dat',status='unknown')

      
c reads the parameter file for the dislocation cell based hardening case
      if(nhard.eq.2) call dislhard

c produces the slip system schmid tensors (mc):
      call miller(file2,nsys,mc,tau0,
     & nab,h0,tsat,a,q,gam0)
c nab: hardening matrix

      fgrad=0.d0
c read input texture data file for normalizing:
      fsum=0.d0
      open(2,file=file1,status='old')

	read(2,'(a)')  title
	read(2,*)
	read(2,*)
	read(2,'(1x,i12)') ngrain
      do 100 ig=1,ngrain
      read(2,*) fi10(ig),fi0(ig),fi20(ig),f(ig)
      do j=1,nsys
	tauc0(ig,j)=tau0(j)
	end do
	do i=1,3
	fgrad(ig,i,i)=1.d0
	end do
      fsum=fsum+f(ig)
100   continue
      fi1=fi10
	fi=fi0		 ! texture for normal simulation (not Lankford)
	fi2=fi20

	fgrad(ig,i,i)=1.d0

c*****************************************************************
c Lankford simulation:
      if(Lankford.eq.0) goto 299
	! loop for the angle from axis 1
	dt=1.d-2  ! time increment

	do 48 iangle=1,NLank
	angle=omega(iangle)
	! rotation of the texture by alfa
	do ig=1,ngrain
	fi1(ig)=fi10(ig)+angle
	end do
	fi=fi0
	fi2=fi20
      igu=0  ! no initial guess for stress states in grains

      call increm(ngrain,nsys,mc,nab,h,gam0,nhard,irel,vel,fi1,fi,fi2,f,tauc0,igu,dt)
      write(3,'(a,f10.3)') ' angle: ',angle
      write(3,'(a,2f10.3)') ' angle and R value: ',angle,(avdef(2,2)/avdef(3,3))
      write(7,'(2f10.3)') angle,(avdef(2,2)/avdef(3,3))
      write(*,'(a,2f10.3)') ' angle and R value: ',angle,(avdef(2,2)/avdef(3,3))
	
      avfgrad=0.d0
	do ig=1,ngrain
	do i=1,3
	do j=1,3
	avfgrad(i,j)=avfgrad(i,j)+fgrad(ig,i,j)*f(ig)/fsum
	end do
	end do
	end do

      write(3,*) ' average strain rate:'
      do i=1,3     
	write(3,'(3f10.5)') (avdef(i,j),j=1,3)
      end do
	write(3,*) ' average slip-made velocity gradient:'
      do i=1,3     
	write(3,'(3f10.5)') (avvelg(i,j),j=1,3)
      end do
	write(3,*) ' average strain gradient:'
      do i=1,3     
	write(3,'(3f10.5)') (avfgrad(i,j),j=1,3)
      end do
	write(3,*) ' average stress state: '
	do i=1,3     
      write(3,'(3f10.5)') (avst(i,j),j=1,3)
      end do

48    continue
	close(3)
	close(7)
	stop
 
      
299   continue
c*********************************************************************
c loop for the increments:

      do 4000 inc=1,nrun               
      write(*,'(a,i5)') ' increment :  ',inc

c von Mises eqv. strain rate:
      vMs=0.d0
	do i=1,3
	do j=1,3
	vMs=vMs+((vel(i,j)+vel(j,i))/2.d0)**2
	end do
	end do
	vMs=dsqrt(2.d0*vMs/3.d0)
c the time increment:
      dt=deps/vMs

c call one increment strain rate sensitive polycrystal routine
      if(inc.eq.1) then 
        igu=0    ! no initial stress state
        else
        igu=1
      end if 

      call increm(ngrain,nsys,mc,nab,h,gam0,nhard,irel,vel,fi1,fi,fi2,f,tauc0,igu,dt)
     
      write(3,'(a,i5)') ' increment: ',inc
      write(9,'(a,i5)') ' increment: ',inc
	
      avfgrad=0.d0
	do ig=1,ngrain
	do i=1,3
	do j=1,3
	avfgrad(i,j)=avfgrad(i,j)+fgrad(ig,i,j)*f(ig)/fsum
	end do
	end do
	end do

      write(3,*) ' average strain rate:'
      do i=1,3     
	write(3,'(3f10.5)') (avdef(i,j),j=1,3)
      end do
	write(3,*) ' average slip-made velocity gradient:'
      do i=1,3     
	write(3,'(3f10.5)') (avvelg(i,j),j=1,3)
      end do
	write(3,*) ' average strain gradient:'
      do i=1,3     
	write(3,'(3f10.5)') (avfgrad(i,j),j=1,3)
      end do
	write(3,*) ' average stress state: '
	do i=1,3     
      write(3,'(3f10.5)') (avst(i,j),j=1,3)
      end do

c ekv. von Mises strain rate
     	vonmis=0.d0
	do i=1,3
	do j=1,3
	vonmis=vonmis+avdef(i,j)*avdef(i,j)
	end do
	end do
	vonmis=dsqrt(2.d0*vonmis/3.d0)
	rslip=0.d0	 ! total slip

c print data on grains into grains.lst file
      do ig=1,ngrain
	write(9,'(a,i5)') 'grain index:',ig
	write(9,*)  ' velocity gradient made by slip:'
	do i=1,3
	write(9,'(3f10.6)') (velg(ig,i,j),j=1,3)
	end do
	write(9,*)  ' grain strain gradient:'
	do i=1,3
	write(9,'(3f10.6)') (fgrad(ig,i,j),j=1,3)
	end do
	write(9,*)  ' Cauchy stress tensor:'
	do i=1,3
	write(9,'(3f10.4)') (st(ig,i,j),j=1,3)
	end do
      write(9,*)  ' slip rates in slip systems, and slip line orientation:'
	do j=1,nsys
      write(9,'(i3,2f10.5)') j,slip(ig,j),beta(ig,j)
	rslip=rslip+f(ig)*dabs(slip(ig,j))
	end do
	write(9,*) ' Euler angles:'
	write(9,'(3f10.3)')  fi1(ig),fi(ig),fi2(ig)
	write(9,*)
	end do

c Taylor factor
      taylor=rslip/vonmis/fsum
	write(*,*) ' Taylor factor on polycyrstal:  ',taylor
	write(3,*) ' Taylor factor on polycyrstal:  ',taylor
4000  continue
c end of increment-loop

c write into output texture files
5000  continue
      open(4,file='temp.smt',status='unknown')   
c temp.smt: ready for ODF program of Paul Van Houtte              
	  write(4,1) title
	  write(4,*)
	  write(4,*)	
      write(4,'(a,i10)') 'B',ngrain 	
      do i=1,ngrain                 
      xfi1=fi1(i)
      xfi=fi(i)
      xfi2=fi2(i)
         if(xfi1.lt.0.d0) xfi1=xfi1+360.d0
         if(xfi1.lt.0.d0) xfi1=xfi1+360.d0
         if(xfi.lt.0.d0) xfi=xfi+360.d0
         if(xfi.lt.0.d0) xfi=xfi+360.d0
         if(xfi2.lt.0.d0) xfi2=xfi2+360.d0
         if(xfi2.lt.0.d0) xfi2=xfi2+360.d0
      write(4,'(4f12.3)') xfi1,xfi,xfi2,f(i)
      end do
301   format(4f10.5,i5,5X,f10.5,f3.1)  
      close (2)
      close (3)
      close (4)
      
      end

c end of governing main program for ECAE



c ****************************
c Produces the Schmid tensors of the slip systems and the hardening matrix

       subroutine miller(file2,isyss,mc,tauc0,
     & nab,h0,tsat,a,q,gam0)
      include 'POLYCR.DIM'
       implicit real*8 (a-h,o-z)
       real*8 mc
       integer*2 nab
       character*16 title
       character*30 sys
       character*20 file2
       dimension tauc0(NSYSMX),nn(NSYSMX,4),nb(NSYSMX,4),mc(3,3,NSYSMX)
       dimension ifam1(10),ifam2(10),q(4),nab(NSYSMX,NSYSMX),sh(NSYSMX)
       dimension sb(NSYSMX,3),ssni(3),ssnj(3),ssbi(3),ssbj(3)
	common/trace/ sn(NSYSMX,3),beta(NGRMX,NSYSMX)  ! slip plane normals in crystal reference system, not normalized, and the orientation of the slip line with respect to the x-axis (the surface plane must be x-y)

       sq3=dsqrt(3.d0)
       if(file2.eq.'cubicsys.dat') then
       open(12,file='cubicsys.dat',status='old')
       ind=3
       end if
       if(file2.eq.'hcpsyst.dat') then
       open(12,file='hcpsyst.dat',status='old')
       ind=4
       end if
       open(13,file='planes.dat',status='unknown')
       open(14,file='systems.lst',status='unknown')
       write(14,55)
55     format(/'****** List of the selected slip sytems'//)

       nfam=0   ! number of system families
1      format(A16)
4      format(A30,i2)
60     format(/
     !'  There is an error in the definition of slip system nr. ',i2)
61     format(/
     !'  ******* Slip system nr ',i2,' is not orthogonal! *******')

       if(ind.eq.3) then
       read(12,1) title
       read(12,*)
       read(12,*)
       end if
       if(ind.eq.4) then
       read(12,1) title
       read(12,'(17x,f10.5)') ca  ! c/a ratio for hexagonal
       read(12,*)
       end if

       read(12,*) gam0,h0,tsat,a
       read(12,*) q
c q(1): hardening parameter for coplanar slip systems
c q(2): hardening parameter for collinear slip systems
c q(3): hardening parameter for slip perpendicular to the other plane
c q(4): hardening parameter for all other
       read(12,*)

       isyss=0
30     read(12,4) sys,nsys
       if(sys.eq.'end of systems') go to 100
       read(12,*)
       read(12,*) crss,ssh
       if(crss.gt.1.d-6) then
          nfam=nfam+1
          ifam1(nfam)=isyss+1
          ifam2(nfam)=ifam1(nfam)+nsys-1
          do 20 i=(isyss+1),(isyss+nsys)
          tauc0(i)=crss
          sh(i)=ssh
          read(12,*) in,(nn(i,j),j=1,ind),(nb(i,j),j=1,ind)
       if(ind.eq.3)
     & write(14,75) i,in,(nn(i,j),j=1,ind),(nb(i,j),j=1,ind),tauc0(i)
       if(ind.eq.4)
     & write(14,76) i,in,(nn(i,j),j=1,ind),(nb(i,j),j=1,ind),tauc0(i)
       nerror=0
       do j=1,ind
       nerror=nerror+nn(i,j)*nb(i,j)
       end do
75     format(2i3,'  (',i2,i3,i3,')  [',i2,i3,i3']  tauc: ',f10.5)
76     format(2i3,'  (',i2,i3,i3,i3')  [',i2,i3,i3,i3']  tauc: ',f10.5)
          if(nerror.ne.0) then
             write(*,60) i
             stop
          end if
20     continue
       else
          do 21 j=1,nsys
21        read(12,*)
       end if
       read(12,*)
       if(crss.gt.1.d-6) isyss=isyss+nsys
       go to 30
100    continue
       close (12)

c check for normality:

       do 80 i=1,isyss
       nerror=0
       do j=1,ind
       nerror=nerror+nn(i,j)*nb(i,j)
       end do
       if(nerror.ne.0) then
          write(*,61) i
          stop
       end if
80     continue

c conversion of hexagonal Miller-Bravais indexes to Miller:
       if(ind.eq.4) then
       do i=1,isyss
       sb(i,1)=sq3*dble(nb(i,1)-nb(i,3))/2.d0
       sb(i,2)=3.d0*dble(nb(i,2))/2.d0
       sb(i,3)=ca*dble(nb(i,4))
       sn(i,1)=(2.d0*dble(nn(i,1))+dble(nn(i,2)))/sq3
       sn(i,2)=dble(nn(i,2))
       sn(i,3)=dble(nn(i,4))/ca
       end do
       end if

c The orientation factors:
      if(ind.eq.3) then
      do 235 i=1,isyss
c length of slip plane normal and slip direction vectors:
      vb=dsqrt(dble(nb(i,1))**2+dble(nb(i,2))**2+dble(nb(i,3))**2)
      vn=dsqrt(dble(nn(i,1))**2+dble(nn(i,2))**2+dble(nn(i,3))**2)
      do 235 j=1,3
      sn(i,j)=dble(nn(i,j))/vn
      do 235 k=1,3
235   mc(j,k,i)=dble(nb(i,j)*nn(i,k))/vb/vn
      end if

c The orientation factors:
      if(ind.eq.4) then
      do 236 i=1,isyss
c length of slip plane normal and slip direction vectors:
      vb=dsqrt(sb(i,1)**2+sb(i,2)**2+sb(i,3)**2)
      vn=dsqrt(sn(i,1)**2+sn(i,2)**2+sn(i,3)**2)
      do 236 j=1,3
      do 236 k=1,3
236   mc(j,k,i)=sb(i,j)*sn(i,k)/vb/vn
      end if

c Construction of the hardening matrix:
      do 450 i=1,isyss
      do 450 j=1,isyss
      nab(i,j)=0
c normalisation:
      if(ind.eq.3) then
      sni=dsqrt(dble(nn(i,1)*nn(i,1)+nn(i,2)*nn(i,2)+nn(i,3)*nn(i,3)))
      snj=dsqrt(dble(nn(j,1)*nn(j,1)+nn(j,2)*nn(j,2)+nn(j,3)*nn(j,3)))
      sbi=dsqrt(dble(nb(i,1)*nb(i,1)+nb(i,2)*nb(i,2)+nb(i,3)*nb(i,3)))
      sbj=dsqrt(dble(nb(j,1)*nb(j,1)+nb(j,2)*nb(j,2)+nb(j,3)*nb(j,3)))
      do 324 k=1,3
      ssni(k)=dble(nn(i,k))/sni
      ssnj(k)=dble(nn(j,k))/snj
      ssbi(k)=dble(nb(i,k))/sbi
324   ssbj(k)=dble(nb(j,k))/sbj
      end if
      if(ind.eq.4) then
      sni=dsqrt(sn(i,1)*sn(i,1)+sn(i,2)*sn(i,2)+sn(i,3)*sn(i,3))
      snj=dsqrt(sn(j,1)*sn(j,1)+sn(j,2)*sn(j,2)+sn(j,3)*sn(j,3))
      sbi=dsqrt(sb(i,1)*sb(i,1)+sb(i,2)*sb(i,2)+sb(i,3)*sb(i,3))
      sbj=dsqrt(sb(j,1)*sb(j,1)+sb(j,2)*sb(j,2)+sb(j,3)*sb(j,3))
      do 325 k=1,3
      ssni(k)=sn(i,k)/sni
      ssnj(k)=sn(j,k)/snj
      ssbi(k)=sb(i,k)/sbi
325   ssbj(k)=sb(j,k)/sbj
      end if

c check for coplanar situation:
      co=ssni(1)*ssnj(1)+ssni(2)*ssnj(2)+ssni(3)*ssnj(3)
      if(dabs(1.d0-dabs(co)).lt.1.d-5) nab(i,j)=1
      if(nab(i,j).ne.0) go to 450
c check for collinear situation:
      co=ssbi(1)*ssbj(1)+ssbi(2)*ssbj(2)+ssbi(3)*ssbj(3)
      if(dabs(1.d0-dabs(co)).lt.1.d-5) nab(i,j)=2
      if(nab(i,j).ne.0) go to 450
c check for perpendicular situation
      co=ssni(1)*ssbj(1)+ssni(2)*ssbj(2)+ssni(3)*ssbj(3)
      if(dabs(co).lt.1.d-5) nab(i,j)=3
      co=ssnj(1)*ssbi(1)+ssnj(2)*ssbi(2)+ssnj(3)*ssbi(3)
      if(dabs(co).lt.1.d-5) nab(i,j)=3
450   if(nab(i,j).eq.0) nab(i,j)=4

c write output file:

      write(*,'(a,i5,a)') ' Numb. of systems:',isyss,'- LSTAYLOR begins'
      write(13,766) isyss
766   format(i5)
      do 743 i=1,isyss
743   write(13,755) i,tauc0(i),((mc(j,k,i),k=1,3),j=1,3)
754   format(2x,i3,4x,f7.4,3x,9(f5.2,2x))
755   format(2x,i3,4x,f7.4,3x,9(f15.12))
      write(13,'(a,i5)') ' number of familie(s): ',nfam
      do 102 i=1,nfam
102   write(13,'(2i5)') ifam1(i),ifam2(i)
      write(14,*)
      write(14,'(a)') ' Hardening parameters:'
      write(14,'(a,f10.4)') ' h0=    ',h0
      write(14,'(a,f10.4)') ' tausat=',tsat
      write(14,'(a,f10.4)') ' a=     ',a
      write(14,*)
      write(14,'(a)') ' Hardening matrix qab:'
      write(14,*)
      write(14,'(a,f7.4)') 'coplanar hard. factor q1:     ',q(1)
      write(14,'(a,f7.4)') 'collinear hard. factor q2:    ',q(2)
      write(14,'(a,f7.4)') 'perpendicular case factor q3: ',q(3)
      write(14,'(a,f7.4)') 'for all other cases q4:       ',q(4)
      write(14,*)
      do 678 i=1,isyss
678   write(14,'(48i2)') (nab(i,j),j=1,isyss)

       close (13)
       close (14)
       return
       end

c************************************
c reads the parameter file for the dislocation cell based hardening case

      subroutine dislhard
	include 'POLYCR.DIM'
      implicit real*8(a-h,o-z)
      real*8 k0
      common/cell1/ res,ressum(NGRMX),row(NGRMX),roc(NGRMX),f1,f2
	&,res_bar,alfa,G,b,hm,hn,holt,alfs,bets,k0,res0,kvol
      common/cell2/tres,f_inf,res1,res2,fw(NGRMX)

c variable names:
c res=resolved shear strain rate (same in walls and in cell)
c res0=reference equivalent resolved shear strain rate
c f= volume fraction of cell walls
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
      open(11,file='cell.ctl',status='old')
      read(11,*) row0
      read(11,*) roc0
      read(11,*) kvol
      read(11,*) f_inf
      read(11,*) res_bar
      read(11,*) f1,res1
      read(11,*) f2,res2
      read(11,*) hm
      read(11,*) hn
      read(11,*) res0
      read(11,*) alfa
      read(11,*) G
      read(11,*) b
      read(11,*) holt
      read(11,*) alfs   ! alfa-star
      read(11,*) bets   ! beta-star
      read(11,*) k0
      close (11)
      return
      end

c   