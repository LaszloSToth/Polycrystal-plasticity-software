texin.smt                     ! input texture file                   
cubicsys.dat                 ! slip system file             
building the program                 ! title of run
  100.                    ! rate sensitivity (1/m)
  1                        ! if 1: Lankford parameter (R value) simulation, sheet normal is axis 3
  3 								 ! If Lankford simulation, for how many positions?
  0 45 90  							 ! list of angles where Lankford is calculated
0                           !0: no hardening, 1:Bronkhorst hard.,2: isotropic cell hardening is activated
1.d-8                        ! required precision in the stress state           
1                           ! number of increments                              
 0.01       increment of strain in von Mises eqv. strain
 1.  0.  0.		  velocity gradient
 0.  -0.5  0.
 0.  0   -0.5
    0    1    0                ! relaxation matrix
    1    1    0
    0    0    1
