!***********************************************************************
!*****               PIC ExB ion motion Simulator                  *****
!*****           Developed by Gen Ito  (Univ of Tokyo)             *****
!*****                    Last Update 01/31/2017                   *****
!***********************************************************************


program main

   implicit none
   
   !Program parameter
   character(len=80)             :: TOPDIR= '/home/ito/'   !Directory for input and output
   !Program parameter
   integer         ,parameter    :: ITM   = 20000                !Maximum number of time steps
   double precision,parameter    :: DTPIC = 5.0d-4  ![s] Time step
   !Program Switch (0:No, 1:Yes)
   integer         ,parameter    :: OINF   = 0                           !Output time step information?
   integer         ,parameter    :: OINFL  = 1                         !Output last time step information?
  
   
   !Physical parameter
   double precision,parameter    :: Q       = 1.0d0 !1.60217733d-19            !Elementary charge[c]
   double precision, parameter   :: ME      = 1.0d0 !9.1093826d-31             !Electron mass[kg]
   !Particle parameter
   integer         ,parameter    :: NMAX  = 1                           !Maximum number of particles
   
   integer :: m,nm=0,it
   double precision :: ttotal,bb,vx0,vy0,vz0,vv0,rL,wc,x0,y0,z0,vxh,vyh,vzh
   double precision,dimension(NMAX) :: &
      x,y,z,vx,vy,vz,x_post,y_post,z_post,vx_post,vy_post,vz_post, &
      ex,ey,ez,bx,by,bz
   
   
   nm = 1
   
   do m = 1,nm
      !electric field [V]
      ex(m) = 1.0d0
      ey(m) = 0.0d0
      ez(m) = 0.0d0
      !magnetic field [T]
      bx(m) = 0.0d0
      by(m) = 0.0d0
      bz(m) = 1.0d0
      bb = dsqrt(bx(m)**2.0d0+by(m)**2.0d0+bz(m)**2.0d0)
      !initial velocity [m/s]
      vx0 = 10.0d0
      vy0 = 0.0d0
      vz0 = 0.0d0 
      vv0 = dsqrt(vx0**2.0d0+vy0**2.0d0+vz0**2.0d0)
      !Larmor radius
      rL = ME*vv0/Q/bb
      !Larmor radius
      wc = vv0/rL  
      !initial 1/2t velocity [m/s] (when Ex,Bz only)
      vxh = vv0*dcos(wc*DTPIC*0.5d0)
      vyh = -vv0*dsin(wc*DTPIC*0.5d0)-ex(m)/bz(m)
      vzh = vz0
      !velocity of first time step[m/s] 
      vx(m) = vxh
      vy(m) = vyh
      vz(m) = vzh
      
      !initial position [m]
      x0 = rL*dsin(wc*DTPIC)
      y0 = rL*dcos(wc*DTPIC)-ex(m)/bz(m)*DTPIC
      z0 = 0.0d0
      
      !position of first time step [m]
      x(m) = x0
      y(m) = y0
      z(m) = z0
      
   enddo
   
   !Output initial 1t particle information
    if(OINFL.eq.1) then
       write(*,*) 'x,y,z,vx,vy,vz'
       do m = 1,nm
          write(*,'(11E15.5)') x(m),y(m),z(m),vx(m),vy(m),vz(m),wc,rL
       enddo
     endif
    
   
    open(unit=22,file=trim(TOPDIR)//'PIC'//'.dat',form='formatted',status='replace')
    close(22)   
    
   
   do it=2,ITM
      do m = 1,nm
         Call BunemanBoris(NMAX,nm,Q,ME,DTPIC, &
            x(m),y(m),z(m),vx(m),vy(m),vz(m),ex(m),ey(m),ez(m),bx(m),by(m),bz(m), &
            x_post(m),y_post(m),z_post(m),vx_post(m),vy_post(m),vz_post(m))

            
         !update position
            
            x(m) = x_post(m)
            y(m) = y_post(m)
            z(m) = z_post(m)
            
            vx(m) = vx_post(m)
            vy(m) = vy_post(m)
            vz(m) = vz_post(m)

            !Output particle information
            if(OINF.eq.1) then
              write(*,'(11E15.5)') x(m),y(m),z(m),vx(m),vy(m),vz(m)
            endif
            
            open(unit=22,file=trim(TOPDIR)//'PIC'//'.dat',form='formatted',status='old',position='append')
                 
              
                        write(22,'(4E15.5)') &
                           x(m),y(m),vx(m),vy(m)
              
               
            close(22)
      enddo
   enddo

   ttotal = ITM*DTPIC
   !Output last time step particle information
   
   
   

   



            
   if(OINFL.eq.1) then
   write(*,*) 'DTPIC,x,y,z,vx,vy,vz,tt'
     write(*,'(18E17.10)') DTPIC,x(nm),y(nm),ttotal
   endif

   write(*,*) 'Preset time step is completed...'
   write(*,*) 'Program end...'
   write(*,'(A20,E12.3,A)') '         Time =',ttotal          ,' [s]'

      
 stop
 end program main




subroutine BunemanBoris(NMAX,nm,q,mass,dt, &
   x_pre,y_pre,z_pre,vx_hpre,vy_hpre,vz_hpre,ex_pre,ey_pre,ez_pre, &
   bx_pre,by_pre,bz_pre,x_post,y_post,z_post,vx_hpost,vy_hpost,vz_hpost)

!reference
!http://www.astro.phys.s.chiba-u.ac.jp/pcans/algorithm.html
!http://center.stelab.nagoya-u.ac.jp/summer-school/pdf/text6.pdfhttp://center.stelab.nagoya-u.ac.jp/summer-school/pdf/text6.pdf
!ignore Lorentz factor

   implicit none
   double precision, intent(in)              :: q,mass,dt
   integer                                   :: m,nm,NMAX
   double precision                          :: mq,ex,ey,ez,bx,by,bz,x,y,z
   double precision                          :: vx0,vxm,vxs,vxp
   double precision                          :: vy0,vym,vys,vyp
   double precision                          :: vz0,vzm,vzs,vzp
   double precision                          :: ttx,tty,ttz,ssx,ssy,ssz
   !double precision                          :: rp                      !Lorentz factor
   !double precision                          :: cc                
   double precision,dimension(NMAX),intent(in) :: &
      ex_pre,ey_pre,ez_pre,bx_pre,by_pre,bz_pre,x_pre,y_pre,z_pre,vx_hpre,vy_hpre,vz_hpre
   double precision,dimension(NMAX),intent(out) :: &
      x_post,y_post,z_post,vx_hpost,vy_hpost,vz_hpost

   !cc =  2.99792458d8                                                   !speed of light[m/s]  
      
do m = 1,nm
      mq = q/mass
      x  = x_pre(m)
      y  = y_pre(m)
      z  = z_pre(m) 
      vx0 = vx_hpre(m)
      vy0 = vy_hpre(m)
      vz0 = vz_hpre(m)
      ex = ex_pre(m)
      ey = ey_pre(m)
      ez = ez_pre(m)
      bx = bx_pre(m)
      by = by_pre(m)
      bz = bz_pre(m)
      !rp = dsqrt(1.0d0+up2/cc**2.0d0)                                   !Lorentz factor
      ttx = mq*dt*bx/2.0d0
      tty = mq*dt*by/2.0d0
      ttz = mq*dt*bz/2.0d0
      ssx = 2.0d0*ttx/(1.0d0+ttx**2.0d0+tty**2.0d0+ttz**2.0d0)
      ssy = 2.0d0*tty/(1.0d0+ttx**2.0d0+tty**2.0d0+ttz**2.0d0)
      ssz = 2.0d0*ttz/(1.0d0+ttx**2.0d0+tty**2.0d0+ttz**2.0d0)
! acceralation by electric field
      vxm = vx0+0.5d0*mq*dt*ex
      vym = vy0+0.5d0*mq*dt*ey
      vzm = vz0+0.5d0*mq*dt*ez
! acceralation by magnetic field
      vxs = vxm+(vym*ttz-vzm*tty)
      vys = vym+(vzm*ttx-vxm*ttz)
      vzs = vzm+(vxm*tty-vym*ttx)
      
      vxp = vxm+(vys*ssz-vzs*ssy)
      vyp = vym+(vzs*ssx-vxs*ssz)
      vzp = vzm+(vxs*ssy-vys*ssx)
! acceralation by electric field      
      vx_hpost(m) = vxp+0.5d0*mq*dt*ex
      vy_hpost(m) = vyp+0.5d0*mq*dt*ey
      vz_hpost(m) = vzp+0.5d0*mq*dt*ez
! move position
      x_post(m)  = x + vx_hpost(m)*dt
      y_post(m)  = y + vy_hpost(m)*dt
      z_post(m)  = z + vz_hpost(m)*dt

enddo

   return
endsubroutine
