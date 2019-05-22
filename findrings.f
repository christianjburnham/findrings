!********************************************************
!**    FIND RINGS
!**
!**    Ring recognition code
!**    Christian J Burnham, UCD, 2019
!**    @christianjburnham@gmail.com
!********************************************************

      implicit none
      real(8), dimension(3,3) :: cvec,cveci,cvec0,cvec0i
      real(8) :: xdif,ydif,zdif,rcut2,rr2
      real(8), dimension(80000,3) :: rox,rhy,rr
      integer natoms,nox,nhy
      integer, dimension(100) :: ringcount
      real(8), dimension(80000,3,3) :: rwater
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      integer, dimension(10,80000,0:10) :: ringlist
      integer :: maxring
      integer ipdc,ncell
      real(8) :: radius
      character(len = 32) :: input_file_name

      open(14,file = 'findrings.control')
      read(14,*) input_file_name
      read(14,*) radius

      open(10,file = input_file_name)

!     periodic boundary conditions, 1 = on, 0 = off
      ipdc = 1
!     maximum number of rings used in the analysis
      maxring = 10

      write(*,*) '********************&
     & ************************************'
      write(*,*) '**    FIND RINGS'
      write(*,*) '**    ring recognition code'
      write(*,*) '**    Christian J Burnham, UCD, 2019'
      write(*,*) '**    @christianjburnham@gmail.com'
      write(*,*) '********************&
     &************************************'

      open(23,file = 'ring3.dat')
      open(24,file = 'ring4.dat')
      open(25,file = 'ring5.dat')
      open(26,file = 'ring6.dat')
      open(27,file = 'ring7.dat')
      open(28,file = 'ring8.dat')
      open(29,file = 'ring9.dat')
      open(30,file = 'ring10.dat')

!     read in the coordinates from the file
      call readcoordinates(radius,nox,rwater,cvec,cveci,cvec0,cvec0i,ipdc,rox,ncell)

!     extract the graph from the coordinates
      call get_hbondlist(nox,rwater,cvec,cveci,ipdc,neighbors,num_neighbors)

!     find the rings in the graph
      call findrings(nox,neighbors,num_neighbors,ringcount,ringlist,maxring)

!     clean up results - remove rings which are formed from pairs of smaller rings.
      call cleanup(ringcount,ringlist,neighbors,num_neighbors,maxring)

!     format the rings, and print the results
      call format_rings(nox,ncell,maxring,rox,ringcount,ringlist,cvec,cveci)

      end

      subroutine readcoordinates(radius,nox,rwater,cvec,cveci,cvec0,cvec0i,ipdc,rox,ncell)
      implicit none
      real(8), dimension(3,3) :: cvec,cveci,cvec0,cvec0i
      real(8), dimension(80000,3,3) :: rwater
      integer, dimension(80000,2) :: water_index
      real(8), dimension(80000,3) :: rox,rhy,rr
      real(8), dimension(80000,3) :: coord
      real(8) :: xdif,ydif,zdif,rr2,uu
      character(len = 2048) :: cvec_buffer
      character(len = 6), dimension(80000) :: atomicname
      integer ipdc,ncell
      integer i,j,k
      integer natoms,ios,nox,nhy
      real(8) :: rcut,rcut2
      real(8) :: x,y,z,radius
      real(8) :: xx,yy,zz
      real(8) :: astar_mag,bstar_mag,cstar_mag
      integer jcell1max,jcell2max,jcell3max
      integer jcell1,jcell2,jcell3
      character(len = 1) :: atname
      logical okflag

!     reads in the coordinates from the input file

      ios = 0 
      
      rcut2 = 1.2d0**2 

!     read in the coordinates and pack O atom coordinates into rox vector.

      read(10,*,iostat = ios) natoms
      read(10,'(A)',iostat = ios) cvec_buffer
      call get_cvec(cvec_buffer,cvec,cveci,uu)
      call mat3inverse(cvec,cveci,okflag)

      cvec0 = cvec
      cvec0i = cveci

      if(ios .ne. 0) then 
         write(*,*) 'ERROR IN READING COORDINATE FILE'
         stop
      endif

      do i = 1,natoms
         read(10,*,iostat = ios) atomicname(i),coord(i,1),coord(i,2),coord(i,3)
      end do 

!     find the size of the supercell contaning the cut-off

      call mat3inverse(cvec,cveci,okflag)

!     magnitude of the reciprocal lattice vectors
      astar_mag = dsqrt(cveci(1,1)**2 + cveci(1,2)**2 + cveci(1,3)**2)
      bstar_mag = dsqrt(cveci(2,1)**2 + cveci(2,2)**2 + cveci(2,3)**2)
      cstar_mag = dsqrt(cveci(3,1)**2 + cveci(3,2)**2 + cveci(3,3)**2)
      
      jcell1max = int(2.0d0 * radius * astar_mag) + 1
      jcell2max = int(2.0d0 * radius * bstar_mag) + 1
      jcell3max = int(2.0d0 * radius * cstar_mag) + 1

      write(*,*) 'convergence radius = ',radius
      write(*,*) 'supercell indices',jcell1max,jcell2max,jcell3max

      ncell = jcell1max * jcell2max * jcell3max

!     now replicate the coordinates

      nox = 0 
      nhy = 0 

      do jcell1 = 0,jcell1max - 1
         do jcell2 = 0,jcell2max - 1
            do jcell3 = 0,jcell3max - 1

               do i = 1,natoms
                  x = coord(i,1)
                  y = coord(i,2)
                  z = coord(i,3)

                  xx =  x + cvec(1,1) * jcell1 + cvec(1,2) * jcell2 + cvec(1,3) * jcell3
                  yy =  y + cvec(2,1) * jcell1 + cvec(2,2) * jcell2 + cvec(2,3) * jcell3
                  zz =  z + cvec(3,1) * jcell1 + cvec(3,2) * jcell2 + cvec(3,3) * jcell3

                  if(atomicname(i) .eq. 'O') then 
                     nox = nox + 1
                     rox(nox,1) = xx
                     rox(nox,2) = yy
                     rox(nox,3) = zz
                     write(17,*) "O ",rox(nox,1),rox(nox,2),rox(nox,3)
                  else if(atomicname(i).eq.'H') then 
                     nhy = nhy + 1
                     rhy(nhy,1) = xx
                     rhy(nhy,2) = yy 
                     rhy(nhy,3) = zz
                     write(17,*) "H ",rhy(nhy,1),rhy(nhy,2),rhy(nhy,3)
                  endif
                     
               end do 
            end do 
         end do 
      end do 

      cvec(1,1) = cvec(1,1) * jcell1max
      cvec(2,1) = cvec(2,1) * jcell1max
      cvec(3,1) = cvec(3,1) * jcell1max

      cvec(1,2) = cvec(1,2) * jcell2max
      cvec(2,2) = cvec(2,2) * jcell2max
      cvec(3,2) = cvec(3,2) * jcell2max

      cvec(1,3) = cvec(1,3) * jcell3max
      cvec(2,3) = cvec(2,3) * jcell3max
      cvec(3,3) = cvec(3,3) * jcell3max

      call mat3inverse(cvec,cveci,okflag)


!     now find all the waters molecule indices
!     NOTE, this only needs to be on the first step, assuming the bonding doesn't change.

      do i = 1,nox
         k = 0 
         do j = 1,nhy
            xdif = rox(i,1) - rhy(j,1)
            ydif = rox(i,2) - rhy(j,2)
            zdif = rox(i,3) - rhy(j,3)

            if(ipdc.eq.1) then 
               call nearest(xdif,ydif,zdif,cvec,cveci)
            endif

            rr2 = xdif**2 + ydif**2 + zdif**2

            if(rr2 < rcut2) then 
               k = k + 1
               water_index(i,k) = j
            endif

         end do 
      end do 

      do i = 1,nox

         rwater(i,1,1) = rox(i,1)
         rwater(i,1,2) = rox(i,2)
         rwater(i,1,3) = rox(i,3)

         j = water_index(i,1)

         rwater(i,2,1) = rhy(j,1)
         rwater(i,2,2) = rhy(j,2)
         rwater(i,2,3) = rhy(j,3)

         j = water_index(i,2)

         rwater(i,3,1) = rhy(j,1)
         rwater(i,3,2) = rhy(j,2)
         rwater(i,3,3) = rhy(j,3)

      end do 


      end subroutine readcoordinates


      subroutine get_hbondlist(nox,rwater,cvec,cveci,ipdc,neighbors,num_neighbors)
      implicit none
      integer nox,i,j,ipdc
      integer, dimension(100) :: ringcount
      real(8), dimension(80000,3,3) :: rwater
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      real(8) :: xdif,ydif,zdif,rcut2,x,y,z,rr2
      real(8), dimension(3,3) :: cvec,cveci
      real(8) :: x1dif,y1dif,z1dif,x2dif,y2dif,z2dif
      real(8) :: r1dis,r2dis,costheta,ff
      real(8) :: rohdis, roodis,theta,dot
      real(8) :: theta1,theta2,theta3,theta4
      real(8) :: angfac,pi
      logical :: hbonded

!     extracts the H-bond graph from the coordinates

      rcut2 = 3.5d0**2
      neighbors = 0 
      num_neighbors = 0 

      ringcount = 0 

      pi = 4.0 * atan(1.0d0) 

      angfac = 180.0d0 / pi 


!     build up an array of all the nearest neighbors of each O site.
      do i = 1,nox
         do j = i+1,nox

            xdif = rwater(i,1,1) - rwater(j,1,1)
            ydif = rwater(i,1,2) - rwater(j,1,2)
            zdif = rwater(i,1,3) - rwater(j,1,3)

            if(ipdc.eq.1) then 
               call nearest(xdif,ydif,zdif,cvec,cveci)
            endif
            rr2 = xdif**2 + ydif**2 + zdif**2
               
            if(rr2 < rcut2) then 
!     nearest neighbor found

!     now check the OHO angle, and only accept it as a H-bond if r_OH ^ r_OO < 30 degrees.

               roodis = dsqrt(rr2)
               hbonded = .false.

               x1dif = rwater(i,2,1) - rwater(i,1,1)
               y1dif = rwater(i,2,2) - rwater(i,1,2)
               z1dif = rwater(i,2,3) - rwater(i,1,3)
               
               if(ipdc.eq.1) then 
                  call nearest(x1dif,y1dif,z1dif,cvec,cveci)
               endif

               rohdis = dsqrt(x1dif**2 + y1dif**2 + z1dif**2)

               dot = xdif * x1dif + ydif * y1dif + zdif * z1dif
               dot = dot / (roodis * rohdis)
               theta1 = 180.0d0 - acos(dot) * angfac

               if(abs(theta1) < 30.0d0) then 
                  
                  hbonded = .true.
               else
                  
                  
                  x2dif = rwater(i,3,1) - rwater(i,1,1)
                  y2dif = rwater(i,3,2) - rwater(i,1,2)
                  z2dif = rwater(i,3,3) - rwater(i,1,3)
                  
                  if(ipdc.eq.1) then 
                     call nearest(x2dif,y2dif,z2dif,cvec,cveci)
                  endif

                  rohdis = dsqrt(x2dif**2 + y2dif**2 + z2dif**2)
                  
                  dot = xdif * x2dif + ydif * y2dif + zdif * z2dif
                  dot = dot / (roodis * rohdis)
                  theta2 = 180.0d0 - acos(dot) * angfac
                  
               if(abs(theta2) < 30.0d0) then 
                  
                  hbonded = .true.
               else
                  
                  x1dif = rwater(j,2,1) - rwater(j,1,1)
                  y1dif = rwater(j,2,2) - rwater(j,1,2)
                  z1dif = rwater(j,2,3) - rwater(j,1,3)
               
                  if(ipdc.eq.1) then 
                     call nearest(x1dif,y1dif,z1dif,cvec,cveci)
                  endif

                  rohdis = dsqrt(x1dif**2 + y1dif**2 + z1dif**2)
                  
                  dot = xdif * x1dif + ydif * y1dif + zdif * z1dif
                  dot = -dot / (roodis * rohdis)
                  theta3 = 180.0d0 - acos(dot) * angfac

                  
                  if(abs(theta3) < 30.0d0) then 
                     
                     hbonded = .true.
                  else
                     
                     x2dif = rwater(j,3,1) - rwater(j,1,1)
                     y2dif = rwater(j,3,2) - rwater(j,1,2)
                     z2dif = rwater(j,3,3) - rwater(j,1,3)
                     
                     if(ipdc.eq.1) then 
                        call nearest(x2dif,y2dif,z2dif,cvec,cveci)
                     endif

                     rohdis = dsqrt(x2dif**2 + y2dif**2 + z2dif**2)
                     
                     dot = xdif * x2dif + ydif * y2dif + zdif * z2dif
                     dot = -dot / (roodis * rohdis)
                     theta4 = 180.0d0 -  acos(dot) * angfac
                     
                     
                     if(abs(theta4) < 30.0d0) then
                        hbonded = .true.
                     endif
                  endif
               endif
            endif
            
            if(hbonded) then
               num_neighbors(i) = num_neighbors(i) + 1
               num_neighbors(j) = num_neighbors(j) + 1
               neighbors(i,num_neighbors(i)) = j 
               neighbors(j,num_neighbors(j)) = i
            endif
            
         endif
            
      end do


      end do 


      end subroutine get_hbondlist


      subroutine findrings(nox,neighbors,num_neighbors,ringcount,ringlist,maxring)
      implicit none
      integer i,j,k,l,m1,n1,m2,n2,m3,n3,m4,n4,m5,n5,m6
      integer i1,i2,i3,i4,i5,i6,i7,natoms,nox,count
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      integer, dimension(100) :: ringcount
      integer, dimension(10,80000,0:10) :: ringlist
      integer depth,ii,iox,npath,maxring
      integer, dimension(1000) :: pathlist

!     finds all the rings in the graph

      do i = 1,maxring
         ringcount(i) = 0 
      end do 

      do iox = 1,nox
         depth = 0 
         pathlist = 0 
         npath = 0 
         call tree(num_neighbors,neighbors,ringcount,nox,depth,iox,pathlist,npath,ringlist,maxring)
      end do 

      end subroutine findrings


      recursive subroutine tree(num_neighbors,neighbors,ringcount,nox,depth,iox,pathlist,npath,ringlist,maxring)
      implicit none
      integer, dimension(80000) :: num_neighbors
      integer, dimension(80000,10) :: neighbors
      integer, dimension(100) :: ringcount
      integer nox,depth,iox,iox2,neighbor,noxlist,npath,dd,count
      integer, dimension(1000) :: oxlist,pathlist
      integer, dimension(10,80000,0:10) :: ringlist
      integer :: maxring,sum
      logical visited
      integer j,k,npath2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!  RECURSIVE SUBROUTINE TO FIND ALL RINGS IN GRAPH
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      if(depth .ge. maxring) return

      pathlist(depth+1) = iox

      do neighbor = 1, num_neighbors(iox)
         iox2 = neighbors(iox,neighbor)

!     check if a ring has been found
         if(iox2 .eq. pathlist(1).and.depth.gt.1) then 

!     the following line ensures that not both directions around the ring are counted
            if(pathlist(depth+1).gt.pathlist(2)) then 

               ringcount(depth + 1) = ringcount(depth + 1) + 1
               count = ringcount(depth + 1)
               sum = 0 
               do k = 1,depth + 1
                  ringlist(depth+1,count,k) = pathlist(k)
                  sum = sum + pathlist(k)
               end do 
               ringlist(depth+1,count,0) = sum
            endif
         endif

!     check if this oxygen is on the list
         visited = .false.
         do j = 1,depth+1
            if(pathlist(j).eq.iox2) then 
               visited = .true.
               exit
            endif
         end do 

!     the iox2 > pathlist(1) ensures that each ring is only counted once

         if(.not.visited.and.iox2 > pathlist(1)) then
            dd = depth + 1
            call tree(num_neighbors,neighbors,ringcount,nox,dd,iox2,pathlist,npath2,ringlist,maxring)
        endif

      end do 

      return

      end subroutine tree



      subroutine cleanup(ringcount,ringlist,neighbors,num_neighbors,maxring)
      implicit none
      integer, dimension(100) :: ringcount,ringcount2
      integer, dimension(10,80000,0:10) :: ringlist,ringlist2
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      integer, dimension(100) :: ring1,ring2
      integer, dimension(1000) :: pathlist
      integer, dimension(100,100) :: shortestpath,shortestpath2
      integer :: ringsize
      integer, dimension(10,80000) :: exists
      integer min,kmin
      integer maxring
      logical shortcutfound,match
      integer maxdepth,l,n1,nn,k,k1,iox,ii,j,i,i1,diff

!     removes all the rings in the graph which can be formed from smaller rings
      
      exists = 1

      do i1 = 1,maxring
         do n1 = 1,ringcount(i1)

               ring1 = 0 
               do k1 = 1,i1
                  ring1(k1) = ringlist(i1,n1,k1)
               end do 

               do k = 1,i1
                  do l = 1,i1
                     shortestpath(k,l) = 1000
                  end do 
               end do 

               do k1 = 1,i1
                  pathlist = 0 

!     a maxdepth of only half the ringsize seems to work.
                  maxdepth = i1/2
                  ringsize = i1
                  iox = ring1(k1)
                  nn = k1
                  call tree2(maxdepth,ring1,ringsize,iox,nn,neighbors,num_neighbors,0,pathlist,shortestpath)
               end do

               shortcutfound = .false.
               do k = 1,i1
                  do l = 1,i1
                     diff = k - l
                     ii = int(abs(diff - i1 * anint(diff / real(i1))))
                     if(ii.ne.shortestpath(k,l)) shortcutfound = .true.
                     shortestpath2(k,l) = ii
                  end do 
               end do 

               if(shortcutfound) then 
                  exists(i1,n1) = 0
               endif

         end do 
      end do 


      ringlist2 = 0 

      do i = 1,maxring
         nn = 0 
         do j = 1,ringcount(i)
            if(exists(i,j).eq.1) then 
               nn = nn + 1
               do k = 0,i
                  ringlist2(i,nn,k) = ringlist(i,j,k)
               end do 
            endif
         end do 
         ringcount2(i) = nn
      end do 
      
      ringlist = ringlist2
      ringcount = ringcount2

      end subroutine cleanup

      recursive subroutine tree2(maxdepth,ring,ringsize,iox,nn,neighbors,num_neighbors,depth,pathlist,shortestpath)
      integer, dimension(100) :: ring
      integer, dimension(80000,10) :: neighbors
      integer, dimension(80000) :: num_neighbors
      integer, dimension(1000) :: pathlist
      integer depth
      integer, dimension(100,100) :: shortestpath
      logical visited
      integer :: ringsize

!     recursive subroutine for finding the shortest path between two nodes on ring

      do k = 1,ringsize
         if(iox.eq.ring(k)) then 
            if(depth.lt.shortestpath(k,nn)) then 
               shortestpath(k,nn) = depth
               shortestpath(nn,k) = depth
            endif
         endif
      end do 


      if(depth.ge.maxdepth) return

      pathlist(depth+1) = iox

      do neighbor = 1,num_neighbors(iox)
         iox2 = neighbors(iox,neighbor)

!     check if this oxygen is on the list
         visited = .false.
         do j = 1,depth+1
            if(pathlist(j).eq.iox2) then 
               visited = .true.
               exit
            endif
         end do 

         if(.not.visited) then
            newdepth = depth + 1
            call tree2(maxdepth,ring,ringsize,iox2,nn,neighbors,num_neighbors,newdepth,pathlist,shortestpath)
         endif

      end do 

      end subroutine tree2


      subroutine get_cvec(buffer,cvec,cveci,uu)
      implicit none
      character(len = 2048) :: buffer
      real(8), dimension(6) :: vec
      real(8) :: amag,bmag,cmag,alpha,beta,gamma,uu
      real(8), dimension(3,3) :: cvec,cveci
      character(len = 8) :: type
      logical okflag

!     returns the cell vectors in cartesian form from the cell-vector line in the inpu

      if(buffer.ne.'    ') then 
         read(buffer,*) type,vec
         if(type.eq.'XYZ') then
            cvec(1,1) = vec(1)
            cvec(1,2) = vec(2)
            cvec(2,2) = vec(3)
            cvec(1,3) = vec(4)
            cvec(2,3) = vec(5)
            cvec(3,3) = vec(6)
         else if(type.eq.'ANGLE') then 
!     if the cell vectors are in angle format then need to convert to cartesians

            amag = vec(1)
            bmag = vec(2)
            cmag = vec(3)
            alpha = vec(4)
            beta = vec(5)
            gamma = vec(6)

!     convert to cartesians
            call get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec)
         else
            write(*,*) 'ERROR: MUST SPECIFY EITHER XYZ&
     & OR ANGLE CVEC TYPE IN INPUT_FILE'
            stop
         endif

         call mat3inverse(cvec,cveci,okflag)
      else 
      endif

      end subroutine get_cvec


      subroutine get_cartesian_cvec(amag,bmag,cmag,alpha,beta,gamma,cvec0)
      implicit none
      real(8) :: alpha,beta,gamma,amag,bmag,cmag
      real(8) :: alpha_rad,beta_rad,gamma_rad
      real(8), dimension(6) :: param
      real(8) :: angfac,pi
      real(8), dimension(3,3) :: cvec0
!     converts cell-vectors from angular format to cartesians

      pi = 4.0d0*datan(1.d0) 
      angfac = 180.0d0 / pi

      alpha_rad = alpha / angfac
      beta_rad = beta / angfac
      gamma_rad = gamma / angfac

      cvec0 = 0.0d0

      cvec0(1,1) = amag
      cvec0(1,2) = bmag * dcos(gamma_rad)
      cvec0(2,2) = dsqrt(bmag**2 - cvec0(1,2)**2)
      cvec0(1,3) = cmag * dcos(beta_rad)
      cvec0(2,3) = (bmag * cmag * dcos(alpha_rad) - cvec0(1,2) * cvec0(1,3))/cvec0(2,2)
      cvec0(3,3) = dsqrt(cmag**2 - cvec0(1,3)**2 - cvec0(2,3)**2)

      end subroutine get_cartesian_cvec

      SUBROUTINE mat3inverse (A, AINV, OK_FLAG)

!***********************************************************************************************************************************
!  M33INV  -  Compute the inverse of a 3x3 matrix.
!
!  A       = input 3x3 matrix to be inverted
!  AINV    = output 3x3 inverse of matrix A
!  OK_FLAG = (output) .TRUE. if the input matrix could be inverted, and .FALSE. if the input matrix is singular.
!***********************************************************************************************************************************


      IMPLICIT NONE

      DOUBLE PRECISION, DIMENSION(3,3), INTENT(IN)  :: A
      DOUBLE PRECISION, DIMENSION(3,3), INTENT(OUT) :: AINV
      LOGICAL, INTENT(OUT) :: OK_FLAG

      DOUBLE PRECISION, PARAMETER :: EPS = 1.0D-10
      DOUBLE PRECISION :: DET
      DOUBLE PRECISION, DIMENSION(3,3) :: COFACTOR


      DET =   A(1,1)*A(2,2)*A(3,3)  &
            - A(1,1)*A(2,3)*A(3,2)  &
            - A(1,2)*A(2,1)*A(3,3)  &
            + A(1,2)*A(2,3)*A(3,1)  &
            + A(1,3)*A(2,1)*A(3,2)  &
            - A(1,3)*A(2,2)*A(3,1)

      IF (ABS(DET) .LE. EPS) THEN
         AINV = 0.0D0
         OK_FLAG = .FALSE.
         RETURN
      END IF

      COFACTOR(1,1) = +(A(2,2)*A(3,3)-A(2,3)*A(3,2))
      COFACTOR(1,2) = -(A(2,1)*A(3,3)-A(2,3)*A(3,1))
      COFACTOR(1,3) = +(A(2,1)*A(3,2)-A(2,2)*A(3,1))
      COFACTOR(2,1) = -(A(1,2)*A(3,3)-A(1,3)*A(3,2))
      COFACTOR(2,2) = +(A(1,1)*A(3,3)-A(1,3)*A(3,1))
      COFACTOR(2,3) = -(A(1,1)*A(3,2)-A(1,2)*A(3,1))
      COFACTOR(3,1) = +(A(1,2)*A(2,3)-A(1,3)*A(2,2))
      COFACTOR(3,2) = -(A(1,1)*A(2,3)-A(1,3)*A(2,1))
      COFACTOR(3,3) = +(A(1,1)*A(2,2)-A(1,2)*A(2,1))

      AINV = TRANSPOSE(COFACTOR) / DET

      OK_FLAG = .TRUE.

      RETURN

      END SUBROUTINE mat3inverse

      subroutine nearest(xdif,ydif,zdif,cmat,cmati)
      implicit none
      real(8),dimension(3,3) :: cmat,cmati
      real(8) :: xdif,ydif,zdif
      real(8),dimension(3) :: vec1,vec2

!     applies the minimum image convention for triclinic cells
!     note, assumes that both cmat and cmati are upper triangular

      vec1(1) = anint(cmati(1,1) * xdif + cmati(1,2) * ydif + cmati(1,3) * zdif)
      vec1(2) = anint(cmati(2,2) * ydif + cmati(2,3) * zdif)
      vec1(3) = anint(cmati(3,3) * zdif)

      xdif = xdif - (cmat(1,1) * vec1(1) + cmat(1,2) * vec1(2) + cmat(1,3) * vec1(3))
      ydif = ydif - (cmat(2,2) * vec1(2) + cmat(2,3) * vec1(3))
      zdif = zdif - (cmat(3,3) * vec1(3))

      end subroutine nearest

      subroutine format_rings(nox,ncell,maxring,rox,ringcount,ringlist,cvec,cveci)
      implicit none
      integer nox,ncell,maxring
      integer m,min,n,nplus,nminus,mm,nox_unit
      real(8) :: xdif,ydif,zdif
      integer, dimension(10) :: rlist,rlist2
      real(8), dimension(100,3) :: ringxyz,ring2xyz,ring3xyz
      real(8), dimension(80000,3) :: rox,rhy
      integer kmin,i,j,k,isign
      integer, dimension(100) :: ringcount
      integer, dimension(10,80000,0:10) :: ringlist
      real(8), dimension(3,3) :: cvec,cveci,cvec0,cvec0i

!     formats the rings such that the indices corresponds to those in the unit cell
!     and they're ordered with the minimum index first. 
!     Also print the results.


!     number of oxygens in unit cell
      nox_unit = nox / ncell

!     print results
      write(*,*) 
      write(*,*) 'number of rings per unit cell'
      write(*,*) '          ring size   #rings'
      do i = 3,maxring
         if(ringcount(i).ne.0) write(*,*) i,ringcount(i)/ncell
         do j = 1,ringcount(i)
            do k = 1,i
               mm = ringlist(i,j,k)
               ringxyz(k,1) = rox(mm,1)
               ringxyz(k,2) = rox(mm,2)
               ringxyz(k,3) = rox(mm,3)
            end do 

!     shift each ring so that the minimum index comes first.

            min = 2000
            do k = 1,i
               m = ringlist(i,j,k)
!     find the index wrt the unit cell, not the supercell
               n = mod(m-1,nox_unit)+1
               rlist(k) = n
               if(n.lt.min) then 
                  min = n
                  kmin = k
               endif
            end do 

!     decide which direction to go in. Choose direction such that next
!     index in ring is smaller than the index before.
            isign = 1
            nplus = mod(kmin + 1 -1,i) + 1
            nminus = mod(kmin -1 -1,i) + 1
            if(nminus.le.0) nminus = nminus + i
            
            if(rlist(nminus).lt.rlist(nplus)) isign = -1

!     get the ring in the new ordering.

            do k = 0,i-1
               if(isign.eq.1)  n = mod(kmin + k -1,i)+1
               if(isign.eq.-1)  n = mod(kmin - k -1,i)+1
               if(n.le.0) n = n + i
               rlist2(k+1) = rlist(n)
               ring2xyz(k+1,1) = ringxyz(n,1)
               ring2xyz(k+1,2) = ringxyz(n,2)
               ring2xyz(k+1,3) = ringxyz(n,3)
            end do 

!     now we're going to print out the ring such that the first xyz is equal 
!     to its position in the unit cell, and the subsequent coordinates describe 
!     the geometry of the ring without regard to periodic boundary conditions.

            xdif = ring2xyz(1,1)
            ydif = ring2xyz(1,2)
            zdif = ring2xyz(1,3)

            call nearest(xdif,ydif,zdif,cvec0,cvec0i)
            
            ring3xyz(1,1) = xdif
            ring3xyz(1,2) = ydif
            ring3xyz(1,3) = zdif

            do k = 2,i
               xdif = ring2xyz(k,1) - ring2xyz(k-1,1)
               ydif = ring2xyz(k,2) - ring2xyz(k-1,2)
               zdif = ring2xyz(k,3) - ring2xyz(k-1,3)
               
               call nearest(xdif,ydif,zdif,cvec,cveci)
               
               ring3xyz(k,1) = ring3xyz(k-1,1) + xdif
               ring3xyz(k,2) = ring3xyz(k-1,2) + ydif
               ring3xyz(k,3) = ring3xyz(k-1,3) + zdif
               
            end  do 

!     print out the ring

            write(20+i,16,advance = 'no') (rlist2(k),k=1,i)
            write(20+i,17) (ring3xyz(k,1),ring3xyz(k,2),ring3xyz(k,3),k=1,i)
 16         format(10i6)
 17         format(30f10.4)
            
         end do 
!     write a single blank line if zero ringcount, to make sure the file is blank.
         if(ringcount(i).eq.0) write(20+i,*)
      end do 

      end subroutine format_rings

