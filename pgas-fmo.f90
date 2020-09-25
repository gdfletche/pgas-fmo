!  Based on, 
!  [1]  D.G. Fedorov, K. Kitaura, J. Chem. Phys. 120, 6832-6840(2004) 
!  and,
!  [2]  T. Nakano, et al, Chemical Physics Letters 351 (2002) 475–480 
!  which describes the electrostatic approximations (long-range and 
!  intermediate-range - eqns. (12) and (13), respectively) 
!
!  Simplest possible implementation (0.0) assumes a single He atom 
!  per fragment 
!  Assumes the atomic basis consists of uncontracted s-type GTOs 
!  where the density is the product of the contraction weights 


 module pgasfmo_params  
 parameter  ( pi = 3.1415926535897931d0 )
 real(8)      sqrpi2
 parameter  ( sqrpi2 = ( pi**( -0.5d0 ) )*2.0d0 )
 real(8)      dtol,          rcut 
 parameter  ( dtol= 1.0d-10, rcut= 1.0d-12 ) 
 real(8)      tobohrs 
 parameter  ( tobohrs = 1.889725987722d0 )
 end module pgasfmo_params 



 module pgasfmo_info 
 implicit   none 
 integer(4) num_proc, myrank, errcon 
 integer(4) mpi_comm_compute  
 integer    num_comp 
 integer,   allocatable ::  cmap( : ) 
 end module pgasfmo_info 



 program   pgas_fmo  
 use       pgasfmo_params 
 use       pgasfmo_info  
 implicit  none 
 integer   nfrag,ngauss,myfrag,ifrag, i,j,k,l 
 real(8),  allocatable ::  xpnt( : )
 real(8),  allocatable ::  coef( : )
 real(8),  allocatable ::  geom( : , : )
 real(8),  allocatable ::  fock( : , : )
 real(8),  allocatable ::  dens( : , : )
 real(8),  allocatable ::  denovl( : )
 real(8)   rmedium,rlong, eri,rsep,mullpop,sum,ovl,ptch,potev  


 call pgasfmo_init 

 if ( myrank == 0 ) read *,  ngauss, nfrag      

 call pgasfmo_bcast( 'i', nfrag, 1  ) 
 call pgasfmo_bcast( 'i', ngauss, 1 ) 

 allocate( xpnt( ngauss ) )  
 allocate( coef( ngauss ) )  
 allocate( geom( 3, nfrag ) )  

 if ( myrank == 0 ) then 
!  read the GTO exponents and contraction coefficients 
 read *,  ( xpnt( i ), coef( i ), i = 1, ngauss ) 

!  read intermediate and long-range inter-fragment separation cutoffs 
!  and the (fake) Mulliken population charge 
 read *,  rmedium, rlong, mullpop

!  read the geometry 
 do  i  =  1,  nfrag 
 read *,  ( geom( j, i ), j = 1, 3 )  
 end do

 end if  !  master input  

!  broadcast input data to all compute ranks 
 call pgasfmo_bcast( 'd', rmedium, 1  ) 
 call pgasfmo_bcast( 'd',   rlong, 1  ) 
 call pgasfmo_bcast( 'd', mullpop, 1  ) 
 call pgasfmo_bcast( 'd', xpnt, ngauss  ) 
 call pgasfmo_bcast( 'd', coef, ngauss  ) 
 call pgasfmo_bcast( 'd', geom, 3*nfrag ) 

!  normalize the primitive GTO weights 
 do i = 1, ngauss 
 coef( i ) = coef( i )*( ( 2.0d0*xpnt( i ) )**0.75d0 )
 end do

!  scale the geometry to Bohrs for energy calculations in AU 
 do  i  =  1,  nfrag 
 geom( 1, i )  =  geom( 1, i )*tobohrs   
 geom( 2, i )  =  geom( 2, i )*tobohrs   
 geom( 3, i )  =  geom( 3, i )*tobohrs   
 end do
 rmedium = rmedium * tobohrs 
 rlong   = rlong   * tobohrs 

 allocate( fock( ngauss, ngauss ) )  
 do  i  =  1,  ngauss  
 do  j  =  1,  ngauss  
 fock( i, j )  =  0.0d0  
 end do
 end do

!  make the density - this could be input 
 allocate( dens( ngauss, ngauss ) )  
 do  i  =  1,  ngauss  
 do  j  =  1,  i  
 dens( i, j ) = coef( i ) * coef( j ) *2.0d0 
 dens( j, i ) = dens( i, j )  
 end do
 end do
 
!  make tr[ density * overlap ] for eq.(12) of [2]  
 allocate( denovl( ngauss ) )  
 do  i  =  1,  ngauss  
 sum = 0.0d0  
 do  j  =  1,  ngauss  
 sum = sum + dens( i, j )*ovl( j,i, xpnt, geom )  
 end do
 denovl( i ) = sum 
 end do
 
!  create the PGAS array to hold the fragment densities 
 call pgasfmo_create( ngauss**2, nfrag )
 
!  fill the PGAS array on each rank 
 call pgasfmo_put( myrank, ngauss**2, dens )   ! local PUT  

!  sync to ensure the array is filled before use 
 call pgasfmo_sync() 

 myfrag = myrank + 1 
!  Main loop over fragments, compute the e-e potential 
 do  ifrag  =  1,  nfrag   

!  compute the inter-fragment separation 
!  we would compute and store the fragment centers-of-mass, above 
!  but if the fragments consist only of single He atoms we can use their coords 
  rsep = sqrt( ( geom( 1, ifrag ) - geom( 1, myfrag ) )**2  &  
       +       ( geom( 2, ifrag ) - geom( 2, myfrag ) )**2  &  
       +       ( geom( 3, ifrag ) - geom( 3, myfrag ) )**2  )

!  begin compute kernel, consider offloading 
  do  i  =  1,  ngauss
  do  j  =  1,  ngauss 
!  options (approximations) for computing the potential 
  if ( rsep > rlong ) then    !  eq.13 of [2] 
   fock( i, j ) = fock( i, j ) + mullpop*ptch( myfrag,ifrag, i,j, xpnt, geom )   
  else  if ( rsep > rmedium ) then    !  eq.12 of [2] 
   do  k  =  1,  ngauss
   fock( i, j ) = fock( i, j ) + denovl( k )*eri( myfrag,ifrag, i,j,k,k, xpnt, geom )   
   end do  
  else                        !    eq.5 (RHS, 2nd term in curly brackets) of [1] 
   call pgasfmo_get( ifrag-1, ngauss**2, dens )  
   do  k  =  1,  ngauss
   do  l  =  1,  ngauss 
   fock( i, j ) = fock( i, j ) + dens( k, l )*eri( myfrag,ifrag, i,j,k,l, xpnt, geom )   
   end do  ;  end do  
  end if  
  end do  ;  end do     ! fock matrix elements
! end of compute kernel 

 end do  ! fragments (monomers)  


!  compute total interaction energy (eq.6 of [1]) - mainly for validation 
 call pgasfmo_get( myrank, ngauss**2, dens )  
 potev  =  0.0d0  
 do  i  =  1,  ngauss
 do  j  =  1,  ngauss 
 potev  =  potev  +  fock( i, j )*dens( i, j )  
 end do  ;  end do  
 call pgasfmo_gsumf( potev, 1 )  
 if ( myrank == 0 ) print *,'V= ', potev  


 deallocate( denovl )  
 deallocate( dens )  
 deallocate( fock )  
 deallocate( geom )  
 deallocate( coef )  
 deallocate( xpnt )  
 call pgasfmo_end()  
 end program pgas_fmo    




 real(8)    function eri( myfrag,ifrag, ib,jb,kb,lb, xpnt, geom )
 use        pgasfmo_params 
 implicit   none 
 integer    myfrag,ifrag, i,j,k,l 
 real(8)    xpnt(*), geom(3,*) 
 integer    ib,jb,kb,lb 
 real(8)    aij,dij,xij,yij,zij, akl,dkl, aijkl,tt,f0t  

 i = myfrag    !  formula is left 
 j = myfrag    !  general for later 
 k = ifrag     !  extension, see also 
 l = ifrag     !  ptch() and ovl() 
 eri = 0.0d0 
 aij = 1.0d0/( xpnt( ib ) + xpnt( jb ) ) 
 dij = exp( -xpnt( ib )*xpnt( jb )*aij*  &  
     ( ( geom( 1, i ) - geom( 1, j ) )**2   &
     + ( geom( 2, i ) - geom( 2, j ) )**2   &
     + ( geom( 3, i ) - geom( 3, j ) )**2  )  )*( aij**1.5d0 )  
 if ( abs( dij ) > dtol ) then      
 xij = aij*( xpnt( ib )*geom( 1, i ) + xpnt( jb )*geom( 1, j ) )  
 yij = aij*( xpnt( ib )*geom( 2, i ) + xpnt( jb )*geom( 2, j ) )  
 zij = aij*( xpnt( ib )*geom( 3, i ) + xpnt( jb )*geom( 3, j ) )  
 akl = 1.0d0/( xpnt( kb ) + xpnt( lb ) ) 
 dkl = dij*exp( -xpnt( kb )*xpnt( lb )*akl*  &  
     ( ( geom( 1, k ) - geom( 1, l ) )**2   &
     + ( geom( 2, k ) - geom( 2, l ) )**2   &
     + ( geom( 3, k ) - geom( 3, l ) )**2  )  )*( akl**1.5d0 )  
 if ( abs( dkl ) > dtol ) then      
 aijkl = ( xpnt( ib ) + xpnt( jb ) )*( xpnt( kb ) + xpnt( lb ) )  &  
       / ( xpnt( ib ) + xpnt( jb )  +  xpnt( kb ) + xpnt( lb ) )  
 tt = aijkl*( ( xij -akl*( xpnt( kb )*geom( 1, k ) + xpnt( lb )*geom( 1, l ) ) )**2  & 
            + ( yij -akl*( xpnt( kb )*geom( 2, k ) + xpnt( lb )*geom( 2, l ) ) )**2  & 
            + ( zij -akl*( xpnt( kb )*geom( 3, k ) + xpnt( lb )*geom( 3, l ) ) )**2  ) 
 f0t  =  sqrpi2 
 if ( tt > rcut )  f0t  =  ( tt**( -0.5d0 ) )*erf( sqrt(tt) ) 
 eri  =  dkl*f0t*sqrt(aijkl)  
 end if 
 end if  
 end function eri 




 real(8)    function ptch( myfrag,ifrag, ib,jb, xpnt, geom )
 use        pgasfmo_params 
 implicit   none 
 integer    myfrag,ifrag, i,j
 real(8)    xpnt(*), geom(3,*) 
 integer    ib,jb 
 real(8)    eij,aij,ovl,xij,yij,zij, tt,f0t  

 ptch  =  0.0d0 
 i = myfrag 
 j = myfrag 
 eij = xpnt( ib ) + xpnt( jb ) 
 aij = 1.0d0/eij 
 ovl = exp( -xpnt( ib )*xpnt( jb )*aij*  &  
     ( ( geom( 1, i ) - geom( 1, j ) )**2   &
     + ( geom( 2, i ) - geom( 2, j ) )**2   &
     + ( geom( 3, i ) - geom( 3, j ) )**2  )  )*( aij**1.5d0 )  
 if ( abs( ovl ) > dtol ) then      
 xij = aij*( xpnt( ib )*geom( 1, i ) + xpnt( jb )*geom( 1, j ) )  
 yij = aij*( xpnt( ib )*geom( 2, i ) + xpnt( jb )*geom( 2, j ) )  
 zij = aij*( xpnt( ib )*geom( 3, i ) + xpnt( jb )*geom( 3, j ) )  
 tt = (  ( xij -geom( 1, ifrag ) )**2  & 
      +  ( yij -geom( 2, ifrag ) )**2  & 
      +  ( zij -geom( 3, ifrag ) )**2  )*eij 
 f0t  =  sqrpi2 
 if ( tt > rcut )  f0t  =  ( tt**( -0.5d0 ) )*erf( sqrt(tt) ) 
 ptch  =  ovl*f0t*sqrt( eij )
 end if  
 end function ptch




 real(8)    function ovl( ib,jb, xpnt, geom )
 use        pgasfmo_params 
 implicit   none 
 integer    myfrag,ifrag, i,j
 real(8)    xpnt(*), geom(3,*) 
 integer    ib,jb 
 real(8)    aij

 i = myfrag 
 j = myfrag 
 aij = 1.0d0/( xpnt( ib ) + xpnt( jb ) ) 
 ovl = exp( -xpnt( ib )*xpnt( jb )*aij*  &  
     ( ( geom( 1, i ) - geom( 1, j ) )**2   &
     + ( geom( 2, i ) - geom( 2, j ) )**2   &
     + ( geom( 3, i ) - geom( 3, j ) )**2  )  )*( aij**1.5d0 )  
 end function ovl




 subroutine  pgasfmo_init
 use         pgasfmo_info 
 implicit    none 
 include      'mpif.h'
 integer     i  
 integer(4)  group_world, group_comp  
 integer(4), allocatable :: ranks( : ) 
 call mpi_init( errcon )
 call mpi_comm_rank( mpi_comm_world, myrank,   errcon )
 call mpi_comm_size( mpi_comm_world, num_proc, errcon ) 
 if ( mod( num_proc, 2 ) .ne. 0 ) stop '!odd rank count!' 
 num_comp = num_proc/2 
 if ( myrank == 0 ) print *,'total (compute) ranks= ',num_comp   
!     create compute process communicator to enable 
!     collective communications between them  
 allocate( ranks( 0 : num_comp - 1 ) ) 
 do i = 0, num_comp - 1  
 ranks(i) = i
 end do
 call mpi_comm_group( mpi_comm_world, group_world, errcon )
 call mpi_group_incl( group_world, num_comp ranks,   &  
                      group_comp, errcon )
 call mpi_comm_create( mpi_comm_world, group_comp,   &  
                       mpi_comm_compute, errcon )
 deallocate( ranks )
 if ( myrank >= num_comp ) then 
 call pgasfmo_server()   
 else
 allocate( cmap( 0 : num_comp ) ) 
 call pgasfmo_sync()   
 end if 
 end subroutine pgasfmo_init  




 subroutine    pgasfmo_bcast( chtype, buff, len )
 use           pgasfmo_info 
 implicit      none
 include      'mpif.h'
 character(1)  chtype
 integer       buff(*), len 
 integer(4)    type, from, len4   
 if ( chtype == 'i' ) then
 type = mpi_integer4
 else if ( chtype == 'd' ) then
 type = mpi_real8
 else
!       error
 end if
 len4 = len  
 from = 0
 call mpi_bcast( buff, len4, type, from, mpi_comm_compute, errcon )
 end subroutine  pgasfmo_bcast  




 subroutine  pgasfmo_gsumf( buff, datalen )
 use         pgasfmo_info 
 implicit    none
 include    'mpif.h'
 real(8)     buff(*)
 integer     datalen 
 integer(4)  len4, type, oper 
 integer     maxrcv
 parameter  (maxrcv = 16384 )
 real(8)     rcv( maxrcv )
 integer     i,j,loc,npass,len 
 type = mpi_real8
 oper = mpi_sum
 npass = 1 + ( datalen - 1)/maxrcv
 loc  = 1
 len4 = maxrcv
 do i = 1, npass
 if ( i == npass ) len4 = datalen - maxrcv*(npass-1)
 call mpi_allreduce( buff( loc ), rcv, len4, type,  &
                     oper, mpi_comm_compute, errcon  )
 do j = 1, len4 
 buff( loc + j - 1 ) = rcv( j )
 end do
 loc = loc + len4  
 end do
 end subroutine pgasfmo_gsumf  




 subroutine pgasfmo_sync() 
 use        pgasfmo_info 
 implicit   none 
 call mpi_barrier( mpi_comm_compute, errcon )
 end subroutine pgasfmo_sync 




 subroutine pgasfmo_create( nrows, ncols ) 
 use        pgasfmo_info 
 implicit   none 
 include    'mpif.h'
 integer    nrows, ncols
 integer    mincol,lftcol,icol,i 
 integer(4) message(3) 
 message(1) = 0          ! = do a CREATE  
 message(2) = nrows 
 message(3) = ncols 
 call pgasfmo_send( message, 12, myrank + num_comp )  
!  create the domain decomposition (map columns to ranks) 
 mincol = ncols/num_comp 
 lftcol = mod(ncols,num_comp) 
 icol = 1
 do i = 0, num_comp 
 cmap( i ) = icol
 icol = icol + mincol
 if ( i .lt. lftcol ) icol = icol + 1
 end do 
 call mpi_barrier( mpi_comm_world, errcon )
 end subroutine pgasfmo_create    




 subroutine pgasfmo_decomp( i, ncols, data_loc, item ) 
 use        pgasfmo_info 
 implicit   none
 integer    i, ncols, data_loc, item 

! return item 'item' on rank 'data_loc' of global addr 'i' 

 end subroutine pgasfmo_decomp




 subroutine pgasfmo_get( data_loc, len, buff )
 use        pgasfmo_info 
 implicit   none
 integer    data_loc, len, item 
 real(8)    buff(*)  
 integer(4) message(3) 
 message(1) = 1          ! = do a GET
 message(2) = len  
 message(3) = item 
 call pgasfmo_send( message, 12, data_loc + num_comp )  
 call pgasfmo_recv( buff, len*8, data_loc + num_comp )  
 end subroutine pgasfmo_get




 subroutine pgasfmo_put( data_loc, len, buff )
 use        pgasfmo_info 
 implicit   none
 integer    data_loc, len, item 
 real(8)    buff(*)  
 integer(4) message(3) 
 message(1) = 2          ! = do a PUT
 message(2) = len  
 message(3) = item 
 call pgasfmo_send( message, 12, data_loc + num_comp )
 call pgasfmo_send( buff, len*8, data_loc + num_comp )  
 end subroutine pgasfmo_put




 subroutine pgasfmo_send( sendbuff, len, to )
 implicit   none
 include    'mpif.h'
 integer    sendbuff(*)  
 integer    len, to 
 integer(4) len4, mpityp, to4, msgtag, errcon 
 msgtag = 270908 
 mpityp = mpi_byte
 len4   = len
 to4    = to
 call mpi_ssend( sendbuff, len4, mpityp, to4, msgtag,  & 
                 mpi_comm_world, errcon )
 end subroutine pgasfmo_send




 subroutine pgasfmo_recv( recvbuff, len, from )
 implicit   none
 include    'mpif.h'
 real(8)    recvbuff(*)  
 integer    len, from 
 integer(4) len4,mpityp,from4,msgtag,status(mpi_status_size),errcon 
 msgtag = 270908 
 mpityp = mpi_byte
 len4   = len
 from4  = from
 call mpi_recv( recvbuff, len4, mpityp, from4, msgtag,  &  
                mpi_comm_world, status, errcon )
 end subroutine pgasfmo_recv




 subroutine pgasfmo_rcvany( recvbuff, len, from )
 implicit   none
 include    'mpif.h'
 integer    recvbuff(*)  
 integer    len, from 
 integer(4) len4,mpityp,msgtag,status(mpi_status_size),errcon 
 msgtag = 270908 
 mpityp = mpi_byte
 len4   = len
 call mpi_recv( recvbuff, len4, mpityp, mpi_any_source, msgtag,  &  
                mpi_comm_world, status, errcon )
 from = status(mpi_source)
 end subroutine pgasfmo_rcvany




 subroutine pgasfmo_server
 implicit   none
 include    'mpif.h'
 logical    serving 
 integer(4) message(3), errcon 
 integer    from,task,len 
 real(8), allocatable :: pgasfmo_data( : )   ! FMO1 density 

 serving = .true.
 do while ( serving )
  call pgasfmo_rcvany( message, 12, from )
  task = message(1)
   len = message(2)  
   len = message(3)  
       if ( task == 0 ) THEN
   allocate( pgasfmo_data( len ) )  !  create local part of the PGAS 
   call mpi_barrier( mpi_comm_world, errcon )
  else if ( task == 1 ) then        !  GET  
   call pgasfmo_send( pgasfmo_data, len*8, from )
  else if ( task == 2 ) then        !  PUT  
   call pgasfmo_recv( pgasfmo_data, len*8, from )
  else if ( task == 3 ) then
   serving = .false.                !  done serving 
  end if 
 end do 
 call mpi_barrier( mpi_comm_world, errcon )
 deallocate( pgasfmo_data )  
 call mpi_finalize( errcon )
 stop   ! server does not return to the main program 
 end subroutine pgasfmo_server




 subroutine pgasfmo_end      !  called by compute ranks 
 use        pgasfmo_info 
 integer(4) message(3) 
 message(1) = 3      ! = terminate server 
 call pgasfmo_send( message, 12, myrank + num_comp )
 call mpi_barrier( mpi_comm_world, errcon )
 call mpi_finalize( errcon )
 end subroutine pgasfmo_end  




