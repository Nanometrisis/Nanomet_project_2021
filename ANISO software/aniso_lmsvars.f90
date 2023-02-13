module lmsvars
    implicit NONE
    integer,parameter:: dp= kind(1.d0)
    
    real(dp),allocatable:: lmsfx(:,:),lmsfy(:,:)
    integer,allocatable:: lmsn(:)
    character(len=30),allocatable:: lmslabel(:)
    type ::animat
     character(len=20):: label
     complex(dp)      :: eps
     complex(dp)      :: mu = cmplx(1.d0,0.d0,dp)
     integer          :: type ! 1 constant, 2: read file 
    end type animat
    
contains
   integer function lmsgetid(label)
! It finds the index where label is stored in arrays lmslabel     
   character(len=30):: label
   integer :: n,i
   n = size(lmslabel)
   lmsgetid = 0
   do i=1,n
      if (trim(label)==trim(lmslabel(i))) then 
         lmsgetid = i
         exit
      end if 
      if (i==n) then 
         write(6,*) 'Label:',label,' NOT FOUND  Stoping.'
         stop 
      end if 
   end do
                 
   end function lmsgetid 

   subroutine readinput(nlmax,nfreq,start_freq,end_freq,theta,nlayers,type, &
        &                 thickness,pitch,linewidth,layermat,layermatlw)
   implicit NONE
   integer:: nlmax
   integer,intent(out):: nfreq,nlayers
   real(dp), intent(out):: start_freq,end_freq,theta
   real(dp),intent(out):: thickness(nlmax)
   real(dp),intent(out):: pitch(nlmax),linewidth(nlmax)        ! pitch and linewidth of gratting
   character(len=20),intent(out):: layermat(nlmax),layermatlw(nlmax)
   integer,intent(out):: type(nlmax)

   integer:: i,itype,read_unit
   
   character(len=20):: LTYPE(3)
   character(len=200):: line
   integer:: ios,ios1
   real(dp):: pi
   pi=4.d0*atan(1.d0)
   ios = 0
   LTYPE(1)='HOMOGENEOUS LAYER  '
   LTYPE(3)='GRATTING    LAYER  '
    read_unit = 9
    open(read_unit,file='anisolayer.in',iostat=ios)
    if (ios/=0 ) stop 'Error opening: anisolayer.in'

   read(read_unit, '(A)', iostat=ios1) line
   read(line,*) nfreq,start_freq,end_freq
   write(6,"('N:',I6,'  START(eV): ',F12.5,'  END(eV): ',F12.5)")  &
       &         nfreq,start_freq,end_freq
   
   read(read_unit,'(A)') line
   read(line,*) theta
   write(6,"('Angle of Incidence (deg):',F12.5)") theta
   theta = theta*pi/180.d0 ! in rad
   
   read(read_unit,'(A)') line
   read(line,*) nlayers
   write(6,*) 'Number of layers :', nlayers
   do i=1,nlayers
      read(read_unit,'(A)') line
      read(line,*) itype
      write(6,*) 'Reading ',LTYPE(itype)
      if (itype==1) then 
   !  write(6,*) 'here 1'
         read(line,*) type(i),thickness(i),layermat(i)
         write(6,"(A20,' Thickness: ',F12.5,' Material: ',A)") &
             &   LTYPE(type(i)),thickness(i),trim(layermat(i))
 
      elseif (itype==3) then
          !   write(6,*) 'here 2' 
         read(line,*) type(i),thickness(i),pitch(i),linewidth(i),   &
         &            layermatlw(i),layermat(i)
 ! matlw = e1, mat= e2
 
         write(6,"(A20,' Thickness(um): ',F12.5,'Pitch(um): ',F12.5, &
                       'Linewidth(nm): ',F12.5,2A)") &
             &      LTYPE(type(i)),thickness(i),pitch(i),linewidth(i),&
             &      trim(layermatlw(i)),trim(layermat(i))
 
      else
         write(6,*) 'Error in reading input 1'
         stop
      end if 
   end do
   close(read_unit)
   write(6,*) ' Reading Control Input Finished'
end subroutine  
 
 
 
subroutine read_file(filename,label,icol)
 ! This is used to read data from file=filename
 ! The data columns icol from the file are assigned
 ! to arrays

!     This version just reads 4 columns of the file
!     Lines starting with !,#,@ are ignored        
implicit NONE

integer,parameter:: nmax=10000
integer,parameter:: colmax=20
character(len=30),intent(in):: filename
character(len=30),intent(in):: label    ! same as lmslabel
integer,intent(in)::icol
!
integer :: read_unit
character(len=30)::label1
real(kind=dp):: mydata(nmax,colmax)  ! local variable
integer:: ndata
          
integer:: ios,line,ios1,n,nc,i,newid
character(len=200) buffer 
       
read_unit = 22
if (trim( filename ).ne.' ' ) then ! read file

   write(6,*)' read_file sub, Opening : ',filename
   ios = 0     
   open(unit=read_unit, file=filename, iostat=ios)
             
   if ( ios /= 0 ) stop "Error opening file"
   line = 0
   ios1 = 0
   n = 0
   do
     read(read_unit, '(A)', iostat=ios1) buffer
     if (ios1 /= 0) exit
     if (buffer(1:1)=='!'.or.buffer(1:1)=='@'.or.buffer(1:1)=='#') cycle
     n = n + 1
     if (n>nmax) STOP 'file has many lines increase NMAX'
     nc = 4
     read(buffer,*) mydata(n,1:4)  ! This version reads 3 columns of the file
   end do
   ! Consistency test for splines 
     if (mydata(3,1)<mydata(2,1)) then 
         write(6,*) 'Error Reading Material File: ',filename
         write(6,*) 'x-range must be accending in this version'
         write(6,*) 'Spline interpolation will fail!'
         STOP 
     end if 
    
   !
   ndata = n
           
   label1 = " "  ! lmslabels are intialized with blanck for this to work

   newid  = lmsgetid(label1)
   if (newid>size(lmslabel)) then 
      write(6,*) 'Increase lmsdata dimensions to use more data'
      write(6,*) 'Error Stop in read_file'
      stop
   end if 
           
   lmslabel(newid)= label
   lmsn(newid) = ndata
   lmsfx(1:ndata,newid) = mydata(1:ndata,1)
   lmsfy(1:ndata,newid) = mydata(1:ndata,icol)
 
   write(6,123) trim(filename),n, nc
123    format('Reading file : ',(A),', contains ',I5,' lines and ',I5,' columns'  )
close(read_unit)
end if
       
end subroutine

subroutine mulinterp(xnew,ynew,nnew,xdata,ydata,ndata)
! Interopate using SPLINES driver
! xdata(ndata), ydata(ndata) are the original
! ynew(xnew) are the interpolation this version
! uses arrays for spead up
!   .           
! CAUTION: Data to interpolate should be accending!   
    implicit none
    integer:: nnew,ndata  
    real(dp) :: xnew(nnew),ynew(nnew)
    real(dp) :: xdata(ndata),ydata(ndata)
    real(dp),allocatable  :: b(:),c(:),d(:)
    integer:: i 

    allocate (b(ndata),c(ndata),d(ndata))
    b= 0.d0
    c= 0.d0
    d= 0.d0
    call spline (ndata,xdata,ydata,b,c,d)
    do i=1,nnew
       ynew(i) = seval(ndata,xnew(i),xdata,ydata,b,c,d)
    end do
    
    deallocate(b,c,d)
end subroutine mulinterp

    FUNCTION SEVAL (N,U,X,Y,B,C,D)
       !------------------------------------------------------------------------
       !     EVALUATE A CUBIC SPLINE INTERPOLATION OF A DISCRETE FUNCTION F(X),
       !     GIVEN IN N POINTS X(I), Y(I). THE B, C AND D COEFFICIENTS DEFINING
       !     THE BEST CUBIC SPLINE FOR THE GIVEN POINTS, ARE CALCULATED BEFORE
       !     BY THE SPLINE SUBROUTINE.
       !
       !     INPUTS:
      !           CAUTION X data SHOULD be INCREASING
       !     N       NUMBER OF POINTS OF CURVE Y = F(X)
       !     U       ABSCISSA OF POINT TO BE INTERPOLATED
       !     X,Y     TABLES OF DIMENSION N, STORING THE COORDINATES
       !             OF CURVE F(X)
       !     B,C,D   TABLES STORING THE COEFFICIENTS DEFINING THE
       !             CUBIC SPLINE
       !
       !     OUTPUTS:
       !     SEVAL   INTERPOLATED VALUE
       !             = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
       !             WITH DX = U-X(I), U BETWEEN X(I) AND X(I+1)
       !
       !     REFERENCE :
       !     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
       !     COMPUTATIONS. PRENTICE-HALL,INC.
       !------------------------------------------------------------------------
             REAL *8 B(N),C(N),D(N),X(N),Y(N),U,DX,SEVAL
             integer:: n
             integer:: i,j,k
             !DATA I/1/
             I = 1
       !     BINARY SEARCH
       
             IF (I.GE.N) I = 1
             IF (U.LT.X(I)) GO TO 10
             IF (U.LE.X(I+1)) GO TO 30
          10 I = 1
             J = N+1
          20 K = (I+J)/2
             IF (U.LT.X(K)) J = K
             IF (U.GE.X(K)) I = K
             IF (J.GT.I+1) GO TO 20
       
       !     SPLINE EVALUATION
       
          30 DX = U-X(I)
             SEVAL = Y(I)+DX*(B(I)+DX*(C(I)+DX*D(I)))
             RETURN
             END
       
             SUBROUTINE SPLINE (N,X,Y,B,C,D)
       !---------------------------------------------------------------------
       !     THIS SUBROUTINE CALCULATES THE COEFFICIENTS B,C,D OF A CUBIC
       !     SPLINE TO BEST APPROXIMATE A DISCRETE FUNCTION GIVEN BY N POINTS
       !
       !     INPUTS:
       !     N       NUMBER OF GIVEN POINTS
       !     X,Y     VECTORS OF DIMENSION N, STORING THE COORDINATES
       !             OF FUNCTION F(X)
       !
       !     OUTPUTS:
       !     A,B,C   VECTORS OF DIMENSION N, STORING THE COEFFICIENTS
       !             OF THE CUBIC SPLINE
       !
       !     REFERENCE:
       !     FORSYTHE,G.E. (1977) COMPUTER METHODS FOR MATHEMATICAL
       !     COMPUTATIONS. PRENTICE-HALL,INC.
       !---------------------------------------------------------------------
             IMPLICIT REAL *8 (A-H,O-Z)
             integer :: n,i,nm1,l
             DIMENSION B(N),C(N),D(N),X(N),Y(N)
             NM1 = N-1
             IF (N.LT.2) RETURN
             IF (N.LT.3) GO TO 50
       
       !     BUILD THE TRIDIAGONAL SYSTEM
       !     B (DIAGONAL), D (UPPERDIAGONAL) , C (SECOND MEMBER)
       
             D(1) = X(2)-X(1)
             C(2) = (Y(2)-Y(1))/D(1)
             DO 10 I = 2,NM1
             D(I) = X(I+1)-X(I)
             B(I) = 2.D0*(D(I-1)+D(I))
             C(I+1) = (Y(I+1)-Y(I))/D(I)
             C(I) = C(I+1)-C(I)
          10 CONTINUE
       
       !     CONDITIONS AT LIMITS
       !     THIRD DERIVATIVES OBTAINED BY DIVIDED DIFFERENCES
       
             B(1) = -D(1)
             B(N) = -D(N-1)
             C(1) = 0.D0
             C(N) = 0.D0
             IF (N.EQ.3) GO TO 15
             C(1) = C(3)/(X(4)-X(2))-C(2)/(X(3)-X(1))
             C(N) = C(N-1)/(X(N)-X(N-2))-C(N-2)/(X(N-1)-X(N-3))
             C(1) = C(1)*D(1)*D(1)/(X(4)-X(1))
             C(N) = -C(N)*D(N-1)**2/(X(N)-X(N-3))
       
       !     FORWARD ELIMINATION
       
          15 DO 20 I = 2,N
             T = D(I-1)/B(I-1)
             B(I) = B(I)-T*D(I-1)
             C(I) = C(I)-T*C(I-1)
          20 CONTINUE
       
       !     BACK SUBSTITUTION
       
             C(N) = C(N)/B(N)
             DO 30 L = 1,NM1
             I = N-L
             C(I) = (C(I)-D(I)*C(I+1))/B(I)
          30 CONTINUE
       
       !     COEFFICIENTS OF 3RD DEGREE POLYNOMIAL
       
             B(N) = (Y(N)-Y(NM1))/D(NM1)+D(NM1)*(C(NM1)+2.D0*C(N))
             DO 40 I = 1,NM1
             B(I) = (Y(I+1)-Y(I))/D(I)-D(I)*(C(I+1)+2.D0*C(I))
             D(I) = (C(I+1)-C(I))/D(I)
             C(I) = 3.D0*C(I)
          40 CONTINUE
             C(N) = 3.D0*C(N)
             D(N) = D(NM1)
             RETURN
       
       !     CAS N = 2
       
          50 B(1) = (Y(2)-Y(1))/(X(2)-X(1))
             C(1) = 0.D0
             D(1) = 0.D0
             B(2) = B(1)
             C(2) = 0.D0
             D(2) = 0.D0
             RETURN
             END 
             
   subroutine readmaterials(nmatmax,mat,nmat)
! This is to read material parameters from a file
!   *SI_ASPNES    ,  aspnes_si.txt
!   *AIR          ,  C             , 1.d0    , 0.d0
!   *ARC          , arc_nk.txt
!   *RES          , resist_nk.txt
!   Typical file structure lines starting with * are only considered
!   a label eg 'SI_ASPNES', 'AIR' is used for id
!   the second column is a filename that has 4 columns for the material
!   ev, n, k, ignore
!   Read in in read_file and assigned to lms* arrays
!   If the 2nd column is C the 2 numbers are n, k of the material                                   
!   Overall:
!      If we have constant parameters values go to mat(imat)%eps
!      If the parameters are freq dependent values go to lmsx, lmsy, 
!      lmslabel etc data structure.     
!                                            3Jan2023  
    integer,intent(in):: nmatmax
    integer,intent(out)::nmat
    type(animat),intent(out):: mat(nmatmax)
!
    character(len=200):: buffer
    integer:: ios,ios1
    integer:: read_unit
    character(len=30) :: label1,type1,label 
    real(dp):: epsr,epsi,ni,kappa
!    
    read_unit = 9
    write(6,*) 'Reading Materials from input'
    write(6,*)
    ios = 0
    open (read_unit,file='anisolayer.in',iostat=ios)
    if (ios/=0) stop ' Error reading inisolayer.in'
    
    ios1 = 0
    nmat = 0
    do
        read(read_unit, '(A)', iostat=ios1) buffer
        !write(6,*)'-->', buffer 
        if (ios1 /= 0) exit
        if ( buffer(1:1) == '*') then 
            nmat = nmat + 1
            write(6,*) ' Material no: ', nmat
            if (nmat>nmatmax) stop 'nmat'
 
            read(buffer,*) label1,type1
    !        write(6,*) 'mat rere',label1,type1
            type1 = trim(adjustl(type1))
    !        write(6,*) 'mat rere',label1,type1
            if (type1=='C') then 
                read(buffer,*) label1,type1,ni,kappa
                mat(nmat)%label=trim(label1(2:))
                mat(nmat)%type=1
                epsr = ni*ni-kappa*kappa
                epsi = 2.d0*ni*kappa
                mat(nmat)%eps=cmplx(epsr,epsi,dp)
                write(6,*) nmat,mat(nmat)

            else 
                mat(nmat)%label=trim(label1(2:))
                mat(nmat)%type=2
                label = trim(label1(2:))//'_NI'
                write(6,*) label
                call read_file(type1,label,2)
                label = trim(label1(2:))//'_K'
                write(6,*) label
                call read_file(type1,label,3)
            end if 
                      
        end if  
                
    end do
  end subroutine readmaterials   

  subroutine getmaterial(nmatmax,label,mat,nmat,eps,zval)
!
! This is used to retrieve material parameters 
!         
    implicit none
    integer,intent(in):: nmatmax ! dimension 
    type(animat) :: mat(nmatmax) ! materials data 
    integer,intent(in):: nmat    ! number of materials read in
    real(dp),intent(in):: zval   ! 
    complex(dp),intent(out) ::eps
    character(*),intent(in)::label 
!    
    real(dp),allocatable :: xnew(:),ni(:),kappa(:)
    integer :: i,id,nn,nnew,i1
    character(len=40) lab
    
    do i=1,nmat
       if (trim(label)==trim(adjustl(mat(i)%label))) then

          if (mat(i)%type==1) then 
             eps = mat(i)%eps
          else
            ! Get ni at frequency = zval
            ! Files have eV
             lab = trim(label)//'_NI'
             id = lmsgetid(lab)
             nn = lmsn(id)        
             nnew = 1
             allocate (xnew(nnew),ni(nnew),kappa(nnew))
             xnew(1) = zval
             call mulinterp(xnew,ni,nnew,lmsfx(1:nn,id),lmsfy(1:nn,id),nn)
            ! Now get Kappa
            ! 
             lab = trim(label)//'_K'
             id = lmsgetid(lab)
             nn = lmsn(id)   
             call mulinterp(xnew,kappa,nnew,lmsfx(1:nn,id),lmsfy(1:nn,id),nn)
             eps  = cmplx(ni(1)**2-kappa(1)**2,    &
              &                 2.d0*ni(1)*kappa(1),dp)
             deallocate(xnew,ni,kappa)
          end if
       end IF
    end do
  end subroutine getmaterial   
  end module lmsvars
