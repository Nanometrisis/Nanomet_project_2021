program anisolayers_main
   !1. Effective-medium model for fast evaluation of scatterometric
   !measurements on gratings
   !A. Weidnera, M. Slodowskib, C. Halmb, C. Schneidera, and L. Pfitznera
   !Metrology, Inspection, and Process Control for Microlithography XVIII, edited by
   !Richard M. Silver, Proceedings of SPIE Vol. 5375, p 232
   ! see Refs also
   ! 2. Metal-coated magnetic nanoparticles in an optically active medium: A nonreciprocal metamaterial
   !   Aristi Christofi and Nikolaos Stefanou, PRB 97, 125129 (2018).
   !    
     use lmsvars
   use globalpar
    implicit none
  integer,parameter:: nlmax=500
  integer,parameter:: nmatmax=100
  !
  !
  !      1     2     3       4       layers
  !   |-----|-----|------|------|
  !   0     1     2      3      4    interfaces
  !
  !real(dp):: filedata(nmax,colmax,maxsplines)
  complex(dp):: eps(nlmax),mu(nlmax),epsanis(3,3,nlmax),muanis(3,3,nlmax)
  real(dp):: thickness(nlmax)
  real(dp):: pitch(nlmax),linewidth(nlmax)        ! pitch and linewidth of gratting
  type(animat):: mat(nmatmax) ! Maximum materials correct
  character(len=20)::layermat(nlmax),layermatlw(nlmax)
  integer:: type(nlmax)
  complex(dp) :: eps1g,eps2g        ! materials in gratting 
  logical:: lisotropic(nlmax)

  !integer:: itype(nlmax)   ! type of layer ilay = 1  slab 
!                                                 2  effective medium gratting      
!                                                 3  interface
!
!
  real(dp)::theta,fi
  integer:: nlay,ilay,ifreq,nfreq,nmat,id,n0,nfuns
  real(dp):: test,omega_start,omega_end,zval,div,lambda
  complex(dp):: eps_or,eps_eo
  complex(dp):: TRA(2,2),REF(2,2)
  real(dp),allocatable::xx(:),yy(:),b1(:),c1(:),d1(:)
  real(dp),allocatable::omega(:),ni(:),ka(:)
  real(dp):: start,end,x0,y0,scaleunits
  integer:: i,nn,n1

   n0=10000
   nfuns = 100
allocate(lmsfx(n0,nfuns),lmsfy(n0,nfuns),lmsn(nfuns),lmslabel(nfuns))
lmslabel(1:nfuns)=' '

  ! Read input for structures
  call readinput(nlmax,nfreq,omega_start,omega_end,theta,nlay,type, &
  &                 thickness,pitch,linewidth,layermat,layermatlw)
  ! Units
  SCALEUNITS=0.197327d0

   thickness = thickness/SCALEUNITS
   pitch = pitch 
   linewidth = linewidth 
write(6,*)'pitch', pitch(1:nlay)
  call readmaterials(nmatmax,mat,nmat)
  write(6,*)' I have read ',nmat,' materials from input'
  ! Geometry
  !nlay = 3  ! number of layers
  ! Initialization
  do ilay=1,nlay
     mu(ilay)  = 1.d0     ! mu of each layer
     epsanis(1:3,1:3,ilay) = 0.d0   ! tensor epsilon for each layer
     muanis(1:3,1:3,ilay)  = 0.d0   ! tensor mu for each layer
     muanis(1,1,ilay)=1.d0
     muanis(2,2,ilay)=1.d0
     muanis(3,3,ilay)=1.d0
  end do

  fi = 0 ! deg  

  call makefunlabels
  allocate (output_data(50000,nmaxfun))

  do ilay = 1,nlay
     lisotropic(ilay) = .true.
     test = sum(abs(epsanis(1:3,1:3,ilay)))
     if ( test> 1.d-10.or.pitch(ilay)>1.d-10 ) lisotropic(ilay) = .false. 
  end do
  if (.not.(lisotropic(1).and.lisotropic(nlay))) then
     write(6,*) 'First and last layers should be isotropic in this version'
     write(6,*) 'STOPING !'
     stop
  end if    
!  
  if (nfreq>1) then 
     div = (omega_end - omega_start)/(nfreq-1)
  else
     div = 0.0d0
  end if
  
  write(6,*) 'Freq(eV)   Angle    TR_S     TR_P    Ref_S     Ref_P  '


  

  do ifreq=1,nfreq ! loop in frequencies
     zval = omega_start + div*(ifreq-1)
!!!!!!!!!!!!!!!!!!!!!!!!!! Set Up Effective Medium Parameters (Wavelenght dependence !!!!!!!!!!!
     id=get_funid('OMEGA')
          if (id>0) output_data(ifreq,id) =zval
     do ilay=1,nlay
      ! Get optical parameters
      lambda = 1.24d0/zval
      call getmaterial(nmatmax,layermat(ilay),mat,nmat,eps(ilay),zval)
      write(6,"('EPS:',I3,3F12.5,' ',A20)") ilay,zval,eps(ilay) ,layermat(ilay)

      if (pitch(ilay)>1.d-8) then 
         call getmaterial(nmatmax,layermatlw(ilay),mat,nmat,eps1g,zval)
         call getmaterial(nmatmax,layermat(ilay),mat,nmat,eps2g,zval)
         write(6,*) 'eps',eps1g,eps2g
         !eps1g = real(eps1g)
         !eps2g = real(eps2g)

         call getOpticalParametrers(theta,pitch(ilay),linewidth(ilay),  &
         &    eps1g,eps2g,eps(1),lambda,eps_or,eps_eo )

         epsanis(1,1,ilay) = eps_eo
         epsanis(2,2,ilay) = eps_or
         epsanis(3,3,ilay) = eps_or 
         write(6,"('aniso',10E16.8)") lambda,eps_eo,eps_or
      else
 
      end if

     end do
!!!!!!!!!!!!!!!!!!!!!!!!!! Finished !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!
! Now run T, R for layers
!     
     !write(6,*)'before TR',eps(1:nlay)
     call TRanisoml(nlay,nlmax,eps,epsanis,mu,muanis,lisotropic, &
     &              thickness,theta,fi,zval ,ifreq)
     !write(16,"(10D16.8)") zval,theta,TRA(1,1),TRA(2,2),REF(1,1),REF(2,2)
     write(6,*) 'output_data',n_output_data
     write(21,"(100E16.8)") output_data(ifreq,1:n_output_data)
   end do    ! end of loop in frequencies


end program anisolayers_main

subroutine getOpticalParametrers(theta,pitch,linewidth,eps1,eps2,eps0,lambda,eps_or,eps_eo )
! This is using  
!   Effective-medium model for fast evaluation of scatterometric
!   measurements on gratings
!   A. Weidnera, M. Slodowskib, C. Halmb, C. Schneidera, and L. Pfitznera      
!   Metrology, Inspection, and Process Control for Microlithography XVIII, edited by
!   Richard M. Silver, Proceedings of SPIE Vol. 5375 (SPIE, Bellingham, WA, 2004)
    
  !  The approximation of a grating as an uniaxial crystal is only 
  !  valid if no higher diffraction orders appear
  ! p: pitch
  ! lambda : wavelength
  !                   eps0: light side
  !               \ theta | 
  !                \      |                              
  !                 \     |                             coordinate system          
  !                  \    |  
  !           < lw   ><-----pitch---->                   -------> x
  !            ------          ------                    |              
  !           |  e1  |  e2    |      |                   | 
  !      -----        --------        ------             z   
  !  
  !                   epss : substrate 
  !
  !  
  ! e1 is the material that has linewidth
  !
  !  
!   eq 3 and 4 of the paper 
!   
!      [eps_eo   0         0  ]
!  e = [  0     eps_or     0  ]
!      [  0      0      eps_or]
!
!   
   
   use globalpar,only:dp
   implicit NONE
   real(dp),intent(in):: theta,pitch,linewidth,lambda
   complex(dp),intent(in):: eps1,eps2,eps0
   complex(dp),intent(out):: eps_or,eps_eo
  ! Local
   real(dp):: f,relp2,sth
   complex(dp):: e0,a0,beta2,Apar,Bpar,Cpar,Dpar,gammatm2
   real(dp):: pi
   !
   pi = 4.d0*datan(1.d0)
   sth = sin(theta)
   f = linewidth/pitch 
   relp2 = (pitch/lambda)**2
   e0 = f*eps2 + (1.d0-f)*eps1 
   a0 = (eps1*f + (1.d0-f)*eps2)/eps1/eps2

   eps_or = e0 + relp2 * (eps2-eps1)**2 *pi*pi/3.d0*(f*f *(1.d0-f)**2)

   beta2 = sqrt(eps0)*sth*sth

   
   Dpar = pi**4 *f*f*(1.d0-f)*(1.d0-f)/3.d0

   Cpar = e0 - beta2 + relp2* beta2*beta2 * (eps2-eps1)**2 * Dpar/pi/pi/e0/e0

   Bpar = 2.d0*beta2*(1.d0/eps2 - 1.d0/eps1)*(eps2-eps1)*Dpar*relp2/pi/pi - e0*a0
   
   Apar = relp2*e0*e0*(1.d0/eps2-1.d0/eps1)**2*Dpar/pi/pi
   
   gammatm2 = (-Bpar-sqrt(Bpar*Bpar-4.d0*Apar*Cpar))/2.d0/Apar
   
   eps_eo = gammatm2*eps_or/(eps_or-beta2)


end subroutine getOpticalParametrers
 

subroutine TRanisoml(nlay,nlmax,eps,epsanis,mu,muanis,lisotropic,thickness,theta, &
 &  fi,zval,izval)
  use globalpar,only: dp,igd
  implicit NONE
   integer,intent(in):: nlmax
  complex(dp),intent(in):: eps(nlmax),mu(nlmax),epsanis(3,3,nlmax),muanis(3,3,nlmax)
  real(dp),intent(in):: thickness(nlmax)
  logical,intent(in):: lisotropic(nlmax)
  real(dp),intent(in)::theta,fi,zval
  integer,intent(in):: nlay
  ! Local variables
   complex(dp) :: QQII(2, 2), QQIII(2, 2)
   complex(dp) :: QIL(2, 2), QIIL(2, 2), QIIIL(2, 2), QIVL(2, 2)
   complex(dp) :: QIR(2, 2), QIIR(2, 2), QIIIR(2, 2), QIVR(2, 2)
   complex(dp):: epst1(3,3),epst2(3,3),mut1(3,3),mut2(3,3)
   complex(dp):: alfa1,alfa2,vita1,vita2,kappa0
   logical:: lisotropic1,lisotropic2 
   real(dp):: g(2,igd)
   complex(dp) :: EINCID(2)
   real(dp) :: dl(3),dr(3),d
   complex(dp):: kapin,kapout,eps1,mu1,eps2,mu2
   real(dp):: ak(2)
   integer:: itype,igmax,igkmax,icomp,izval
   real(dp):: trans,refle,absor,R11,R12,R21,R22
 

   
   eincid(1)=1.d0
   eincid(2) = 0.d0
   kapin = zval*sqrt(eps(1)*mu(1))
   kapout = zval*sqrt(eps(nlay)*mu(nlay))
   kappa0 = zval !............................................
   AK(1)=REAL(KAPIN, kind=DP)*SIN(THETA)*COS(FI)  
   AK(2)=REAL(KAPIN, kind=DP)*SIN(THETA)*SIN(FI)
   IGMAX = 1
   G(1,1)=0.d0
   G(2,1)=0.d0
   
   
   dl(1:3) = 0.d0
   dr(1:3) = 0.d0
   igmax= 1
   igkmax  = 2
   
 loop : do icomp = 1,nlay-1

    lisotropic1 = lisotropic(icomp)
    epst1(1:3,1:3)= epsanis(1:3,1:3,icomp)
    mut1(1:3,1:3)=muanis(1:3,1:3,icomp)
    eps1 = eps(icomp)
    mu1 = mu(icomp)
    alfa1 = 0.d0
    vita1 = 0.d0
    dl(3)  = thickness(icomp)
    !
    lisotropic2= lisotropic(icomp+1)
    epst2 = epsanis(1:3,1:3,icomp+1)
    mut2  = muanis(1:3,1:3,icomp+1)
    eps2  = eps(icomp+1)
    mu2   = mu(icomp+1)
    alfa2 = 0.d0
    vita2 = 0.d0
    dr(3) = 0.d0
    d = 0.d0
    itype=3
 
    call GetLayerQmat(itype,ak,kappa0,QIR,QIIR,QIIIR,QIVR,   &
        &             zval,lisotropic1,epst1,mut1,eps1,mu1,alfa1,vita1,dl,     &
        &                  lisotropic2,epst2,mut2,eps2,mu2,alfa2,vita2,dr,d)
           
            if (icomp==1) then 
                 QIL   = QIR
                 QIIL  = QIIR
                 QIIIL = QIIIR 
                 QIVL  = QIVR
            else
               CALL PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
            end if
         end do loop
      !         [Rpp  Rps ]  
      !     R = [Rsp  Rss ]
      !    
         R11 = QIIIL(1,1)*dconjg(QIIIL(1,1))
         R12 = QIIIL(1,2)*dconjg(QIIIL(1,2))
         R21 = QIIIL(2,1)*dconjg(QIIIL(2,1))
         R22 = QIIIL(2,2)*dconjg(QIIIL(2,2))

         
         call SCAT(nlay-1,IGMAX,ZVAL,AK,KAPIN,KAPOUT,                &
     &                EINCID,QIL,QIIIL,QIIL,QIVL,igkmax,                        &
     &                TRANS,REFLE,ABSOR,izval)

end subroutine TRanisoml

   subroutine GetLayerQmat(itype,ak,kappa0,QI,QII,QIII,QIV,   &
        &                        zval,lisotropic1,epst1,mut1,eps1,mu1,alfa1,vita1,dl,     &
        &                             lisotropic2,epst2,mut2,eps2,mu2,alfa2,vita2,dr,d)
    use globalpar,only:dp,igd   
    implicit none
    integer,intent(in):: itype 
    real(kind=dp),intent(in):: zval
    complex(kind=dp),intent(in)::kappa0
    real(dp),intent(in)::ak(2),d,dl(3),dr(3)
    COMPLEX(kind=DP),intent(out) :: QI(2,2),QII(2,2)
    complex(kind=dp),intent(out) :: QIII(2,2),QIV (2,2)
    complex(kind=dp),intent(inout) :: epst1(3,3),epst2(3,3)
    complex(kind=dp),intent(inout) :: eps1,mu1,eps2,mu2,alfa1,alfa2,vita1,vita2
    logical,intent(in):: lisotropic1,lisotropic2
    real(dp) :: G(2,igd)
    complex(kind=dp):: kappa,kappal,kappar,kappasl,kapl,kapr,mu3,eps3
    complex(kind=dp),intent(inout) :: mut1(3,3),mut2(3,3)
    character(len=6) :: styp
    real(dp):: emach
    integer:: igkmax,igmax,i
    emach = 1.d-16
    igkmax = size( QI, 1 ) 
!      
! Interface with an anisotropic material
!
    igmax = 1
    g(1,1)=0.d0
    g(2,1)=0.d0
    if (itype==3) then 
 
       if (lisotropic1) then 
          do i=1,3
             EPST1(i,i) = EPS1
             MUT1(i,i) = MU1
          end do
          alfa1 = 0.d0    
          vita1 = 0.d0   
       end if
       if (lisotropic2) then 
          do i=1,3
             EPST2(i,i) = EPS2
             MUT2(i,i) = MU2
          end do
          alfa2 = 0.d0   
          vita2 = 0.d0  
       end if

!
! Calculate Q-matrices
!       
     call AINTERF(ZVAL,AK,EPST1,MU1,dl,lisotropic1, &
    &             ALFA1,VITA1,EPST2,MU2,DR,           &
    &             lisotropic2,ALFA2,VITA2,QI,QII,QIII,QIV) 
      else IF(itype.EQ.1) THEN
!
!  Calculate for a slab (all isotropic materials)
!  

         KAPPAL = SQRT( MU1 * EPS1 ) * KAPPA0
         if (Real(MU1,dp)<0.and. Real(EPS1,dp)<0) KAPPAL = -KAPPAL
         KAPPASL= SQRT( MU2 * EPS2 ) * KAPPA0
         if (Real(MU2,dp)<0.and. Real(EPS2,dp)<0) KAPPASL = -KAPPASL  
         KAPPAR = SQRT( MU3 * EPS3 ) * KAPPA0  
         if (Real(MU3,dp)<0.and. Real(EPS3,dp)<0) KAPPAR = -KAPPAR
         
         KAPR = KAPPAR  
         kapl = kappal
         CALL HOSLAB(IGMAX, KAPPAL, EPS1 , MU1,  &
     &                 KAPPASL, EPS2, MU2,        &
     &                 KAPPAR, EPS3, MU3,         &
     &                 AK, G, DL, DR,             &
     &                 D, QI,QII,QIII,QIV, EMACH)
          


      end if
    end subroutine GetLayerQmat

    SUBROUTINE HOSLAB(IGMAX,KAPPA1,EPS1,MU1,KAPPA2,EPS2,MU2,KAPPA3,         &
     &                  EPS3,MU3,AK,G,DL,DR,D,QI,QII,QIII,QIV,EMACH)
      use globalpar, ONLY: IGD,DP,czero,cone,ci,ctwo
      !use constants 
      IMPLICIT NONE 
!C-----------------------------------------------------------------------  
!C     THIS SUBROUTINE CALCULATES THE  Q-MATRICES FOR A HOMOGENEOUS  
!C     PLATE  '2' OF THICKNESS 'D', HAVING THE SEMI-INFINITE MEDIUM  
!C     '1' ON ITS LEFT AND THE SEMI-INFINITE MEDIUM '3' ON ITS RIGHT
!                                    Modified on 22.2.2008 to deal with 
!                                    mu <> 1
!                                    Updated on 13Oct2020       
!C     ------------------------------------------------------------------  
!C                                  
!C  
!C  .. SCALAR ARGUMENTS ..  
!C  
      INTEGER,intent(in)          :: IGMAX  
      REAL(kind=DP),intent(in)    :: EMACH,D   
      COMPLEX(kind=DP),intent(in) :: KAPPA1,KAPPA2,KAPPA3 
      COMPLEX(kind=DP),intent(inout) :: EPS1,MU1,EPS2,MU2,EPS3,MU3 
!C  
!C  .. ARRAY AGUMENTS ..  
!C  
      REAL(kind=DP),intent(in)     :: AK(2),G(2,IGD),DL(3),DR(3)  
      COMPLEX(kind=DP),intent(out) :: QI(:,:),QII(:,:),QIII(:,:)  
      COMPLEX(kind=DP),intent(out) :: QIV(:,:)  
  
!C  
!C  .. LOCAL SCALARS ..  
!C   
      INTEGER          :: I,J,IA,IB,JA,IG1,IGKMAX,it  
      REAL(kind=DP)    :: GKKPAR 
      COMPLEX(kind=DP) :: GKKZ1,GKKZ2,GKKZ3,Z1,Z2,Z3,CQI,CQII  
      COMPLEX(kind=DP) :: CQIII,CQIV,DENOMA,DENOMB,GKKDUM,mudum 
      
!C  
!C  .. LOCAL ARRAYS ..  
!C  
      COMPLEX(kind=DP) :: T(4,2),R(4,2),X(4),P(4,2),e(4),m(4) 
!C  
!C  .. INTRINSIC FUNCTIONS ..  
!C  
      INTRINSIC SQRT,EXP  

!C     -----------------------------------------------------------------  
!C 
      it = size(QI,1) 
      if (it<2*igmax) stop 'Error 1 hoslab : Dimensions' 
      IGKMAX=2*IGMAX  
      DO 1 IB=1,IGKMAX  
      DO 1 IA=1,IGKMAX  
      QI  (IA,IB)=CZERO  
      QII (IA,IB)=CZERO  
      QIII(IA,IB)=CZERO  
      QIV (IA,IB)=CZERO  
    1 CONTINUE  
      X(1)=KAPPA1/KAPPA2
      X(2)=CONE/X(1)  
      X(3)=KAPPA2/KAPPA3
      X(4)=CONE/X(3)  
!
      e(1)=cdsqrt(eps1/eps2)
      e(2)=cone/e(1)
      e(3)=cdsqrt(eps2/eps3)
      e(4)=cone/e(3)
      m(1)=cdsqrt(mu1/mu2)
      m(2)=cone/m(1)
      m(3)=cdsqrt(mu2/mu3)
      m(4)=cone/m(3)
!    
      DO 3 IG1=1,IGMAX  
      GKKPAR=SQRT((AK(1)+G(1,IG1))*(AK(1)+G(1,IG1))+                                   &
     &            (AK(2)+G(2,IG1))*(AK(2)+G(2,IG1)))  
      GKKZ1=SQRT(KAPPA1*KAPPA1-GKKPAR*GKKPAR)  
      GKKZ2=SQRT(KAPPA2*KAPPA2-GKKPAR*GKKPAR)  
      GKKZ3=SQRT(KAPPA3*KAPPA3-GKKPAR*GKKPAR)  

      DO 9 J=1,2 
      DENOMA=X(J)*X(J)*GKKZ2+GKKZ1 
      DENOMB=     GKKZ2+GKKZ1

      IF(ABS(DENOMA).LT.EMACH.OR.ABS(DENOMB).LT.EMACH) GO TO 20 

      R(J,1)=(GKKZ1-X(J)*X(J)*GKKZ2)/DENOMA 
      R(J,2)=           (GKKZ1-GKKZ2)/DENOMB
      T(J,1)=CTWO*X(J)*GKKZ1/DENOMA 
      T(J,2)=CTWO*GKKZ1/DENOMB

      GKKDUM=GKKZ1 
      GKKZ1 =GKKZ2 
      GKKZ2 =GKKDUM 
!
      mudum=mu2
      mu2=mu1
      mu1=mudum
!
 9    CONTINUE 
      DO 10 J=3,4 
      DENOMA=X(J)*X(J)*GKKZ3+GKKZ2 
      DENOMB=          GKKZ3+GKKZ2

      IF(ABS(DENOMA).LT.EMACH.OR.ABS(DENOMB).LT.EMACH) GO TO 20 
      R(J,1)=(GKKZ2-X(J)*X(J)*GKKZ3)/DENOMA 
      R(J,2)=           (GKKZ2-GKKZ3)/DENOMB 
      T(J,1)=CTWO*X(J)*GKKZ2/DENOMA 
      T(J,2)=CTWO*GKKZ2/DENOMB 
      GKKDUM=GKKZ2 
      GKKZ2 =GKKZ3 
      GKKZ3 =GKKDUM

      mudum=mu3
      mu3=mu2
      mu2=mudum

 10   CONTINUE 
      Z1=EXP(CI*GKKZ2*D)  
      Z2=Z1*Z1  
      DO 5 I=1,2  
      Z3=CONE/(CONE-Z2*R(2,I)*R(3,I))  
      P(1,I)=T(3,I)*Z3*Z1*T(1,I)  
      P(2,I)=R(4,I)+T(4,I)*R(2,I)*T(3,I)*Z2*Z3  
      P(3,I)=R(1,I)+T(2,I)*R(3,I)*T(1,I)*Z2*Z3  
      P(4,I)=T(2,I)*Z3*Z1*T(4,I)  
    5 CONTINUE  
      CQI  =EXP(CI*((AK(1)+G(1,IG1))*(DL(1)+DR(1))+                                    & 
     &              (AK(2)+G(2,IG1))*(DL(2)+DR(2))+                                    &
     &	             GKKZ1*DL(3)+GKKZ3*DR(3)))  
      CQII =EXP(CTWO*CI*GKKZ3*DR(3))  
      CQIII=EXP(CTWO*CI*GKKZ1*DL(3))  
      CQIV =EXP(-CI*((AK(1)+G(1,IG1))*(DL(1)+DR(1))+                                   &
     &               (AK(2)+G(2,IG1))*(DL(2)+DR(2))-                                   & 
     &	             GKKZ1*DL(3)-GKKZ3*DR(3)))  
      DO 7 JA=1,2  
      IA=2*IG1-2+JA  
      QI  (IA,IA)=CQI  *P(1,JA)  
      QII (IA,IA)=CQII *P(2,JA)  
      QIII(IA,IA)=CQIII*P(3,JA)  
      QIV (IA,IA)=CQIV *P(4,JA)  
    7 CONTINUE  
    3 CONTINUE
      RETURN 
   20 STOP 'FATAL ERROR IN HOSLAB' 
      END subroutine 
SUBROUTINE AINTERF(ZVAL,AK,EPS1,MU1,DL,Lisotropic1,ALFA1,VITA1, &
                  & EPS2,MU2,DR,Lisotropic2,ALFA2,VITA2,QI,QII,QIII,QIV)
  ! This subroutine returns the transmission and reflection matrices 
  ! of an interface between two homogeneous, in general bi-anisotropic
  ! media, characterized by epsilon, mu, alfa and vita tensor. The 
  ! transmission and reflection matrices refer to a center at dl(dr)
  ! of the interface an the 1(2) medium.
  !                                        
  use globalpar, only:DP,igd,czero
 
 IMPLICIT NONE
 INTEGER         ::IGMAX
 REAL(kind=DP)   ,INTENT(IN) ::ZVAL,AK(2),DL(3),DR(3)!,G(2,igd) ! 2,igd solve this
 COMPLEX(kind=DP),INTENT(IN) ::EPS1(3,3),EPS2(3,3),MU1,MU2,alfa1,alfa2,vita1,vita2
 LOGICAL         ,INTENT(IN) ::LISOTROPIC1,LISOTROPIC2
 COMPLEX(kind=DP),INTENT(OUT)::QI(2,2),QII(2,2),QIII(2,2),QIV(2,2)  
 
 COMPLEX(kind=DP) ::GKKZ1(4),GKKZ2(4),bndmatrix1(4,4),bndmatrix2(4,4)
 COMPLEX(kind=DP) ::SI(2,2),SII(2,2),SIII(2,2),SIV(2,2),ep1(3,4),ep2(3,4)
 INTEGER          ::IA,IB,IG1,JA,JB,IG,i,j,IG2,igkmax
 REAL(kind=DP)    ::GKKX,GKKY,GKKPAR(2)
 real(kind=dp)    :: g(2,1)
 
 igkmax = 2       
 QI  =CZERO  
 QII =CZERO  
 QIII=CZERO  
 QIV =CZERO   

  IGMAX = 1
  G(1,1)= 0.d0
  G(2,1)= 0.d0

 DO IG1=1,IGMAX
    GKKX=AK(1)+G(1,IG1)
    GKKY=AK(2)+G(2,IG1)
    GKKPAR(1)=AK(1)+G(1,IG1)
    GKKPAR(2)=AK(2)+G(2,IG1)
    call emvec6x6(zval,GKKPAR,eps1,mu1,lisotropic1,alfa1,vita1,GKkz1,ep1,bndmatrix1)
    call emvec6x6(zval,GKKPAR,eps2,mu2,lisotropic2,alfa2,vita2,GKkz2,ep2,bndmatrix2)
    CALL SMATRIXES (IGMAX,GKKX,GKKY,bndmatrix1,bndmatrix2,GKkz1,GKkz2, &
                    & dL,dR,si,sii,siii,siv)
    IG=2*IG1-2
    DO JA=1,2
       IA=IG+JA
       DO JB=1,2
          IB=IG+JB
          QI  (IA,IB) =SI(JA,JB)
          QII (IA,IB) =SII(JA,JB)
          QIII(IA,IB) =SIII(JA,JB)
          QIV (IA,IB) =SIV(JA,JB)
       END DO
    END DO
 END DO
END SUBROUTINE AINTERF
!=====================================================================================
Subroutine Smatrixes(IGMAX,GKKX,GKKY,bndmatrix1,bndmatrix2&
                     &,kz01,kz02,dL,dR,si,sii,siii,siv)
  ! This subroutine computes the transmission and reflection matrices 
  ! of an interface between two homogeneous, in general bi-anisotropic
  ! media, characterized by epsilon, mu, alfa and vita tensor, after
  ! applying the boundary conditions of continuity of the tangential
  ! components of the EM wave.
 use globalpar, only: DP,pi,czero,cone,ci
 implicit none  
 
 INTEGER         ,INTENT(IN) ::IGMAX
 real(kind=DP)   ,intent(in) ::dL(3),dR(3),GKKX,GKKY
 complex(kind=DP),intent(in) ::bndmatrix1(4,4),bndmatrix2(4,4)
 complex(kind=DP),intent(in) ::kz01(4),kz02(4)
 complex(kind=DP),intent(out)::si(2,2),sii(2,2),siii(2,2),siv(2,2)
 
 complex(kind=DP) ::rtmatrix(4,4),rtvector(4,4),EXP_GKKPAR_P,EXP_GKKPAR_M
 integer          ::i,j,intar(4),info
 real(kind=DP)    ::EMACH

 
 rtvector=czero
 rtmatrix=czero
 rtmatrix(1:4,1:2)= bndmatrix1(1:4,1:2)
 rtmatrix(1:4,3:4)=-bndmatrix2(1:4,3:4)
 rtvector(1:4,1:2)=-bndmatrix1(1:4,3:4)
 rtvector(1:4,3:4)= bndmatrix2(1:4,1:2)

 call zgetrf(4,4,rtmatrix,4,intar,info) 
 call zgetrs('N',4,4,rtmatrix,4,intar,rtvector,4,info) 
 
 SI  (1:2,1:2)=rtvector(3:4,1:2)
 SIII(1:2,1:2)=rtvector(1:2,1:2)
 SII (1:2,1:2)=rtvector(3:4,3:4)
 SIV (1:2,1:2)=rtvector(1:2,3:4)
 
 EXP_GKKPAR_P=EXP( CI*(GKKX*(DR(1)+DL(1))+GKKY*(DR(2)+DL(2))))
 EXP_GKKPAR_M=EXP(-CI*(GKKX*(DR(1)+DL(1))+GKKY*(DR(2)+DL(2))))
 
 DO I=1,2
    si(1,i)  =SI(1,i)  *exp( ci*( kz02(3)*dr(3)+kz01(2+I)*dl(3) ))*EXP_GKKPAR_P
    si(2,i)  =SI(2,i)  *exp( ci*( kz02(4)*dr(3)+kz01(2+I)*dl(3) ))*EXP_GKKPAR_P
    siii(1,i)=siii(1,i)*exp(-ci*( kz01(1)*dl(3)-kz01(2+I)*dl(3) ))
    siii(2,i)=siii(2,i)*exp(-ci*( kz01(2)*dl(3)-kz01(2+I)*dl(3) ))
    siv(1,i) =sIV(1,i) *exp(-ci*( kz01(1)*dl(3)+kz02(I)*dr(3) ))*EXP_GKKPAR_M
    siv(2,i) =SIV(2,i) *exp(-ci*( kz01(2)*dl(3)+kz02(i)*dr(3) ))*EXP_GKKPAR_M
    sii(1,i) =SII(1,i) *exp( ci*( kz02(3)*dr(3)-kz02(I)*dr(3) ))
    sii(2,i) =SII(2,i) *exp( ci*( kz02(4)*dr(3)-kz02(I)*dr(3) ))
 END DO
  
End subroutine smatrixes
!=====================================================================================
Subroutine emvec6x6(omega,ak,eps_mat,mu,lisotropic,alfa,vita,qz0,ep0,boundmatrix)
 ! This subroutine computes the eigenvectors and eigenvalues of an infinite optically 
 ! in general bi-anisotropic medium given an eigenvalue-eigenvector system of 
 ! the Maxwell Equations. Then the boundary matrixes are created from the x and y components 
 ! of the EM wave.
 ! . . . . . . Eigenvalue-Eigenvector System . . . . . .
 !   Maxwell Equations, if plane wave solutions characterized by Q and w are considered,
 !   in the form 
 !   
 !   ( 0    -Qx ) (E)      ( ε   ζ ) (Ε)
 !   | =     =  | |-| -- w | =   = | |-|
 !   | Qx    0  | |H| -- c | ξ   μ | |Η|   , or
 !   ( =     =  ) (-)      ( =   = ) (-)
 !   
 !    (Qpar_mat +qz*Qz_mat)*EH_vec=w/c*(EM_mat)*EH_vec
 !
 !   where == denotes matrix, - denotes vector, and Qx corresponds to cros product of Q vector
 !   with E or H vector. Assuming that the stratified structure is along the z axis, we write
 !   Q = Q_par + q_z*z, where Q, Q_par, z are vectors and q_z is the wanted eigenvalue.
 !   Finally the eigenvalue-eigenvector system is of the form
 !   
 !   [w/c*(EM_mat)-Qpar_mat]^(-1)*Qz_mat*EH_vec=qz^(-1)*EH_vec
 !   10/06/2017 by pepantaz
  use globalpar, only: DP ,ci,czero,cone,pi
  
 implicit none
 
 real(kind=DP)   ,intent(in) ::ak(2),omega
 complex(kind=DP),intent(in) ::eps_mat(3,3),alfa,vita,mu
 logical         , intent(in)::lisotropic
 complex(kind=DP),intent(out)::qz0(4),ep0(3,4),boundmatrix(4,4)
 
 integer::i,j,m,intar(6),info,icum,ifail,info1
 real(kind=DP)::emach,rwork(2*6),epnorm,akpar,waveimp
 complex(kind=DP)::mu_mat(3,3),ksi_mat(3,3),zeta_mat(3,3),hp0(3,4),work(2*6)
 complex(kind=dp)::qnorm,st,cf,ct,e1(3),e2(3),e0(3),xq,swap,emswap(3),dcheck
 complex(kind=DP)::qz_mat(3,3),qpar_mat(3,3),dp0(3,4),bp0(3,4),xep(3),xhp(3)
 complex(kind=DP)::qz_lrg(6,6),qz(6),em_lrg(6,6),qpar_lrg(6,6),eh(6,6),dum(6,6)
 
 real(kind=DP),external::dznrm2

 !WaveImp=376.73d0
 WaveImp = 376.730313668d0 ! Wikipedia
   !print*,"eps",eps_mat
If (.not.lisotropic) then

  !print*,"eps",eps_mat
   mu_mat=czero
   ksi_mat=czero
   zeta_mat=czero
   do i=1,3
      mu_mat(i,i)=cone*mu
      ksi_mat(i,i)=ci*alfa
      zeta_mat(i,i)=ci*vita
   end do
   
   qz_mat=czero
   qpar_mat=czero
   
   qz_mat(1,2)=-cone
   qz_mat(2,1)= cone
   
   qpar_mat(1,3)= ak(2)
   qpar_mat(2,3)=-ak(1)
   qpar_mat(3,1)=-ak(2)
   qpar_mat(3,2)= ak(1)
   
   qz_lrg=czero
   qpar_lrg=czero
   em_lrg=czero
   
   qz_lrg(1:3,4:6)=-qz_mat(1:3,1:3)
   qz_lrg(4:6,1:3)= qz_mat(1:3,1:3)
   
   qpar_lrg(1:3,4:6)=-qpar_mat(1:3,1:3)
   qpar_lrg(4:6,1:3)= qpar_mat(1:3,1:3)
   
   em_lrg(1:3,1:3)= eps_mat(1:3,1:3)
   em_lrg(1:3,4:6)=-zeta_mat(1:3,1:3)
   em_lrg(4:6,1:3)= ksi_mat(1:3,1:3)
   em_lrg(4:6,4:6)= mu_mat(1:3,1:3)
   
   em_lrg=omega*em_lrg
   em_lrg=em_lrg-qpar_lrg
   
   EMACH=1.D-16
   
   call zgetrf(6,6,em_lrg,6,intar,info1) 
   call zgetrs('N',6,6,em_lrg,6,intar,qz_lrg,6,info1) 
   call zgeev ('N','V',6,qz_lrg,6,qz,dum,6,eh,6,work,2*6,rwork,info)
   if (info.ne.0) print*,"Error at zgeev; 6x6"
    icum=0
    do i=1,6
       if (cdabs(qz(i))>1.d-12) then
          icum = icum + 1
          qz0(icum) = cone/qz(i)
          ep0(1:3,icum) = eh(1:3,i)
          hp0(1:3,icum) = eh(4:6,i)
        end if
    end do
    
    if (icum.ne.4) then 
        write(6,*) 'ERROR at emvec6x6: Eigenvalues are not 4 please check'
        print*, "qx=",ak(1)
        print*, "qx=",ak(2)
        print*, "w=",omega
        STOP 
    end if   
    do i=1,4
       epnorm=1.d0/dznrm2(3,ep0(1,i),1)
       ep0(1:3,i)=ep0(1:3,i)*epnorm
       hp0(1:3,i)=hp0(1:3,i)*epnorm/waveimp
    end do
    
    !qz0 sorting based on the real part of it. so as abs(kz0(1)>kz0(2))
    do j=2,4
        Xq=qZ0(J)
        XEP(1:3)=EP0(1:3,J)
        xhp(1:3)=hp0(1:3,j)
        if (dabs(DREAL(qZ0(j))).gt.1.D-8) then
           do m=j-1,1,-1
              if (DREAL(qZ0(m)).le.DREAL(Xq)) goto 1
              qZ0(m+1)=qZ0(m)
              EP0(1:3,m+1)=EP0(1:3,m)
              hP0(1:3,m+1)=hP0(1:3,m)          
           end do 
           m=0
 1         qZ0(m+1)=Xq
           EP0(1:3,M+1)=XEP
           hP0(1:3,M+1)=XhP
        else 
           do m=j-1,1,-1
              if (DIMAG(qZ0(m)).le.DIMAG(Xq)) goto 2
              qZ0(m+1)=qZ0(m)
              EP0(1:3,m+1)=EP0(1:3,m)
              hP0(1:3,m+1)=hP0(1:3,m)
           end do 
           m=0
 2         qZ0(m+1)=Xq
           EP0(1:3,M+1)=XEP
           hP0(1:3,M+1)=XhP
        end if
    end do
    swap=qz0(4)
    qz0(4)=qz0(3)
    qz0(3)=swap
    
    emswap(1:3)=ep0(1:3,4)
    ep0(1:3,4)=ep0(1:3,3)
    ep0(1:3,3)=emswap(1:3)
    
    emswap(1:3)=hp0(1:3,4)
    hp0(1:3,4)=hp0(1:3,3)
    hp0(1:3,3)=emswap(1:3)
    
    boundmatrix(1:2,1:4)=ep0(1:2,1:4)
    boundmatrix(3:4,1:4)=hp0(1:2,1:4)
 
    do i=1,4
       dp0(1:3,i)=matmul(eps_mat,ep0(1:3,i))+ci*matmul(-zeta_mat,hp0(1:3,i))
       bp0(1:3,i)=matmul(mu_mat,hp0(1:3,i))+ci*matmul(ksi_mat,ep0(1:3,i))
    end do   
 
    do I=1,4
       dcheck=ak(1)*dp0(1,i)+ak(2)*dp0(2,i)+qz0(i)*dp0(3,i)
       if (cdabs(dcheck).gt.1.d-8) print*,">eps=anis-6x6;wave not transverse at w=",omega,"dcheck=",cdabs(dcheck)
    end do
else
   qz0(1:2)=- cdsqrt(omega*omega*eps_mat(1,1)*mu-ak(1)*ak(1)-ak(2)*ak(2))
   qz0(3:4)=- qz0(1:2)
   qnorm= omega*cdsqrt(eps_mat(1,1)*mu)
   AKPAR=SQRT(ak(1)*ak(1)+ak(2)*ak(2))  
   ST=AKPAR/qnorm  
   CF=CONE  
   IF(AKPAR.GT.1.D-8) CF=CMPLX(ak(1)/AKPAR,ak(2)/AKPAR,kind=dp) 
   do i=1,4
      CT=qz0(i)/qnorm  
      e0(1)=st*dreal(cf)
      e0(2)=st*dimag(cf)
      e0(3)=dreal(cf)
      e1(1)=ct*dreal(cf)
      e1(2)=ct*dimag(cf)
      e1(3)=-st
      e2(1)=dimag(cf)
      e2(2)=-dreal(cf)
      e2(3)=czero
      if (i.eq. 1 .or. i.eq. 3) then ! tm polarisation
         ep0(1:3,i)=e1(1:3)
         dp0(1:3,i)=matmul(eps_mat,e1)
      else ! te polarisation
         ep0(1:3,i)=e2(1:3)
         dp0(1:3,i)=matmul(eps_mat,e2)
      end if
      hp0(1,i)=(ak(2)*ep0(3,i)-qz0(i)*ep0(2,i))/(omega*mu*WaveImp)
      hp0(2,i)=(-ak(1)*ep0(3,i)+qz0(i)*ep0(1,i))/(omega*mu*WaveImp)
      hp0(3,i)=(ak(1)*ep0(2,i)-ak(2)*ep0(1,i))/(omega*mu*WaveImp) 
      dcheck=ak(1)*dp0(1,i)+ak(2)*dp0(2,i)+qz0(i)*dp0(3,i)
      if (cdabs(dcheck).gt. 1.d-8) print*,"6x6 wave not transverse at w="&
                                          &,omega,i,lisotropic,dcheck,eps_mat
      boundmatrix(1:2,i) = ep0(1:2,i)
      boundmatrix(3:4,i) = hp0(1:2,i)
   end do
 end if    
end Subroutine emvec6x6


SUBROUTINE SCAT(ncomp,IGMAX,ZVAL,AK,KAPIN,KAPOUT,                &
     &                EINCID,QI,QIII,QII,QIV,igkmax,                        &
     &                TRANS,REFLE,ABSOR,izval)
      use globalpar   
      IMPLICIT NONE 

!C     ------------------------------------------------------------------  
!C     THIS SUBROUTINE CALCULATES THE REFLECTIVITY, TRANSMITTANCE AND  
!C     ABSORBANCE OF A FINITE SLAB, CHARACTERIZED BY TRANSMISSION AND 
!C     REFLECTION MATRICES QI AND QIII, RESPECTIVELY
!C     ------------------------------------------------------------------  
!C ..  SCALAR ARGUMENTS  ..  
!C      
      INTEGER         ::  IGMAX,izval,ncomp 
      REAL(kind=DP)   ::  ZVAL,KAPIN,KAPOUT   
!C  
!C ..  ARRAY ARGUMENTS ..  
!C  
      REAL(kind=DP)     :: AK(2),G(2,IGd),AL(3),AR(3)  
      COMPLEX(kind=DP),intent(in)  :: QI(igkmax,igkmax),QIII(igkmax,igkmax),EINCID(igkmax)
      COMPLEX(kind=DP),intent(in)  :: QII(igkmax,igkmax),QIV(igkmax,igkmax)  
!C  
!C ..  LOCAL SCALARS  ..  
!C  
      INTEGER           :: IGK1,IG1,K1,IGK2,ig2,IGKMAX,igone,ifile,ier,igtwo,ii
      integer           :: info,i,igk,icomp,k2,ig0,ig10,ig01,igm10,ig0m1,id
      REAL(kind=DP)     :: DOWN,REFLE,TRANS,ABSOR,GKZIN,GKZOUT,TES1,tr2,rf2,DOWNN
      REAL(kind=DP)     :: slabthick,scaleunits ,deltan,dos,ddos,ddos1,ddos2,ddos3 
      COMPLEX(kind=DP)  :: tr,rf,inc,t,r,epseff,mueff,neff,zeff
      COMPLEX(kind=dp) :: dot1,dot2,cqi,cqii,cqiii,cqiv   
      character(len=100):: uuio
      logical           :: run
      real(kind=dp)     :: phi,chi,kapa,thetaF,chir 
      COMPLEX(kind=dp) ::  chi1,ETRANSL,ETRANSR
     real(kind=dp) :: Transm_x_pol,Transm_y_pol,normEx 
     real(kind=dp) ::  normEy,phase_x,phase_y,phasedif,tr_Lpol,tr_Rpol,TRANS2 
     character(len=4):: s1,s2
     character(len=20) :: lab0,lab
     integer:: ig,itime
!C  
!C ..  LOCAL ARRAYS  ..  
!C
      real(kind=dp),allocatable:: Tgv(:),tgvn(:), Rgv(:), RgvN(:) 
      complex(kind=dp),allocatable::ETRANS(:),EREFLE(:),RC(:),gkk(:,:),gkkzl(:),gkkzr(:)
      complex(kind=dp),allocatable::SMATRIX(:,:),COMVEC(:,:)
      complex(kind=dp),allocatable::work(:),dummyl(:,:) 
      real(kind=dp),allocatable:: rwork(:)
      complex(kind=dp),allocatable::EINCIDn(:),ETRANSn(:),EREFLEn(:)    
      REAL(kind=DP)     ::  REFLEn,TRANSn,ABSORn ,chiim, eliptangle
      
      
!C  
!C ..  INTRINSIC FUNCTIONS ..  
!C  
      INTRINSIC SQRT,CONJG    
!C     ------------------------------------------------------------------ 
      
      G(1,1)= 0.d0
      G(2,1)= 0.d0
      ifile = 10
      allocate ( Tgv(igkmax), tgvn(igkmax), Rgv(igkmax), RgvN(igkmax) )  
      allocate ( ETRANS(IGKMAX), EREFLE(IGKMAX) )
      allocate( EINCIDn(IGKMAX))
      allocate ( ETRANSn(IGKMAX),EREFLEn(IGKMAX) )
      itime = izval
      DOWN=0.D0  
      DOWNN = 0.d0
      REFLE=0.D0  
      TRANS=0.D0  
      REFLEn=0.D0  
      TRANSn=0.D0 
      IGK1=0 
      TgV = 0.d0
      TgvN = 0.d0
      Rgv = 0.d0
      RgvN = 0.d0
      DO  IG1=1,IGMAX  
         TES1=(AK(1)+G(1,IG1))*(AK(1)+G(1,IG1))+(AK(2)+G(2,IG1))*(AK(2)+                   &
              &            G(2,IG1))
         GKZIN =0.D0  
         GKZOUT=0.D0
         IF(( KAPIN*KAPIN -TES1).GT.0.D0) GKZIN =SQRT( KAPIN*KAPIN -TES1) 
         IF((KAPOUT*KAPOUT-TES1).GT.0.D0) GKZOUT=SQRT(KAPOUT*KAPOUT-TES1)
         do i=1,igmax
             EINCIDn(2*I-1)=-EINCID(2*I)
             EINCIDn(2*I)  = EINCID(2*I-1)
          end do
          DO K1=1,2  
             IGK1=IGK1+1  
             ETRANS(IGK1)=CZERO  
             EREFLE(IGK1)=CZERO 
             ETRANSn(IGK1)=CZERO
             EREFLEn(IGK1)=CZERO     

             DO IGK2=1,IGKMAX
                ETRANS(IGK1)=ETRANS(IGK1)+QI  (IGK1,IGK2)*EINCID(IGK2)  
                EREFLE(IGK1)=EREFLE(IGK1)+QIII(IGK1,IGK2)*EINCID(IGK2)
                ETRANSn(IGK1)=ETRANSn(IGK1)+QI  (IGK1,IGK2)*EINCIDn(IGK2)  
                EREFLEn(IGK1)=EREFLEn(IGK1)+QIII(IGK1,IGK2)*EINCIDn(IGK2) 
             end do

             DOWN =DOWN +EINCID(IGK1)*DCONJG(EINCID(IGK1))*GKZIN
             TRANS=TRANS+ETRANS(IGK1)*DCONJG(ETRANS(IGK1))*GKZOUT  
             REFLE=REFLE+EREFLE(IGK1)*DCONJG(EREFLE(IGK1))*GKZIN
       
             DOWNN = DOWNN +EINCIDN(IGK1)*DCONJG(EINCIDN(IGK1))*GKZIN
             TRANSn=TRANSn+ETRANSn(IGK1)*DCONJG(ETRANSn(IGK1))*GKZOUT  
             REFLEn=REFLEn+EREFLEn(IGK1)*DCONJG(EREFLEn(IGK1))*GKZIN
     
             Tgv(igk1)=ETRANS(IGK1)*DCONJG(ETRANS(IGK1))*GKZOUT 
             TgvN(igk1)=ETRANSn(IGK1)*DCONJG(ETRANSn(IGK1))*GKZOUT  
             Rgv(igk1) =EREFLE(IGK1)*DCONJG(EREFLE(IGK1))*GKZIN
             RgvN(igk1)=EREFLEn(IGK1)*DCONJG(EREFLEn(IGK1))*GKZIN
 
          end do
       end do
      
          TRANS=TRANS/DOWN  
          REFLE=REFLE/DOWN  
          ABSOR=1.D0-TRANS-REFLE 
     
          TRANSn=TRANSn/DOWNN  
          REFLEn=REFLEn/DOWNN  
          ABSORn=1.D0-TRANSn-REFLEn
      
          Tgv(:)=Tgv(:)/DOWN
          TgvN(:)=TgvN(:)/DOWNN
          Rgv(:)=Rgv(:)/DOWN
          RgvN(:)=RgvN(:)/DOWNN
          
          ig0 = 1
       
          WRITE(6,101)  ZVAL,TRANS,REFLE,ABSOR 

          id=get_funid('TRANS_P')
          if (id>0) output_data(itime,id) =TRANS
          id=get_funid('TRANS_S')
          if (id>0) output_data(itime,id) = transn
          id=get_funid('REFLE_P')
          if (id>0) output_data(itime,id) = REFLE
          id=get_funid('REFLE_S')
          if (id>0) output_data(itime,id) = reflen
          id=get_funid('ABSOR_P')
          if (id>0) output_data(itime,id) = ABSOR
          id=get_funid('ABSOR_S')
          if (id>0) output_data(itime,id) = absorn
         
          deallocate ( Tgv, tgvn, Rgv, RgvN )  
          deallocate (ETRANS,EREFLE)
          deallocate( EINCIDn)
          deallocate ( ETRANSn,EREFLEn)
          
          RETURN    
101       FORMAT(3F16.10,25E20.10)
201       format(25E16.8) 
        END subroutine SCAT

     SUBROUTINE PAIR(IGKMAX,QIL,QIIL,QIIIL,QIVL,QIR,QIIR,QIIIR,QIVR)
      use globalpar, ONLY:DP,czero,cone 
      IMPLICIT NONE 
!C     ------------------------------------------------------------------  
!C     THIS SUBROUTINE CALCULATES SCATTERING Q-MATRICES FOR A  DOUBLE  
!C     LAYER, FROM THE CORRESPONDING MATRICES OF THE INDIVIDUAL, LEFT  
!C     (L) AND RIGHT (R), LAYERS. THE RESULTS ARE STORED IN Q*L.  
!                                                   Updated 9Dec21
!C     -----------------------------------------------------------------  
!C  
!C ..  PARAMETER STATEMENTS  ..  
!C    
!C  
!C ..  SCALAR ARGUMENTS  ..  
!C  
      INTEGER,intent(in):: IGKMAX  
!C  
!C ..  ARRAY ARGUMENTS  ..  
!C  
      COMPLEX(kind=DP),intent(inout):: QIL (igkmax,igkmax),QIIL(igkmax,igkmax),QIIIL(igkmax,igkmax)  
      COMPLEX(kind=DP),intent(inout):: QIVL(igkmax,igkmax)  
      COMPLEX(kind=DP),intent(inout):: QIR (igkmax,igkmax),QIIR(igkmax,igkmax),QIIIR(igkmax,igkmax)  
      COMPLEX(kind=DP),intent(inout):: QIVR(igkmax,igkmax)  
!C  
!C ..  LOCAL SCALARS  ..  
!C  
      INTEGER  ::  IGK1,IGK2,IGK3,info
      REAL(kind=DP),parameter::    EMACH=1.d-8    
!C  
!C ..  LOCAL ARRAYS  ..  
!C  
      INTEGER,allocatable::    INT(:),JNT(:)  
      COMPLEX(kind=DP),allocatable:: QINV1(:,:),QINV2(:,:),W1(:,:)  
      COMPLEX(kind=DP),allocatable:: W2(:,:),W3(:,:),W4(:,:) 
!C-----------------------------------------------------------------------  

     
      allocate (QINV1(igkmax,igkmax),QINV2(igkmax,igkmax))
      allocate (W1(igkmax,igkmax),W2(igkmax,igkmax),W3(igkmax,igkmax),W4(igkmax,igkmax))
      allocate (INT(igkmax),jnt(igkmax))
      QINV1 = QIL
      QINV2 = QIVR
      w2 = - matmul(QIIL,QIIIR)
      w3 = - matmul(QIIIR,QIIL)
      do igk1=1,igkmax
         W2(igk1,igk1) = cone +  W2(igk1,igk1) 
         W3(igk1,igk1) = cone +  W3(igk1,igk1)
      end do
      call zgetrf(IGKMAX,IGKMAX,W2,IGKMAX,INT,info)
      call zgetrf(IGKMAX,IGKMAX,W3,IGKMAX,JNT,info)
      CALL ZGETRS('N',IGKMAX,IGKMAX,W2,IGkmax,INT,Qinv1,IGkmax,INFO) 
      CALL ZGETRS('N',IGKMAX,IGKMAX,W3,IGkmax,JNT,Qinv2,IGkmax,INFO) 
      w1 = matmul(QIR,QINV1)
      w2 = matmul(QIIL,QINV2)
      w3 = matmul(QIIIR,QINV1)
      w4 = matmul(QIVL,QINV2)
      qinv1 = QIIR + matmul(QIR,W2)
      qinv2 = QIIIL + matmul(QIVL,W3)
      QIL   = W1    
      QIIL  = QINV1  
      QIIIL = QINV2  
      QIVL  = W4  
      deallocate (Qinv1,qinv2,w1,w2,w3,w4)
      deallocate(int,jnt)
      RETURN  
   END 
