module reflexMod
implicit none
!dllModule_working6.f90 + testdll_working8.f90  (finale e funzionante)
!TODO: mettere controllo dimensioni e layer dispari, con scrittura su file di 
!log in caso di errore.
!mettere materiali come argomento.
!gestione overcoating-> se layer dispari tiene sempre il materiale leggero
!in superficie.

!2011/10/29 Vincenzo
!compiled with:
! f2py -c dllModule_working6.f90 -m reflexf90
!It creates the file reflexf90.so, which can be
! imported in python with 
!import reflexf90
!The reflectivity routine can be called with
!


integer, parameter:: pmaxEn=1000
real*4 ::  ener(pmaxEn)
complex,parameter :: ci=(0.0,1.0)
real(8),parameter::  H_C = 12.39841857D0 
real(8),parameter::  TORAD = 0.017453292519943295769237d0
real(8), parameter ::PI=3.141592653589793238462643D0 !pi greek

contains

SUBROUTINE LOCATE(XX,N,X,j)
	IMPLICIT none
        INTEGER j,n
	REAL(8) x,xx(n)
	INTEGER jl,jm,ju
        
	! Dato un vettore xx(1:n) e dato un valore x ritorna un valore j tale che x sia
	! tra xx(j) e xx(j+1). xx(1:n) deve essere monotonica (crescente o decrescente)
	! j=1 e j=n is ritornata per indicare che x e' fuori dal range
	JL=0
	JU=N+1
	10	IF ((JU-JL).GT.1)THEN
	JM=(JU+JL)/2
	IF((XX(N).GT.XX(1)).EQV.(X.GT.XX(JM)))THEN
		JL=JM
	ELSE
		JU=JM
	ENDIF
	GO TO 10
	ENDIF
	j=JL
	RETURN
END	subroutine

subroutine reflex(dSpacing,n_layers,ener,nener,angle,rough,rs2,re2,ro2,ref)
!versione semplificata della routine originaria: una sola polarizzazione
!solo strati alternati con mat1odd sullo strato esterno.
!dSpacing fornito a partire dal substrato.

implicit none
integer*4, intent(in):: nener,n_layers
real*4,intent(out):: ref(nener)
real*4, intent(in)::rough,ener(nener),angle
integer*4 n_bilayers
integer i,j,t
real*4 dSpacing(n_layers),de(n_layers/2),do(n_layers/2) !n_layers deve essere pari
complex*16,intent(in):: re2(nener), rs2(nener), ro2(nener)

complex*16 fo,fe,fv,ffe,ffv,ffo	
complex*16 ao,ae,r,fs,ffs
complex*16 FcostO,FcostE,costO,costE
complex*16 ro2j, re2j, rs2j
real(8) SIGMA_O2,SIGMA_E2,SIGMA_S2
real(8) seno,coseno,refv,prefact
real(8) pp,qq
real(8) arg_e,arg_o,arg_v,arg_s,fnevot_e,fnevot_o,fnevot_v,fnevot_s
real(8) VAL,r_l
real(8) oddrough,evenrough,subrough

	n_bilayers=int(n_layers/2)
	de=0
	do=0
	DO I=1,N_biLAYERS !i va dal basso			
		DE(I)=Dspacing(i*2-1)
		DO(I)=Dspacing(i*2)
	ENDDO
	
	open(1,file="sonoUnaDll.txt",position='append')
	DO I=1,N_biLAYERS !i va dal basso			
		write(1,*) DE(I),DO(I)
	ENDDO
	close(1)

	oddrough=rough
	evenrough=rough
	subrough=rough
	coseno=cos(angle)
	seno=sin(angle)
	refv=dabs(seno)
	fv = Dcmplx(refv)
	ref=0

	DO J=1,nener
  																																					 			
		r=(0.0,0.0)
		val=ener(j)
		r_l=H_C/ener(j)
		ro2j=ro2(j)
		re2j=re2(j)
		rs2j=rs2(j)
		
		prefact=2*((2*(PI))/(r_l))**2

		fo = ro2j - coseno**2 !sin theta_inc - sin theta_ critical
		fo=cDsqrt(fo)

		fe = re2j - coseno**2
		fe=cDsqrt(fe)

		fs = rs2j - coseno**2
		fs=cDsqrt(fs)

		ffe=(fe-fo)/(fe+fo)! Fresnel formula "S" (in function of incidence                         
		ffo=-ffe							!  angle and critical angle)
		ffs=(fe-fs)/(fe+fs)

		!------------------------------------------------------------------------------
		! Nevot-Croce roughness
		!-----------------------------------------------------------------------------
		sigma_e2=evenrough**2 !roughn. even layer
		sigma_o2=oddrough**2 !roughn. odd layer
		sigma_s2=subrough**2 !roughn. substrate
		!----
		arg_e=FO*FE*sigma_e2
		arg_O=FO*FE*sigma_o2
		arg_v=FO*seno*sigma_o2
		arg_s=FE*FS*sigma_s2
		!---
		FNEVOT_E=DEXP(-prefact*arg_e)
		FNEVOT_O=DEXP(-prefact*arg_o)
		FNEVOT_V=DEXP(-prefact*arg_v)
		FNEVOT_S=DEXP(-prefact*arg_s)

		FcostO=ffe*fnevot_o
		FcostE=ffo*fnevot_e
		costO=-4*ci*pi*fo/r_l
		costE=-4*ci*pi*fe/r_l

		ao=costo*do(1)
		ae=coste*de(1)
		ao=cDexp(ao)
		ae=cDexp(ae)

		r=ae*(r+ffs*fnevot_s)/(r*ffs*fnevot_s+1.0)
		r=ao*(r+ffo*fnevot_e)/(r*ffo*fnevot_e+1.0)

		do T=2,N_BILAYERS
			ao=costO*do(T)
			ae=costE*de(T)
			ao=cDexp(ao)
			ae=cDexp(ae)
       
			r=ae*(r+fcosto)/(r*fcosto+1.0)
			r=ao*(r+fcoste)/(r*fcoste+1.0)
		end do

		ffv=(fv-fo)/(fv+fo)
		r=(r+ffv*fnevot_v)/(r*ffv*fnevot_v+1.0)
		
		pp = Dimag(r)
		qq = Dreal(r)
		ref(j)=PP**2+QQ**2

	ENDDO

!	open(1,file="D:\workOpt\simulatore\sonounadll.txt",position='append')
!		write(1,*) "scritto da reflex"
!		write(1,*) "-----------------"
!		write(1,*) "nener ener reflex"
!		do i=	1,nener
!			write(1,*) i,ener(i),ref(i)
!		end do

!		write(1,*) "-----------------"
!		write(1,*) "nbilayer de do (from substrate)"
!		do i=	1,n_bilayers			
!			write(1,*) i,de(i),do(i)
!		end do
!		write(1,*) "-----------------"
!		write(1,*) "angle=",angle
!	close(1)

	RETURN
END	subroutine


end module 
