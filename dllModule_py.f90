module reflexMod
implicit none
!dllModule_working6.f90 + testdll_working8.f90  (finale e funzionante)
!TODO: mettere controllo dimensioni e layer dispari, con scrittura su file di 
!log in caso di errore.
!mettere materiali come argomento.
!gestione overcoating-> se layer dispari tiene sempre il materiale leggero
!in superficie.

!2011/11/07: running on (new) remote machine at itc, same error as below.
! N.B.: f2py works in compiling fortran.
! try to solve the problem. It seems the problem is with the return value
! being an array.
! 
! try to transform in a subroutine (copied previous)
! changes marked with !++
! changes:
! -renamed the internal variable reflex (previously function value) 
! to ref, to avoid name conflict between the function and the variable.

!N.B.; 2011/10/30: ubuntu: it works with python, it doesn't work with ipython (segmentation fault).
!works also in ipython if reflexf90 is imported before launching testReflex.pymodule
!(or the corresponding caller program).
!
!after that, not able to make it work with any system,
! in python or in ipython (when able to prevent crashing, sometimes activating pdb in ipython
! does it), it gives error:
!failed to create intent(cache|hide)|optional array-- must have defined
!    dimensions but got (-1,)

!changes 2011/10/27 Vincenzo
!-made all real (8) and complex(8), f2py does not complain
!-TODO: there is some uncertain point in the nevot-croce calculation:
!arg_e,f,s,o were real obtained as a function of complex.
!I changed them to complex (and the corresponding Dexp to Cdexp),
!but it does not solve the problem, now it is on fnevot_e,f,s,o.
!Check the original programs and paper formula.
!- TODO: reflex has dimensions NENER (correct), ener has dimension PMAXEN,
!- ener was dimensioned with ener(pmaxener), changed to ener(nener)
!- reflex defined also as module variable reflex(pamxen), removed.


complex,parameter :: ci=(0.0,1.0)
real(8),parameter::  H_C = 12.39841857D0 
real(8),parameter::  TORAD = 0.017453292519943295769237d0
real(8), parameter ::PI=3.141592653589793238462643D0 !pi greek

contains

!++function reflex(dSpacing,n_layers,ener,nener,angle,rough,rs2,re2,ro2)
subroutine reflex(dSpacing,n_layers,ener,nener,angle,rough,rs2,re2,ro2,ref)
!versione semplificata della routine originaria: una sola polarizzazione
!solo strati alternati con mat1odd sullo strato esterno.
!dSpacing fornito a partire dal substrato.

implicit none
integer*4, intent(in):: nener,n_layers
!++real(8)  reflex(nener)
real(8), intent(out):: ref(nener)
real(8), intent(in):: dSpacing(n_layers)
real(8), intent(in):: ener(nener),angle,rough
integer*4 n_bilayers
integer i,j,t
real(8) de(n_layers/2),do(n_layers/2) !n_layers deve essere pari
complex(8),intent(in):: re2(nener), rs2(nener), ro2(nener)

complex(8) fo,fe,fv,ffe,ffv,ffo    
complex(8) ao,ae,r,fs,ffs
complex(8) FcostO,FcostE,costO,costE
complex(8) ro2j, re2j, rs2j
real(8) SIGMA_O2,SIGMA_E2,SIGMA_S2
real(8) seno,coseno,refv,prefact
real(8) pp,qq
complex(8) arg_e,arg_o,arg_v,arg_s
real(8) fnevot_e,fnevot_o,fnevot_v,fnevot_s
real(8) VAL,r_l
real(8) oddrough,evenrough,subrough

    n_bilayers=int(n_layers/2)
    de=0
    do=0
    DO I=1,N_biLAYERS !i va dal basso            
        DE(I)=Dspacing(i*2-1)
        DO(I)=Dspacing(i*2)
    ENDDO
    
    !open(1,file="sonoUnaDll.txt",position='append')
    !DO I=1,N_biLAYERS !i va dal basso            
    !    write(1,*) DE(I),DO(I)
    !ENDDO
    !close(1)

    oddrough=rough
    evenrough=rough
    subrough=rough
    coseno=cos(angle)
    seno=sin(angle)
    refv=dabs(seno)
    fv = Dcmplx(refv)
    ref=0
    
    !!test--
    !ref=ener
    !return
    !!---
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
        ffo=-ffe                            !  angle and critical angle)
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
        FNEVOT_E=cDEXP(-prefact*arg_e)
        FNEVOT_O=cDEXP(-prefact*arg_o)
        FNEVOT_V=cDEXP(-prefact*arg_v)
        FNEVOT_S=cDEXP(-prefact*arg_s)

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
    
!    open(1,file="D:\workOpt\simulatore\sonounadll.txt",position='append')
!        write(1,*) "scritto da reflex"
!        write(1,*) "-----------------"
!        write(1,*) "nener ener ref"
!        do i=    1,nener
!            write(1,*) i,ener(i),ref(i)
!        end do

!        write(1,*) "-----------------"
!        write(1,*) "nbilayer de do (from substrate)"
!        do i=    1,n_bilayers            
!            write(1,*) i,de(i),do(i)
!        end do
!        write(1,*) "-----------------"
!        write(1,*) "angle=",angle
!    close(1)

    !++RETURN
!++ END FUNCTION
END subroutine


end module 
