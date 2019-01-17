!MS$ATTRIBUTES DLLEXPORT :: READINDEX
!MS$ATTRIBUTES ALIAS: 'readindex' :: readindex
!MS$ATTRIBUTES DLLEXPORT :: REFLEX
!MS$ATTRIBUTES ALIAS: 'reflex' :: reflex


!dllModule_working6.f90 + testdll_working8.f90 (finale e funzionante)

SUBROUTINE readindex(argc, argv) !Called by IDL
	!readindex(ener)
	use reflexMod
	use portlib
	implicit none
	INTEGER*4 argc, argv(*) 
	integer*4 j
	
	!open(1,file="D:\workOpt\simulatore\sonounadll.txt")
	!write(1,*) fdate()
	!close(1)
	j = LOC(argc) !Obtains the number of arguments (argc)

	CALL readIndexCore(%VAL(argv(1)), %VAL(argv(2)), &
		%VAL(argv(3)), %VAL(argv(4)), %VAL(argv(5)), &
		%VAL(argv(6)), %VAL(argv(7)), %VAL(argv(8)))
		

END


SUBROUTINE REFLEX(argc, argv) !Called by IDL
	! reflexCore(dSpacing,n_bilayers,reflex,nener,angle,rough)
	use reflexMod									
	implicit none

	INTEGER*4 argc, argv(*) !Argc and Argv are integers
	integer*4 j
	
	j = LOC(argc) 

	CALL reflexCore (%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)), %VAL(argv(4)), %VAL(argv(5)), %VAL(argv(6)))

	RETURN
END



!INTEGER*4 argc, argv(*) !Argc and Argv are integers
!i vettori si passano con array e numero di elementi, che sono contenuti in argv
!arv(1) = array, argv(2) = n elementi

!	j = LOC(argc) !Obtains the number of arguments (argc)
!	!Because argc is passed by VALUE.
!
!	! Call subroutine SUM_ARRAY1, converting the IDL parameters
!	! to standard FORTRAN, passed by reference arguments:
!
!	CALL myRoutine (%VAL(argv(1)), %VAL(argv(2)), %VAL(argv(3)))

!il blocco interfaccia non dovrebbe serbire in quanto in modulo
!interface 
!SUBROUTINE myRoutine(vec1,n,vec2)
!	real*4 vec1(n),vec2(n)
!end subroutine
!end interface
!________________________

! vettori di lunghezza varibile non si possono usare in common blocks.
! possibili soluzioni:
!	1) definirli come variabili di lunghezza costante e usare solo la 
! parte fino ad un dato indice nener
!	2) come 1, ma definendo i vettori a livello di modulo.
! 3) usare array allocatable definite come variabili di modulo. l'array
! indici si alloca alla lettura (ma non si sa quando deallocarlo...)
! reflex viene deallocato e riallocato prima del calcolo della riflettivita'
! ma anche qui non si sa come deallocarlo
! 4) puntatori, ma non esageriamo con la fantasia.

! uso la 2 che mi sembra la piu' facile da implementare
