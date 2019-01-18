program  createGeo
!230209 geoinp 6: routine per calcolo coefficienti polinomio a e b
!per calcolare accuratamente la geometria con un generico profilo.
!spostati i calcoli di coefficienti e profili in routine esterne.
!L'inversione dei profili funziona correttamente per i profili implementati 
!(doppio cono e Wolter)
!#230209# cambiato segno coefficiente a2 in tipo di profilo. inserito rtoll per
!tolleranza in procedura iterativa. Aggiunti a1 e b1.	Rivisto calcolo volume.
!corretto calcolo degli spessori per prima shell (teneva conto di diametro 
!parete esterna per la prima shell)
!inserito minthickness, prima non funzionava (era mancante)
!OSSERVAZIONE: il metodo ricorsivo per la generazione delle shell si basa
!su wolter anche quando altri profili sono assunti. da' risultato impreciso in tal caso.
!in piu' lo spessore viene calcolato sul diametro massimo. 
!26-11-2008 baltimora: il meccanismo con linearthick non e' sicuro se si 
!sbagliano le impostazioni (comunque ora non mi ricordo come funziona,
!quindi controllare istruzioni su quaderno). minthickmm non funziona
!(vedi caso dithering simbolx65). 
!14-10-2008 changed the parametrization of exit pupil, according to conconi
!campana paper. changed setting	of createfileshells here and in subplot
!to accept diameter file.
!12-9-2008 corretta impostazione spessori tramite thickExtsh che dava
!spessore doppio
!21-1-08 cambridge: geo_inp_3,changed calculation of central block,
!now calculated as a further shell
!14-11-07 cambridge: added 'profileType' in namelist 'creatediams'  
!for the case of different profiles for the two suefaces
!2 letters string: c->cono, p->parabola, h->hyperbole
!6-11-07 aggiunto supporto per 2 semialtezze diverse, da OA3
!13-4-07 traie4 tentativo di fare veramente possibilita' una shell

!22-12-06 da palermo
!rispetto a geo_inp cambia due cazzate nel caso di una shell
!ancora incompleto, modificare originale per gestire caso 1 shell

!9-3-2006
!aggiunto scrittura pesi alla fine

!prepara il file col calcolo della geometria se non gia' pronto
!nella routine chiamante cretefileshells e namefileshells vengono impostati in
!modo che:
!se createfileshells
!=0 questa routine non viene chiamata,
!		il programma non fa niente (file gia' esistente, nome in namelist impraytrace)
!=1 crea da impostazioni	(ma leggendo diam max da file dei diametri),
!		usa come nome transRayTrace.dat
!=2 crea da file dei diametri	(25-5-05 non implementato)

!----------------

!carica i valori delle variabili puo' essere diverso per ogni dirfom)

!Programma di Conconi modificato da G. Pareschi per calcolare la geometria 
!di un telescopio X a doppio cono (!!!!!)
!NB NBN NB NB: questa versione fornisce l'input per il ray-tracing
!--
! Sistema di riferimento: tre assi cartesiani con origine nel centro
! della circonferenza intersezione fra la front e la rear surface;
! asse F_HEIGHT=asse ottico  positivo verso la rear surface.
! Definizione della superficie della shell tramite:
!   
! Rf*Rf=R0*R0*(1.+2*tan(theta)*zn+a2*zn**2+a3*zn**3+a4*zn**4
!
! Rr*Rr=R1*R1*(1.+2*tan(thetb)*zw+b2*zw**2+b3*zw**3+b4*zw**4
!
! con Rf= raggio della front surface alla coordinata zn
!     Rr= raggio della rear surface alla coordinata zw
!     R0= raggio della front surface a z=0.
!     R1= raggio della rear surface  a Z=0.
!     zn= z/R0 con z coordinata cartesiana lungo l'asse ottico
!     zw= z/R1 con z coordinata cartesiana lungo l'asse ottico
!     a2,a3,a4= coefficienti del polinomio per la front surface
!     b2,b3,b4= coefficienti del polinomio per la rear surface
!     theta= angolo della front surface
!     thetb= angolo della rear surface  in teoria thetb=3*theta
!
!     Focale nominale all'infinito = R0/tan(4.*theta)
!     Focale nominale a distanza d = R0/tan(4.*theta-arctan(R0/d))
!
!      questo programma calcola i valori di rmax,rmed,rmin
!      di un sistema Wolter 1 nota la focale f_length il valore misurato
!      di r a distanza z sulla parabola e la lunghezza proiettata
!      F_HEIGHT assunta  eguale per la iperbole e la parabola.
!
!      Fornisce poi l'equazione della parabola e della iperbole
!      con la rappresentazione ed il sistema di riferimento 
!      descritto in precedenza
!
!	implicit real*8 (a-h,o-z)


implicit none

integer, parameter:: pmaxEn=150	!massimo numero di punti per calcolo energie
integer, parameter:: maxThroughputEn=10	!massimo numero di energie per scrittura throughput
integer, parameter:: nMaxSh=400 !massimo numero di shell
real(8),parameter::H_C = 12.39841857D0 
real(8),parameter::TORAD = 0.017453292519943295769237d0
real(8), parameter ::PI=3.141592653589793238462643D0 !pi greek


!---- namelist ---------!
character(40) dirrisult,dirorigin,workdir
character(120) psfdir
namelist /impdir/ workdir,dirorigin,dirrisult,psfdir
!-----
integer numeropar,controllaener
namelist /somma/ numeropar,controllaener ! controllaener non serve
!-----------------------!

integer io,io2
character(40) afflst
logical (4) result
character(100) rec
integer aeff_flag,iresult

open (10,file='geopars.txt',status='old',iostat=io)
if (io/=0) then 
	print *,"Error in opening file imp_OffAxis.txt"
	print*,"iostat= ",io
	pause ""
	stop
endif
read (10,impdir)
read(10,somma)
close(10)

afflst="af_files"
affl=trim(workdir)//"\"//trim(afflst)

io=0
result = MAKEDIRQQ(trim(workdir))
!print*,"copy "//trim(afflst)//" "//trim(workdir)//"\"//trim(afflst)//"\"
!pause ""

!-- Common Blocks
!------------
!character (20) df
integer nw
logical diamDecresce
real(8) F_LENGTH,F_HEIGHT1_mm,F_HEIGHT2_mm
common /to_geoinp/ nw,F_Length,F_Height1_mm,F_Height2_mm,diamdecresce !HHH 
!f_length e' in cm
!numero di shell,nome del file dei diametri da passare a geoinp
!flag per la creazione da file dei diametri,diamdecresce=t se il file dei diametri e' decrescente(calcolato in geoinp)
!------------
character(40) affl
common /affilesloc/ affl
!--------
character(100) namefileshells !indica se creare da file dei diametri piuttosto che da impostazioni
common /transfile/ namefileshells
!---
!logical useStartDirFlag
!common /csdf/ useStartDirFlag
!-----

character(*)  dirfom !dirfomor e' quella di origine dirfom e' dei risultati
integer i,io,iresult

integer k,createFileShells 
!real (8) angle_two

real(8)r,alkw,spes
real(8)r0,r1
real(8)peso,area,thick
real(8)diam_max,rmax1,diam_min,rmin1,diam_med,rmed1,theta
real(8)rstart,weight,acoll,temp(nMaxSh)

dimension rmin1(nw+1),rmed1(nw+1),rmax1(nw+1),acoll(nw+1)
dimension thick(nw+1),theta(nw+1),weight(nw+1)

real(8) a1(nw+1),a2(nw+1),a3(nw+1),a4(nw+1),a5(nw+1),b1(nw+1),b2(nw+1),b3(nw+1),b4(nw+1),b5(nw+1)
real(8) tpeso,tarea
real(8)	FieldOfViewDeg,ThickExtSh,LinearThick,WallDensity
real(8) massimo
real(8) ThickM,ThickQ,minThickmm,F_LENGTHMM
character(2) profileType

real(8) rtoll !tollerance in iterative procedure	 !#230209#
real(8) apot_p,apot_h,volume !#230209#
real(8), external :: polyProfile  !,polyInvert 

interface 
	function polyInvert(z,r1,c1,c2,c3,c4,c5)	
		real(8), intent(in):: z,r1,c1,c2 
		real(8), intent(in),optional :: c3,c4,c5
		real(8) polyInvert
	end function 
end interface

!----- Namelists			
namelist /CreateDiams/ FieldOfViewDeg,ThickExtSh,LinearThick, &
& ThickM,ThickQ,minThickmm,WallDensity,profileType,massimo
!
!      intodurre f_length,F_HEIGHT,r,al definite come sopra
!      Es: f_length=3500. F_HEIGHT=300. r=138.188 al=300
!      dati relativi al mandrel n 3 di jet-x
!      questa versione non usa al, ma considera al= f_height
!      per la versione con le due variabili separate vedere geo_inpAl
!-----
!il programma calcola solo l'area geometrica (collecting area)
!-----

! questo pezzo era in subplot, spostato
!0 non fa niente (file gia' esistente,nome in namelist impraytrace)
!1 crea da impostazioni	(ma leggendo diam max da file dei diametri)
!2 crea da file dei diametri

temp=0
rmed1=0
rtoll=1d-10
diamDecresce=.true.	
minThickmm=0

if(CreateFileShells==0) then
	iresult =system("copy "//trim(affl)//"\"//trim(Namefileshells)//" "//trim(dirfom)//"\"//trim(Namefileshells))
	namefileshells=dirfom//"\"//trim(namefileshells)
	return		
else if(CreateFileShells==2) then		
	print*, "CreateFileShells=2,"
	print*, "Read diameters from file ", trim(affl),"\"//trim(Namefileshells)
	open (1,file=trim(affl)//"\"//trim(Namefileshells),status='old',iostat=io)
	if (io/=0) then
		print*,"Error nr. ",io," in reading diameters file",trim(affl),"\"//trim(Namefileshells)
		pause ""	
		stop
	end if
	do i =1,nw
		read(1,*,iostat=io) rmed1(i)
		if (io /=0) then
			print*,"Error nr. ",io," in reading diameters from ",	trim(affl),"\"//Namefileshells		
			print*,"at shell  nr. ",nw
			print*,"last diameter read (",i,")= ",rmed1(i)
			print*,"previous= ",rmed1(i-1)
			print*,"expected number of shells: ",nw		 
			pause ""
			stop
		end if
	end do
	close(1)
	rmed1=rmed1/2
	if (rmed1(1) < rmed1(nw)) then
		diamdecresce=.false.
		temp=rmed1
		do i =1,nw
			rmed1(i)=temp(nw+1-i)
		end do
	end if	
	Namefileshells=dirfom//"\transRayTrace.dat"		 	
else if (CreateFileShells==1) then
	Namefileshells=dirfom//"\transRayTrace.dat"
else 
  !if CreateFileShells==1 the structure file will be created 	
	print *,"The value for CreateFileShells is not valid =",CreateFileShells
	pause ""
	stop
end if

open(25,file=trim(dirfom)//"\transRayTrace.dat", status='new',iostat=I)
IF (I.GT.0) THEN
	PRINT*,"il file ",trim(dirfom)//"\transRayTrace.dat"," esiste gia':ERRORE!"
	PAUSE ''
	stop
END IF

open (10,file='imp_OffAxis.txt')
!default values
linearThick=2																																											 ! 
massimo=0	!mancata lettura
read (10,CreateDiams,iostat=io)
if (io/=0) then
	print*
	print*,"It looks like the namelist CreateDiams does not exists, "
	print*, "or has an error: default values assigned:"
	print*,"FieldOfViewDeg=0.15,ThickM=0.0016,ThickQ=-0.02,minThickmm=0.2,WallDensity=8.8"
	print*,"profile= double cone"
	pause "press a key to go on."
	FieldOfViewDeg=0.15
	linearThick= 2
	ThickExtSh=0
	ThickM=0.0016d0				!mm
	ThickQ=-0.02d0
	minThickmm=0.2d0
	WallDensity=8.8d0					!g/cm**3
	profileType='cc'
end if
if (CreateFileShells==2) then
	rstart=maxval(rmed1)
	print*,"Maximum diameter (at IP) as read from starting diameter file:", rstart*2
else if (CreateFileShells==1) then
	if(massimo/=0 .and. io==0)then
		print*,"Maximum diameter as read from setting file: ", massimo
		rstart=massimo/2
	else
		print*,"Maximum diameter not found in namelist CreateDiams in impOffaxis.txt,"
		print*,"insert the maximum diameter in mm"
		read(*,*) rstart
		rstart=rstart/2
	end if
else 
  !if CreateFileShells=0 it should be returned on the first condition test,
  !it should never happen
	print*,"The value for CreateFileShells is not valid (2) =",CreateFileShells
	pause ""
	stop
end if
close (10)

print*
print*,"-----------------------------"
print*,'Max aperture radius[mm]= ',rstart			


tpeso=0. !inizializza peso
tarea=0. !inizializza area

alkw=F_HEIGHT1_mm*tan(fieldOfViewDeg*PI/180.) !distanza tra uno specchio e l'altro per
! evitare obreggiamenti fuori asse

if (LinearThick==0) then		!cambiato	prima era 1
	thickQ=ThickExtSh
	thickM=0  
else  if (LinearThick==1) then 	!cambiato prima era 2
	thickM=ThickExtSh/(2*rstart)
	thickQ=0
else  if (LinearThick==2) then 	!cambiato prima era 2
	if (thickExtSh*ThickQ/=0) then
		print*,"either thickExtSh or thickQ, must be"
		print*,"set to 0"
		print*,"ThickExtSh= ",thickExtSh
		print*,"ThickQ= ",thickQ
		print*,"press return to go on"
		print*,"(using thickExtSh)"
		pause ""
	end if
	if (thickExtSh/=0) thickQ=thickExtSh-2*thickM*rstart !calcola con thickExtSh
	!senno' usa i thickQ e thickM impostati
endif

F_LENGTHmm=F_LENGTH*10
r=rstart
do i=1,nw+1 !ciclo do per costruzione degli specchi	e calcolo del peso
!gliene faccio calcolare uno in piu' per usarlo per blocco centrale
!*******************************************************		
	
	!il raggio usato per il calcolo degli spessori e' quello esterno della
	!pupilla d'uscita
	!l'ultima shell fittizia ha spessore 0
	if (i== nw+1) then
		spes=0
	else
		spes = (ThickQ+ThickM*2*r)/(1-thickM*2)		 !il denominatore si rende necessario
		spes= max(spes,minThickmm)
	end if											!in quanto r e' diametro esterno (comprende spessore)	

	if (CreateFileShells==1 .or. i==nw+1 ) then
	!iterate to find the intersection plane radius starting from the
	!entrance pupil radius (external)
		r=r-spes !radius at the mirror surface	

		!recursive procedure for r0 calculation from r at the entrance pupil		
		call polyCoeff(profileType,r,F_lengthmm,a1(i),a2(i),a3(i), &
					a4(i),a5(i),b1(i),b2(i),b3(i),b4(i),b5(i))
		r0=polyInvert(-F_HEIGHT1_mm,r,a1(i),a2(i))
	
		do
			call polyCoeff(profileType,r0,F_lengthmm,a1(i),a2(i),a3(i), &
					a4(i),a5(i),b1(i),b2(i),b3(i),b4(i),b5(i)) 
			r1=polyInvert(-F_HEIGHT1_mm,r,a1(i),a2(i))
			
			if(dabs(r1-r0)<=rtoll) exit			
			r0=r1
		end do
		rmed1(i)=r0
	end if
	r=rmed1(i) - alkw !external radius at exit pupil, used in next cycle
	
	theta(i)=datan2(rmed1(i),F_LENGTHmm)/4 !refinement. Assert theta(i)=th (se createFileShells=1)
  call polyCoeff(profileType,rmed1(i),F_lengthmm,a1(i),a2(i),a3(i), &
					a4(i),a5(i),b1(i),b2(i),b3(i),b4(i),b5(i))

	!calculation of the accurate value (2nd order) of radii 
	!at entrance and exit pupils corresponding to rmed
	!assert non dovrebbe cambiare se la procedura iterativa e' accurata.	
	!servono se calcola geometria a partire da diametri medi (createFileShells=2)
	!l'asse z e' orientato verso il piano focale, H e' negativo sulla parabola.
	!non e' necessario cambiare il segno di a2 e b2 perche' al quadrato

	rmax1(i)=polyProfile(-F_HEIGHT1_mm,rmed1(i),a1(i),a2(i),a3(i),a4(i),a5(i))		
	rmin1(i)=polyProfile(F_HEIGHT2_mm,rmed1(i),b1(i),b2(i),b3(i),b4(i),b5(i))		

	!la formula calcola lo spessore in base al diametro alla pupilla di ingresso
	!il denominatore (1-2m) serve perche' si calcola in base al raggio esterno
	!mentre la linearita' e' con il raggio interno.	Anche questo non dovrebbe cambiare
	!if (CreateFileShells==2) thick(i)=(ThickQ+ThickM*2*rmax1(i))/(1-ThickM*2)
	spes=(ThickQ+ThickM*2*rmax1(i))/(1-ThickM*2)
	spes=max(spes,minThickmm)
	thick(i)=spes

	!check for shell not overlapping at the exit pupil 
    if (i/=1)then
		if (rmin1(i)+spes > rmin1(i-1)) then
			print *, "Shell overlapping (at the exit pupil) for shell nr.",i
			print *,"exit pupil radii:"
			print *,"shell nr.",i-1," inner: ",rmin1(i-1)
			print *,"shell ",i," outer: ",rmin1(i)+spes
			pause ""
			stop
		end if
 	end if 

	!peso=(rmax+rmed+rmin)*6.283*spes*2*F_HEIGHTmm*WallDensity/3000000.
	!la superficie laterale del tronco di cono e' pi*(r+R)*a, dove a e'
	!l'apotema
	if (i <= nw) then
		!peso= pi*((rmax1(i)+rmed1(i))/dcos(theta(i))*F_HEIGHT1_mm+&
		!			(rmed1(i)+rmin1(i))/dcos(3*theta(i))*F_HEIGHT2_mm)&
		!			*thick(i)*WallDensity/1000000.
		apot_p=F_HEIGHT1_mm/dcos(theta(i))
		apot_h=F_HEIGHT2_mm/dcos(3*theta(i))
		volume=spes*pi*((rmax1(i)+rmed1(i))*apot_p+(rmed1(i)+rmin1(i))*apot_h)
		peso=volume*WallDensity/1000000.			
		weight(i)=peso
		area=pi*(rmax1(i)**2-rmed1(i)**2)/100 ! collecting area
		acoll(i)=area
		tarea=area+tarea
		tpeso=tpeso+peso
	end if
	
!*******************************************************
end do     ! termine ciclo do per calcolo geometria

!scrittura su file 25
	write(25,51) 
 51	FORMAT(2x,'wall thickness [mm]')
	do I=1,nw
	   write(25,*) thick(i)
	enddo
!---
	write(25,56)
 56	FORMAT(2x,'max diameters of the entrance pupil [mm]')
	do I=1,nw
	   diam_max=rmax1(i)*2.
	   write(25,*) diam_max
	enddo
!---
	write(25,57)
 57	FORMAT(2x,'min diameters of the exit pupil [mm] for next inner shell')
	
	!if (nw>1)	then
	do I=2,nw+1
		 diam_min=rmin1(i)*2.
		 write(25,*) diam_min
	enddo
	!	diam_min=(2*rmin1(nw)) - (2*(rmin1(1)-rmin1(nw)))/float(nw-1)
	!else
	!	dr_perno=alkw+F_HEIGHT1_mm*dtan(theta(1))+thick(1) 
	!	diam_min=	2*rmin1(nw)	- 2*dr_perno  !!16/1/2008
	!end if
	!write(25,*) diam_min

	write(25,58)
 58	FORMAT(2x,'min diameters of the intermediate pupil [mm] for next inner shell')
	do I=2,nw+1
	   diam_med=rmed1(i)*2.
	   write(25,*) diam_med
	enddo
	!if (nw>1)	then
	!	diam_med=(2.*rmed1(nw))-(2.*(rmed1(1)-rmed1(nw)))/float(nw-1)
	!else 
	!	diam_med=2*rmed1(nw)	- 2*2*dr_perno
	!end if
	!write(25,*) diam_med
	open(26,file=dirfom//"\diam.dat")
	!write(26,*) diam_med corretta 9/10/2006 prima partiva da shell n+1 e saltava la prima
	do I=nw,1,-1	 !do I=nw,2,-1
		 write(26,*) rmed1(i)*2
	enddo
	close(26)

	write(25,74)
 74	FORMAT(2x,'radius at the intermediate pupil [mm]')
	do I=1,nw
	   write(25,*) rmed1(i)
	enddo

	write(25,59)
 59	FORMAT(2x,'radius at the intermediate pupil bis [mm]')
	do I=1,nw
	   write(25,*) rmed1(i)
	enddo

	write(25,60)
 60	FORMAT(2x,'theta angle first cone in rads')
	do I=1,nw
	   write(25,*) -theta(i)
	enddo

	write(25,61)
 61	FORMAT(2x,'theta angle second cone in rads')
	do I=1,nw
	   write(25,*) -(3*theta(i))
	enddo
	
	write(25,63)
 63	FORMAT(2x,'Valori di a2,a3,a4,a5')
	do I=1,nw		
	  write(25,*) a2(i),a3(i),a4(i),a5(i)
	enddo

	write(25,64)
  64	FORMAT(2x,'Valori di b2,b3,b4,b5')

	do I=1,nw
	   write(25,*) b2(i),b3(i),b4(i),b5(i)
	enddo
	close(25)
	
	print*,'peso totale?',' ',tpeso
	print*,'area totale',' ',tarea


open (111,file=trim(dirfom)//"\pesi.txt")
write (111,*) 'peso totale',' ',tpeso
write (111,*) 'area totale',' ',tarea
write (111,*)
write (111,*)

do i=1,nw
write (111,*)	weight(i)
end do

close (111)

!	scrittura dei parametri geometrici della shell
open (112,file=trim(dirfom)//"\shellStruct.txt")
area=0
write (112,'(a100)') "Nshell   Dmax(mm)   Dmid(mm)   Dmin(mm)   thickness(mm)   Angle(rad)   Area(cm^2)  Mass(kg)"
11 format  (i4,7f15.6)
do k=1,nw
	write (112,11) k, rmax1(k)*2,rmed1(k)*2,rmin1(k)*2,thick(k),theta(k), acoll(k),weight(k)
end do
close(112)

end subroutine

subroutine polyCoeff(profileType,radiusIP,F_LENGTHmm,a1,a2,a3,a4,a5,b1,b2,b3,b4,b5)	
implicit none	 
character(2), intent (IN):: profileType
real(8), intent (IN):: radiusIP,F_LENGTHmm
real(8), intent (OUT):: a1,a2,a3,a4,a5,b1,b2,b3,b4,b5
real(8) angle_one,angle_two

	!set a and b cofficients according to the profile type
	!attenzione al segno di theta, che e' positivo, ma viene scritto sul file
	!cambiato di segno, a2 invece e' negativo in caso di doppiocono
	angle_one=datan2(radiusIP,F_LENGTHmm)/4
	a1=-2*dtan(angle_one)
	a2=0
	a3=0.
	a4=0. 
	a5=0.
	if (profileType(1:1)=='c' .or. profileType(1:1)=='C') then
		a2=dtan(angle_one)**2 !-(dtan(theta(i)))**2		!#230209#
	else if (profileType(1:1)=='p' .or. profileType(1:1)=='P') then
		a2=0
	else if (profileType(1:1)=='u' .or. profileType(1:1)=='U') then
		!polynomial profile per ora li metto a mano qui
		a2=0.0002005d0	
		a3=0.0003318d0
		a4=	0.0002711d0	
		a5=0.00009652d0	
	else if (profileType(1:1)=='h' .or. profileType(1:1)=='H') then
		print*,"ERROR: THE hyperbolic surface for the primary mirror"
		print*, "is not implemented yet.."
		pause ""
		stop			
	else
		print*, "profile type for the first surface (",profileType(1:1),")" 
		print*, "not recognized (valid: c,h,p)"
		stop
		pause ""
	end if

	angle_two=3*angle_one
	b1=-2*dtan(angle_two)
	b2=0
	b3=0.
	b4=0. 
	b5=0.
	if (profileType(2:2)=='c' .or. profileType(2:2)=='C') then
		b2=dtan(angle_two)**2
	else if (profileType(2:2)=='p' .or. profileType(2:2)=='P') then
		b2=0
	else if (profileType(1:1)=='u' .or. profileType(1:1)=='U') then
		!polynomial profile per ora li metto a mano qui
		b2=0.001d0	
		b3=0.0004591d0
		b4=	-0.0004248d0	
		b5=0.0001651d0	
	else if (profileType(2:2)=='h' .or. profileType(2:2)=='H') then
		b2=2*(radiusIP*dtan(angle_two))/(F_lengthmm+radiusIP*dcotan(2*angle_one))
	else
			print*, "profile type for the second surface (",profileType(2:2),")" 
			print*, "not recognized (valid: c,h,p)"
			stop
			pause ""
	end if


!	!a2 and b cofficients according to the profile type
!	!attenzione al segno di theta, che e' positivo, ma viene scritto sul file
!	!cambiato di segno
!	if (profileType(1:1)=='c' .or. profileType(1:1)=='C') then
!		a2(i)=(dtan(theta(i)))**2 !-(dtan(theta(i)))**2		!#230209#
!	else if (profileType(1:1)=='p' .or. profileType(1:1)=='P') then
!		a2(i)=0
!	else if (profileType(1:1)=='h' .or. profileType(1:1)=='H') then
!		print*,"ERROR: THE hyperbolic surface for the primary mirror"
!		print*, "is still not implemented.."
!		pause ""
!		stop			
!	else
!		print*, "profile type for the first surface (",profileType(1:1),")" 
!		print*, "not recognized (valid: c,h,p)"
!		stop
!		pause ""
!	end if
!	a1(i)=-2*dtan(theta(i))
!	a3(i)=0.
!	a4(i)=0. 
!	a5(i)=0.
!	angle_two=3*theta(i)
!	if (profileType(2:2)=='c' .or. profileType(2:2)=='C') then
!		b2(i)=(dtan(angle_two))**2
!	else if (profileType(2:2)=='p' .or. profileType(2:2)=='P') then
!		b2(i)=0
!	else if (profileType(2:2)=='h' .or. profileType(2:2)=='H') then
!		b2(i)=2*(Rmed1(i)*dtan(angle_two))/(F_lengthmm+Rmed1(i)*dcotan(2*theta(i)))
!	else
!			print*, "profile type for the second surface (",profileType(2:2),")" 
!			print*, "not recognized (valid: c,h,p)"
!			stop
!			pause ""
!	end if
!	b1(i)=-2*dtan(3*theta(i))
!	b3(i)=0.
!	b4(i)=0. 
!	b5(i)=0.

end subroutine

function polyInvert(z,r1,c1,c2,c3,c4,c5)	
implicit none
real(8), intent(in):: z,r1,c1,c2 
real(8), intent(in),optional :: c3,c4,c5
real(8) polyInvert 
real(8) deltaroot
!return the radius at intersection plane, given the five coefficients of the polynomial profile
!(due per ora) and a value r1 of the radius at a given z.
deltaroot=dsqrt((c1*z)**2-4*c2*z**2+4*R1**2)
polyInvert=(-c1*z+deltaroot)/2	  
!se il delta non e' uguale a zero (doppio cono), prendo la soluzione col +,
!l'altra sara' negativa, infatti Rz(R0)=R0**2+a1*R0*z+a2*z**2
!e' parabola rivolta verso l'alto con vertice con x negativo -a1*z/2,
!quindi per ogni Rz ci saranno due possibili soluzioni R0, una >0, l'altra <0.
	 
end function

function polyProfile(z,r0,c1,c2,c3,c4,c5)	
implicit none
real(8), intent(in):: z,r0,c1,c2,c3,c4,c5
real(8) polyProfile,t 
!return the radius for a given z, from the five coefficients of the polynomial profile
t=z/r0
polyProfile=r0*sqrt(1+c1*t+c2*t**2+c3*t**3+c4*t**4+c5*t**5)
	 
end function



end program


