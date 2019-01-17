program  createGeo
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


	!2012/03/01 Program created for Web interface.
	!Create a geometry file (shellStruct.dat) starting from a geometrical
	!	parameters file (geoPars.txt).

	! coefficient and profiles are calculated in external routines.
	! Profile inversion (from max diameter to mid diameter) is tested and works
	! for double cone and wolter.
	!The parametrization of exit pupil follows conconi campana paper. 
	!N.B.: the inversion is based on a second order approximation
	!	(Wolter profile) even for different profiles. It gives slightly inaccurate
	!	results in that case.
	!N.B.: If a linear relation is used for thickness, the shell thickness is
	!	calculated as a function of the maximum diameter.  

	!Dimensions of central block are calculated as a further shell, but values are not saved 
	!'profileType' in namelist 'creatediams' for the case of different 
	!	profiles for the two suefaces
	!	2 letters string: 
	! c->cono, p->parabola, h->hyperbole, u->user-defined (coefficients of polynomial profile)

	!cosa succede nel caso di una sola shell?

	!26-11-2008 baltimora: il meccanismo con linearthick non e' sicuro se si 
	!sbagliano le impostazioni (comunque ora non mi ricordo come funziona,
	!quindi controllare istruzioni su quaderno). minthickmm non funziona
	!(vedi caso dithering simbolx65). 

	implicit none

	real(8),parameter::H_C = 12.39841857D0 
	real(8),parameter::TORAD = 0.017453292519943295769237d0
	real(8), parameter ::PI=3.141592653589793238462643D0 !pi greek
	integer,parameter::	nMaxSh=1000
	logical (4) result
	character(100) diamFile

	integer nw,nshdaimp
	logical diamDecresce
	real(8) F_LENGTH,F_HEIGHT1_mm,F_HEIGHT2_mm,f_height
	!f_length e' in cm
	!numero di shell,nome del file dei diametri da passare a geoinp
	!flag per la creazione da file dei diametri,diamdecresce=t se il file dei diametri e' decrescente(calcolato in geoinp)
	
	integer i,io
	real(8) F_lengthDaImp_m,f_heightDaImp_cm
	real(8) F_HEIGHT1_cm,F_HEIGHT2_cm  !HHH
	integer k
	real(8)r,alkw,spes,r0,r1
	real(8)peso,area,thick
	real(8)diam_max,rmax1,diam_min,rmin1,diam_med,rmed1,theta
	real(8)rstart,weight,acoll,temp(nMaxSh+1)

	dimension rmin1(nMaxSh+1),rmed1(nMaxSh+1),rmax1(nMaxSh+1),acoll(nMaxSh+1)
	dimension thick(nMaxSh+1),theta(nMaxSh+1),weight(nMaxSh+1)

	real(8) a1(nMaxSh+1),a2(nMaxSh+1),a3(nMaxSh+1),a4(nMaxSh+1),a5(nMaxSh+1)
	real(8) b1(nMaxSh+1),b2(nMaxSh+1),b3(nMaxSh+1),b4(nMaxSh+1),b5(nMaxSh+1)
	real(8) tpeso,tarea
	real(8)	FieldOfViewDeg,ThickExtSh,LinearThick,WallDensity
	real(8) massimo
	real(8) ThickM,ThickQ,minThickmm,F_LENGTHMM
	character(2) profileType

	real(8) rtoll !tollerance in iterative procedure
	real(8) apot_p,apot_h,volume
	real(8), external :: polyProfile  !,polyInvert 
    character(8)  :: date
    character(10) :: time
    character(5)  :: zone
	integer :: nargs

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
	!      questa versione non usa al, ma considera al= f_height
	!      per la versione con le due variabili separate vedere geo_inpAl
	!-----
	!il programma calcola solo l'area geometrica (collecting area)
	!-----

	namelist /GeoPars/ nshDaImp,F_lengthDaImp_m,f_heightDaImp_cm,f_height1_cm, &
	f_height2_cm

!------------------------------------------------------------------------------------
	print*,"Starting fortran program creaGeo on:"
	CALL DATE_AND_TIME(DATE, TIME, ZONE)
	print*,date,time,zone 
	!default values
	linearThick=2																																											 ! 
	massimo=0	!mancata lettura
	io=0	  
	temp=0
	rmed1=0
	rtoll=1d-10
	diamDecresce=.true.	
	minThickmm=0
	nshDaImp=0
	
	!The following group will be not compiled by MS compiled.
	!MS$IF .false.
	nargs = command_argument_count()
	if (nargs /=0) call get_command_argument(1, diamFile) 	
	!!nargs = iargc()
	!!if (len(trim(diamFile))/=0) nargs=1
	!!call get_command_argument (1, diamFile)
	!!call getarg (1, diamFile)
	!MS$ENDIF

	open (10,file='geoPars.txt',status='old',iostat=io)
	if (io/=0) then 
		print *,"Error in opening file geoPars.txt"
		print*,"iostat= ",io
		stop
	endif
	read (10,CreateDiams)
	read(10,geoPars)
	close(10)

	if (f_heightdaimp_cm==0) then
		if (f_height1_cm * f_height2_cm /= 0) then	!usa f_height1_cm e 2
			f_height1_mm= f_height1_cm*10							
			f_height2_mm= f_height2_cm*10
			f_height=0
		else
			f_height1_mm=f_height*10				!usa quello eventualmente letto se usestartdirflag=1
			f_height1_mm=f_height*10				!senno' resta 0, da' errore successivamente	(?)
		end if
	else
		if (f_height1_mm/=0 .or. f_height2_mm/=0 ) then		!impostazioni non chiare, da' errore
			print *, "ERRORE sia F_heightDaImp_cm che F_height(1|2)_cm "
			print *, "sono impostati diversi da 0:"
			print *, "f_heightdaimp_cm: ", f_heightdaimp_cm
			print *, "f_height1_cm: ", f_height1_cm
			print *, "f_height2_cm: ", f_height2_cm
			print *, "imposta a zero o togli dalla namelist quelli che non interessano"
			stop
		else
			f_height1_mm=f_heightdaimp_cm*10					!usa f_heightdaimp_cm
			f_height2_mm=f_heightdaimp_cm*10
		end if
	end if
	!HHH fine
	if (f_lengthdaimp_m/=0) f_length=f_lengthdaimp_m*100
	print*,"F_length (cm) = ",F_length	!HHH
	print*,"f_height1 (mm) = ",f_height1_mm	!HHH
	print*,"f_height2 (mm) = ",f_height2_mm	!HHH

	!fissa nw (numero di shell)	da passare a geoinp
	if  (nshDaImp/=0)	then
		if (nshDaImp > nMaxSh) then
			print*,"The maximum number of shells is set to ",nmaxsh,","
			print*,"modify the program."  
			stop
		endif
		nw=nshDaImp
		print*,"Shell number:",	nw
	else
		print*,"ERROR: namelist geopars not found in geoPars.txt"
		stop
	end if

	if(nargs/=0) then  !create geometry file from list of diameter	
		print*, "Read diameters from file ", diamFile
		open (1,file=diamFile,status='old',iostat=io)
		if (io/=0) then
			print*,"Error nr. ",io," in reading diameters file",diamFile
			stop
		end if
		do i =1,nw
			read(1,*,iostat=io) rmed1(i)
			if (io /=0) then
				print*,"Error nr. ",io," in reading diameters from ",diamFile		
				print*,"at shell  nr. ",nw
				print*,"last diameter read (",i,")= ",rmed1(i)
				print*,"previous= ",rmed1(i-1)
				print*,"expected number of shells: ",nw		 
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
		rstart=maxval(rmed1)
		print*,"Maximum diameter (at IP) as read from starting diameter file:", rstart*2		 	
	else if(massimo/=0)then
			print*,"Maximum diameter as read from parameters file: ", massimo
			rstart=massimo/2
	else
		print*,"Maximum diameter not found in namelist CreateDiams in impOffaxis.txt,"
		stop
	end if
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
			print*,"set to 0 or not present."
			print*,"ThickExtSh= ",thickExtSh
			print*,"ThickQ= ",thickQ
			stop
		end if
		if (thickExtSh/=0) thickQ=thickExtSh-2*thickM*rstart !calcola con thickExtSh
		!senno' usa i thickQ e thickM impostati
	endif

	F_LENGTHmm=F_LENGTH*10
	r=rstart
	do i=1,nw+1 !ciclo do per costruzione degli specchi	e calcolo del peso
	!gliene faccio calcolare uno in piu' per usarlo per blocco centrale
	!*******************************************************		
		
		if (nargs==0 .or. i==nw+1 ) then
			!il raggio usato per il calcolo degli spessori e' quello esterno della
			!pupilla d'uscita, l'ultima shell fittizia ha spessore 0.
			if (i== nw+1) then
				spes=0
			else
				spes = (ThickQ+ThickM*2*r)/(1-thickM*2)		 !il denominatore si rende necessario
				spes= max(spes,minThickmm)
			end if											!in quanto r e' diametro esterno (comprende spessore)	
			
			!calculate the sequence of rmed1.
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
	
	print*,'peso totale?',' ',tpeso
	print*,'area totale',' ',tarea

	!	scrittura dei parametri geometrici della shell
	open (112,file="shellStruct.dat")
	area=0
	write (112,'(a100)') "Nshell   Dmax(mm)   Dmid(mm)   Dmin(mm)   thickness(mm)   Angle(rad)   Area(cm^2)  Mass(kg)"
	11 format  (i4,7f15.6)
	do k=1,nw
		write (112,11) k, rmax1(k)*2,rmed1(k)*2,rmin1(k)*2,thick(k),theta(k), acoll(k),weight(k)
	end do
	close(112)

	print*,'Geometry file successfully created on:'
	CALL DATE_AND_TIME(DATE, TIME, ZONE)
	print*,date,time,zone 

	print*
end program

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
		stop			
	else
		print*, "profile type for the first surface (",profileType(1:1),")" 
		print*, "not recognized (valid: c,h,p)"
		stop
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
		b2=2*(radiusIP*dtan(angle_two))/(F_lengthmm+radiusIP/tan(2*angle_one))
		!b2=2*(radiusIP*dtan(angle_two))/(F_lengthmm+radiusIP*dcotan(2*angle_one))
	else
			print*, "profile type for the second surface (",profileType(2:2),")" 
			print*, "not recognized (valid: c,h,p)"
			stop
	end if

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





