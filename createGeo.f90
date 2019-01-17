program  createGeo
	
	
	! create a geometry file

	!2019/01/17 fix area value, rewrite documentation
	!2012/06/25 from v 1.1

	!2012/03/28 from v1.0
	! single namelist.
	! changed variable names.
	! Added (unused) routine LINEARIZE.							   	
	! TODO: implement linear mechanism for 2 shell heights.
	! TODO: add error codes in messages.
	! TODO: add value check
	! In this version thickness must be defined only by m and q.
	! A unique value is accepted for both shell lengths, this is
	!	made to avoid poblems if e.g. the secondary is set to a shorter
	!	value. The length of the secondary indeed does not enter in the
	!	calculation of the area, but only for mass.  
	!	
	!	t=m*D+q
	!
	!2012/03/01 Program created for Web interface.
	!Create a geometry file (shellStruct.dat) starting from a geometrical
	!	parameters file (geoSettings.txt).
	!A commandline argument can be passed, in that case, it is considered a
	!	filename from which to read diameters at intersection plane.
	!	N.B.: this part is excluded by MS FORTRAN, since the functions 
	!	command_argument_count() and get_command_argument() are not recognized. 

	!-Coefficients and profiles are calculated in external routines.
	!-Profile inversion (from max diameter to mid diameter) is tested and works
	! for double cone and wolter.
	!-The parametrization of exit pupil follows conconi campana paper. 
	!N.B.: the inversion is based on a second order approximation
	!	(Wolter profile) even for different profiles. It gives slightly inaccurate
	!	results in that case.
	!N.B.: If a linear relation is used for thickness, the shell thickness is
	!	calculated as a function of the maximum diameter.  

	!Dimensions of central block are calculated as a further shell, but values are not saved 
	!'profileType' for the case of different profiles for the two surfaces
	!	2 letters string: 
	! c->cono, p->parabola, h->hyperbole, u->user-defined (coefficients of polynomial profile)
	!N.B.:2012/03/28 the polynomial profile is not really implemented, it uses fixed
	!	numerical values (add e.g. reading from file)
	!
	!cosa succede nel caso di una sola shell?

	!--
	! The block in italian is from the original documentation. Verify if still valid.
	!
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

	!26-11-2008 baltimora: il meccanismo con linearthick non e' sicuro se si 
	!sbagliano le impostazioni (comunque ora non mi ricordo come funziona,
	!quindi controllare istruzioni su quaderno). minthickness_slopem non funziona
	!(vedi caso dithering simbolx65). 

	implicit none

	real(8),parameter::H_C = 12.39841857D0 
	real(8),parameter::TORAD = 0.017453292519943295769237d0
	real(8), parameter ::PI=3.141592653589793238462643D0 !pi greek
	integer,parameter::	nMaxSh=1000
	!logical (4) result
	character(100) diamFile

	integer nw,nShells
	logical diamDecresce
	real(8) F_LENGTH,shell_length1,shell_length2,shell_length	!,f_height
	!real(8) H1_slope,H1_intercept,H2_slope,H2_intercept,H1_start,H2_start

	integer i,io
	real(8) focal_length	!,f_heightDaImp_cm
	!real(8) F_HEIGHT1_cm,F_HEIGHT2_cm  !HHH
	integer k
	real(8)r,alkw,spes,r0,r1
	real(8)peso,area,thick
	real(8)rmax1,rmin1,rmed1,theta
	real(8)rstart,weight,acoll,temp(nMaxSh+1)

	dimension rmin1(nMaxSh+1),rmed1(nMaxSh+1),rmax1(nMaxSh+1),acoll(nMaxSh+1)
	dimension thick(nMaxSh+1),theta(nMaxSh+1),weight(nMaxSh+1)

	real(8) a1(nMaxSh+1),a2(nMaxSh+1),a3(nMaxSh+1),a4(nMaxSh+1),a5(nMaxSh+1)
	real(8) b1(nMaxSh+1),b2(nMaxSh+1),b3(nMaxSh+1),b4(nMaxSh+1),b5(nMaxSh+1)
	real(8) tpeso,tarea
	real(8)	angular_shell_separation_deg,WallDensity	 !,LinearThick,ThickExtSh
	real(8) maxDiam
	real(8) thickness_slope,thickness_intercept,minThick,F_LENGTHMM
	character(2) profileType
	
	real(8) rtoll !tollerance in iterative procedure
	real(8) apot_p,apot_h,volume
	real(8), external :: polyProfile  !,polyInvert 
	real(8) shell_H1_start,shell_H2_start, r_ip1, r_ip2, rip2, rip1
	real(8) obstr_ip, obstr_edge, h1_start, h2_start !, rip2, rip1
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
	!namelist /CreateDiams/ 2012/03/23 merged to geosettings. 
	!-----
	!il programma calcola solo l'area geometrica (collecting area)
	!-----

	character(*),parameter :: nameSettingsFile='geoSettings.txt'
    
	namelist /geoSettings/  &
		nShells,maxDiam,focal_length,angular_shell_separation_deg,&
		shell_length,thickness_slope,thickness_intercept,minThick,WallDensity,profileType, &
		shell_H1_start,shell_H2_start
	!t=m*D+q	 

!--------------------------------------------2----------------------------------------
	print*,"Starting fortran program creaGeo on time:"
	CALL DATE_AND_TIME(DATE, TIME, ZONE)
	print*,date,time,zone 
	!default values
	!linearThick=2																																											 ! 
	maxDiam=0	!mancata lettura
	io=0	  
	temp=0
	rmed1=0
	rtoll=1d-10
	diamDecresce=.true.	
	minThick=0
	nShells=0
	
	!The following group will be not compiled by MS compiler.
	!MS$IF .false.
	nargs = command_argument_count()
	if (nargs /=0) call get_command_argument(1, diamFile) 	
	!!nargs = iargc()
	!!if (len(trim(diamFile))/=0) nargs=1
	!!call get_command_argument (1, diamFile)
	!!call getarg (1, diamFile)
	!MS$ENDIF

	open (10,file=nameSettingsFile,status='old',iostat=io)
	if (io/=0) then 
		print *,"Error in opening file ",nameSettingsFile
		print*,"iostat= ",io
		stop
	endif
	read (10,geoSettings)
	!read(10,geoPars)
	close(10)

	!if (focal_length/=0) f_length=focal_length*100
	!print*,"F_length (cm) = ",F_length	!HHH
	!HHH fine
	shell_length1=shell_length
	shell_length2=shell_length
	if (focal_length/=0) f_length=focal_length*100
	print*,"F_length (cm) = ",F_length	!HHH
	print*,"f_height1 (mm) = ",shell_length1	!HHH
	print*,"f_height2 (mm) = ",shell_length2	!HHH
	F_LENGTHmm=F_LENGTH*10

	!fissa nw (numero di shell)	da passare a geoinp
	if  (nShells/=0)	then
		if (nShells > nMaxSh) then
			print*,"The maximum number of shells is set to ",nmaxsh,","
			print*,"modify the program."  
			stop
		endif
		nw=nShells
		print*,"Shell number:",	nw
	else
		print*,"ERROR: namelist geopars not found in ",nameSettingsFile
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
	else if(maxDiam/=0)then
			print*,"Maximum diameter as read from parameters file: ", maxDiam
			rstart=maxDiam/2
	else
		print*,"Maximum diameter not found in namelist CreateDiams in ", nameSettingsFile
		stop
	end if
	print*,'Max aperture radius[mm]= ',rstart			

	tpeso=0. !inizializza peso
	tarea=0. !inizializza area

	r=rstart
!	shell_length1=H1_slope*r+H1_intercept
!	shell_length2=H2_slope*r+H2_intercept
!	H1_start=H1_start_slope*r+H1_start_intercept
!	H2_start=H2_start_slope*r+H2_start_intercept
!	!print*,"f_height1 (mm) = ",shell_length1	!HHH
	alkw=shell_length1*tan(angular_shell_separation_deg*PI/180.) !distanza tra uno specchio e l'altro per
		! evitare obreggiamenti fuori asse


	do i=1,nw+1 !ciclo do per costruzione degli specchi	e calcolo del peso
	!gliene faccio calcolare uno in piu' per usarlo per blocco centrale
	!*******************************************************		
		
		if (nargs==0 .or. i==nw+1 ) then
			!Set spes.
			!il raggio usato per il calcolo degli spessori e' quello esterno della
			!pupilla d'uscita, l'ultima shell fittizia ha spessore 0.
			if (i== nw+1) then
				spes=0
			else
				spes = (thickness_intercept+thickness_slope*2*r)/(1-thickness_slope*2)		 !il denominatore si rende necessario
				spes= max(spes,minThick)
			end if											!in quanto r e' diametro esterno (comprende spessore)	
			
			!calculate the sequence of rmed1.
			r=r-spes !radius at the mirror surface	

			!recursive procedure for r0 calculation from r at the entrance pupil		
			call polyCoeff(profileType,r,F_lengthmm,a1(i),a2(i),a3(i), &
						a4(i),a5(i),b1(i),b2(i),b3(i),b4(i),b5(i))
			r0=polyInvert(-shell_length1,r,a1(i),a2(i))
		
			do
				call polyCoeff(profileType,r0,F_lengthmm,a1(i),a2(i),a3(i), &
						a4(i),a5(i),b1(i),b2(i),b3(i),b4(i),b5(i)) 
				r1=polyInvert(-shell_length1,r,a1(i),a2(i))
				
				if(dabs(r1-r0)<=rtoll) exit			
				r0=r1
			end do
			rmed1(i)=r0
		end if
		r=rmed1(i) - alkw !external radius at exit pupil, used in next cycle
		
		theta(i)=datan2(rmed1(i),F_LENGTHmm)/4 !refinement. Assert theta(i)=th (se createFileShells=1)
		call polyCoeff(profileType,rmed1(i),F_lengthmm,a1(i),a2(i),a3(i), &
						a4(i),a5(i),b1(i),b2(i),b3(i),b4(i),b5(i))

!		shell_length1=H1_slope*r+H1_intercept
!		shell_length2=H2_slope*r+H2_intercept
!		H1_start=H1_start_slope*r+H1_start_intercept
!		H2_start=H2_start_slope*r+H2_start_intercept

		!calculation of the accurate value (2nd order) of radii 
		!at entrance and exit pupils corresponding to rmed
		!assert non dovrebbe cambiare se la procedura iterativa e' accurata.	
		!servono se calcola geometria a partire da diametri medi (createFileShells=2)
		!l'asse z e' orientato verso il piano focale, H e' negativo sulla parabola.
		!non e' necessario cambiare il segno di a2 e b2 perche' al quadrato
		rmax1(i)=polyProfile(-shell_length1,rmed1(i),a1(i),a2(i),a3(i),a4(i),a5(i))		
		rmin1(i)=polyProfile(shell_length2,rmed1(i),b1(i),b2(i),b3(i),b4(i),b5(i))	
		!radii at intermediate pupils	
		rip1=polyProfile(h1_start,rmed1(i),a1(i),a2(i),a3(i),a4(i),a5(i))
		rip2=polyProfile(-h2_start,rmed1(i),a1(i),a2(i),a3(i),a4(i),a5(i))		
		!obstr ip and obstr_edge are the extremes reflection heights on surface 1
		obstr_ip=max(h1_start,h2_start*(dtan(theta(i))-dtan(2*theta(i)))/(tan(2*theta(i))-dtan(3*theta(i))))
		obstr_edge=min(shell_length1,shell_length2*(dtan(theta(i))-dtan(2*theta(i)))/(tan(2*theta(i))-dtan(3*theta(i))))
		!volume=spes*pi*((rmax1(i)+rmed1(i))*apot_p+(rmed1(i)+rmin1(i))*apot_h)
		!r_ip1=rmed1(i)+h1_start/dtan(theta(i))
		!r_ip2=rmed1(i)-h2_start/dtan(3*theta(i))   

		!la formula calcola lo spessore in base al diametro alla pupilla di ingresso
		!il denominatore (1-2m) serve perche' si calcola in base al raggio esterno
		!mentre la linearita' e' con il raggio interno.	Anche questo non dovrebbe cambiare
		!if (CreateFileShells==2) thick(i)=(thickness_intercept+thickness_slope*2*rmax1(i))/(1-thickness_slope*2)
		spes=(thickness_intercept+thickness_slope*2*rmax1(i))/(1-thickness_slope*2)
		spes=max(spes,minThick)
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
			!peso= pi*((rmax1(i)+rmed1(i))/dcos(theta(i))*shell_length1+&
			!			(rmed1(i)+rmin1(i))/dcos(3*theta(i))*shell_length2)&
			!			*thick(i)*WallDensity/1000000.
			apot_p=shell_length1/dcos(theta(i))
			apot_h=shell_length2/dcos(3*theta(i))
			!volume=spes*pi*((+(rmed1(i)+obstr_ip/dtan(theta(i))))*apot_p+(rmed1(i)+rmin1(i))*apot_h)
			volume=spes*pi*((rmed1(i)+rmax1(i))*apot_p+(rmed1(i)+rmin1(i))*apot_h)
			peso=volume*WallDensity/1000000.			
			weight(i)=peso
			!area=pi*(rmax1(i)**2-rmed1(i)**2)/100 ! collecting area
			!area=pi*((rmax1(i)-obstr_edge/dtan(theta(i)))**2-(rmed1(i)+obstr_ip/dtan(theta(i)))**2)/100 ! collecting area
			!2019/01/17 it seems that previous version was calculating only double reflection area(maybe?)
			!do projected area
			area=pi*(rmax1(i)**2-max(rmed1(i),rmax1(i+1)+spes)**2)/100 ! collecting area, approx in spes, should be next shell
		
			!print *,obstr_ip,obstr_edge
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
	write (112,'(a120)') "Nshell   Dmax(mm)   Dmid(mm)"// &
		"   Dmin(mm)   thickness(mm)   Angle(rad)   Area(cm^2)"// &
		"  Mass(kg)	Shell_length1	Shell_length2"
	11 format  (i4,9f15.6)
	do k=1,nw
		write (112,11) k, rmax1(k)*2,rmed1(k)*2,rmin1(k)*2,thick(k),theta(k), acoll(k),weight(k),shell_length1,shell_length2
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

	!set a and b cofficients according to the profile type.
	!theta is positive.
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

!if delta is not zero (double cone), take positive solution,
! the other will be negative. Indeed Rz(R0)=R0**2+a1*R0*z+a2*z**2
! is parabola with vertex down with negative x -a1*z/2,
! then for each Pz there are two possible solutions R0, > and < of 0.
	 
end function

function polyProfile(z,r0,c1,c2,c3,c4,c5)	
implicit none
real(8), intent(in):: z,r0,c1,c2,c3,c4,c5
real(8) polyProfile,t 
!return the radius for a given z, from the five coefficients of the polynomial profile
t=z/r0
polyProfile=r0*sqrt(1+c1*t+c2*t**2+c3*t**3+c4*t**4+c5*t**5)
	 
end function

function linearize(X,y0,yN,m,q,minvalue,maxvalue)
!Return an array of values obeying to a linear relation. The relation can
!	be set by a couple of parameters among y0,yN,m,q.

!X		X is a list of values for the x in the
!x0		first value
!xN		last value
!m		linear coefficient
!q		linear coefficient
!min	if provided all values lower than this are set to min
!max	if provided all values higher than this are set to max
implicit none
real(8) X(:)  !assumed-shape array
real(8) linearize(size(X)),ymin,ymax
real(8),optional::y0,yN,m,q,minvalue,maxvalue
integer argcount,npoints,i

!check number of arguments
argcount=0
if (present(y0)) argcount=argcount+1
if (present(yN)) argcount=argcount+1
if (present(m)) argcount=argcount+1
if (present(q)) argcount=argcount+1
if (argcount /= 2) then
	print*, "Wrong number of argument provided to function LINEARIZE (nargs=",argcount,","
	print*, "must be two!)"
	print*,"y0=",y0
	print*,"yN=",yN
	print*,"m=",m
	print*,"q=",q
	stop
endif 

npoints=size(x)
!calculate m and q (if m and q are provided proceeds).
if (present(y0)) then
	if (present(yn)) then
		m=(y0-yn)/(X(npoints)-X(1))
	else if (present(q)) then
		m=(y0-q)/X(1)
	else 
		if (.not. present(m)) then
			print*,"Y0 only argument: this should never happen, check LINEARIZE ROUTINE!!"
			stop
		endif
	endif
	q=y0-x(1)*m
else if (present(yN)) then
	if (present(q)) then
		m=(yN-q)/X(Npoints)
	else if (present(m)) then
		q=yN-x(npoints)*m
	else	
		print*,"YN only argument: this should never happen, check LINEARIZE ROUTINE!!"
	endif
endif

linearize=X*m+q
if (.not. present(minvalue)) then
	ymin=minval(linearize)
else 
	ymin=minvalue
endif
if (.not. present(maxvalue)) then
	ymax=maxval(linearize)
else
	ymax=maxvalue
endif
if (ymin > ymax) then
	print*, "ERROR: MINVALUE>MAXVALUE in function LINEARIZE"
	stop
endif 

do i=1,npoints
	if (linearize(i) >	ymax) linearize(i)=ymax
	if (linearize(i) <	ymin) linearize(i)=ymin
end do

end function



