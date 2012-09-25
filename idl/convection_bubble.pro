;----------------------------------------------------------------------------->

function get_path,solmodel=solmodel

if keyword_set(solmodel) then retval='~/science/procedures/gen/modelling/bubble/solar_model/'




return,retval

end

;----------------------------------------------------------------------------->

function get_param, hphot=hphot, boxsz=boxsz, dxbin=dxbin, dybin=dybin, bubrxy=bubrxy, $
	bubmass=bubmass, radsun=radsun, masssun=masssun, hinit=hinit, tempinit=tempinit, $
	gravconst=gravconst, boltz=boltz, mu=mu, $
	nstep=nstep, gama=gama, dtbin=dtbin, hmin=hmin

;Resolution of box
if keyword_set(dxbin) then retval=1d8 ;(1Mm/px) ;bin size in cm/px
if keyword_set(dybin) then retval=1d8 ;(1Mm/px) ;bin size in cm/px

;Size of box
if keyword_set(boxsz) then retval=[100.,300.] ;in [px] = (cm/DXBIN)

;Min height of box
if keyword_set(hmin) then retval=500. ;in [px] = (cm/DXBIN)

if keyword_set(nstep) then retval=100. ;number of time steps
if keyword_set(dtbin) then retval=3600. ;in seconds

if keyword_set(hphot) then retval=800. ;in [px] = (cm/DXBIN)

if keyword_set(radsun) then retval=695.5*1d8 ;in [cm]
if keyword_set(masssun) then retval=1.9889*1d33 ;in [g]

if keyword_set(gravconst) then retval=6.674*1d-8 ;in [cm^3 / g s^2]
if keyword_set(boltz) then retval=1.38d-16 ;in [erg/K]
if keyword_set(mu) then retval=1.66*1d-24 ;in [g/particle] <- assume mass of hydrogen
if keyword_set(gama) then retval=5./3. ;to do with the degrees of freedom for gas particles

if keyword_set(bubrxy) then retval=[20.,20.] ;[x,y] radii - elipse parameters for initial bubble in [cm*DXBIN]

if keyword_set(bubmass) then retval=60d6 ;grams (estimated from data)

if keyword_set(hinit) then retval=200. ;in [px] = (cm/DXBIN)

if keyword_set(tempinit) then retval=20d3 ;K 


return,retval

end

;----------------------------------------------------------------------------->

pro get_profiles, radsun, tempsun, rhosun, presssun, nocal=nocal

path=get_path(/solmodel)

readcol,path+'bp2004stdmodel.txt',dm,dr,dt,drho,dp,skip=23,format='F,F,F,F,F,F,F,F,F,F,F,F',delim=' '
if keyword_set(nocal) then begin
	radsun=dr
	tempsun=dt
	rhosun=drho
	presssun=dp

	return
endif

rad=dr*get_param(/radsun) ;in Mm

radsun=(findgen((get_param(/boxsz))[1])+get_param(/hmin))*get_param(/dybin) ;in cm
radsun=radsun/1d8  ;in Mm

;Fit end of data to edge of Sun
plot,dr,dt,/ylog,chars=2
print,radsun[0]/get_param(/radsun)
print,get_param(/radsun)







;estimated from data
tempsun=-6359.71*radsun+950000. ;in K
tempsun[where(radsun gt 695.)]=0.

rho=drho; [g/cm^3]   /1d6*(1d8)^2 ; [Mg/Mm^2] = ( g / cm^2 ) * ( 1E8 cm/Mm )^2 / ( 1E6 g/Mg)

rhosun=-0.00025*radsun+0.18 ;in g/cm^2
rhosun[where(radsun gt 695.)]=0. ;in g/cm^2

press=dp

;tempsun=interpolate(rad,)

;density profile

;temperature profile

;pressure profile (derived from other two)



stop



end

;----------------------------------------------------------------------------->

function get_bubble

;get bubble parameters
bubrxy=get_param(/bubrxy)

;generate elipse for bubble
;ellipse, x, [x0,y0,a,b], y, dyda



return,1

end

;----------------------------------------------------------------------------->

pro box_init



if not keyword_set(hphot) then hphot=get_param(/hphot)

if not keyword_set(boxsz) then sz=get_param(/boxsz)




;density, pressure, height, temperature

;DENSITY


;fltarr(sz)


end

;----------------------------------------------------------------------------->

;calculate the acceleration of gravity for a given solar radius
;subtract the mass above H from the mass of the Sun.
;Calculate mass above assuming a density profile
;H - height from center of sun in cm
;RADSUN - array of radii from solar center in cm
;RHOSUN - matching array of densities from solar center in g/cm^3
;IRREGULAR - set if the grid spacing of the radial elements is irregular.

function calc_gravity, h, radsun, rhosun, irregular=irregular

wh=(where(abs(radsun-h) eq min(abs(radsun-h))))[0]
if radsun[wh] gt h then hmass=radsun[wh:*] $
	else hmass=radsun[wh+1:*]

msun=get_param(/masssun)
grav=get_param(/gravconst)

nh=n_elements(hmass)
if keyword_set(irregular) then begin
	mtot=4./3.*!pi*(hmass[1]^3.-hmass[0]^3.)*rhosun[0] ;since i starts at 1, estimate initial shell
	for i=1,nh-1 do begin
		mshell=4./3.*!pi*(hmass[i]^3.-hmass[i-1]^3.)*rhosun[i-1]
		mtot=mtot+mshell
	endfor
endif else begin
	dh=hmass[1]-hmass[0]
	mtot=total(4./3.*!pi*(hmass^3.-(hmass-dh)^3.)*rhosun) ;shells calculated using h-dh instead of h+dh
endelse

mbelow=msun-mtot ;mass below H contributing to g(h)

gh=grav*mbelow/(h^2.) ;acceleration of gravity at specified radius from sun center.

return, gh
end

;----------------------------------------------------------------------------->


pro convection_iterate

nstep=get_param(/nstep)
dtbin=get_param(/dtbin)

bubmass=get_param(/bubmass) ;in g
cgas=get_param(/boltz)/get_param(/mu) ;k/mu
gama=get_param(/gama)
radinit=(get_param(/bubrxy))[0]*get_param(/dxbin) ;in cm
hinit=get_param(/hinit)*get_param(/dybin) ;in cm
tempinit=get_param(/tempinit) ;in K
densinit=bubmass*(4./3.*!pi*(radinit)^3.) ;in g/cm^3
pressinit=densinit*cgas*tempinit

get_profiles, radsun, tempsun, rhosun,presssun0 ;in cgs
presssun=rhosun*cgas*tempsun

;TEMP!!!! compare gas law to modelled pressure.
plot,rhosun,presssun
oplot,rhosun,presssun0,color=!red
stop

tarr=-dtbin ;initialise so that at first iteration, t=0
accarr=0.
velarr=0.
harr=hinit
temparr=tempinit
densarr=densinit
pressarr=pressinit
radarr=radinit

for i=1,nstep-1 do begin
	
	thist=tarr[i-1]+dtbin
	
	lasthe=harr[i-1]
	lasttemp=temparr[i-1]
	lastdense=densarr[i-1]
	lastpress=lastdense*cgas*lasttemp	
	
	lastind=where(abs(radsun-lasthe) eq min(abs(radsun-lasthe)))
	lastdensee=rhosun[lastind]
	lasttempe=tempsun[lastind]
	
	
stop	
	;calculate the acceleration of gravity
	thisgh=calc_gravity(lasthe, radsun, rhosun)
	
	
	thisacc=(1./bubmass-1./(lastdensee*4./3.*!pi*(radarr[i-1])^3.))*(1./(lastdensee*4./3.*!pi*(radarr[i-1])^3.)-1./bubmass)
	
	thisvel=velarr[i-1]+thisacc*dtbin
	
	deltah=thisvel*dtbin
	thishe=lasthe+deltah
	
	thisind=where(abs(radsun-thishe) eq min(abs(radsun-thishe)))
	thisdensee=rhosun[thisind]
	thistempe=tempsun[thisind]
	thispress=thisdensee*cgas*thistempe
	
	thisrad=(bubmass*lasttemp*cgas/(4./3.*!pi*lastpress))^(1./3.)*(lastpress/thispress)^((1.-gama)/2./gama)
	
	thistemp=lasttemp*(lastpress/thispress)^((1.-gama)/gama)
	
	tarr=[tarr,thist]
	accarr=[accarr,thisacc]
	velarr=[velarr,thisvel]
	harr=[harr,thishe]
	temparr=[temparr,thistemp]
	densarr=[densarr,4./3.*!pi*thisrad^3.]
	pressarr=[pressarr,thispress]
	radarr=[radarr,thisrad]

;TEMP!!!!!!!	
	print,thisacc,thisvel,deltah*get_param(/dybin)
	wait,.2
endfor






stop
end

;----------------------------------------------------------------------------->

pro convection_bubble

;initialise parameters



;initialise bubble


;initialise box




;iterate model
convection_iterate






end

;----------------------------------------------------------------------------->