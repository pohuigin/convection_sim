""" 
File:	convect.py
Author: Sean Blake
Date: 	December 2012 - January 2013

This code uses data from the 'Model S' solar data set. This can be found at: http://users-phys.au.dk/jcd/solar_models/

The code can be (roughly) separated into 3 parts. 

1) Data extraction + definition of the functions used.
-lines 18-340

2) The 'For' loop used for the motion of the cell.
-lines 340-462

3) The various plots used for my project. These plots are all commented out. In order to use a particular plot, comment that plot back into the code, and run again.
-lines 462-930
"""

import pdb
import math
import scipy
from scipy import interpolate
import matplotlib.font_manager
from time import *

t0=clock()
rc('text', usetex=True)												#Allows 'tex' in my graphs
rc('font',**{'family':'times','sans-serif':['times']})				#Uses times font in my plots.

############################################## Data Extraction ##############################################
a1, a2, a3, a4, a5 = np.loadtxt('bleh.txt', usecols = (0, 1, 2, 3, 4), unpack=True, skiprows = 0)		#extracting data from main GONG file.

adgrad = a1[2::5]
opacity = a3[1::5]
opacity = opacity[::-1]
opacity= list(opacity)
del opacity[-1]
opacity = array(opacity)

t1 = a3[0::5] 						#temperature
t1 = t1[::-1]
t1 = t1[1::]

r1 = a1[0::5]						#radius
r1 = r1[::-1]
r1 = ((r1/100)/6.96342E8)
r1 = r1[1::]

m1 = a2[0::5]						#mass
m1 = m1[::-1]
m1 = exp(m1)
m1 = m1
m1 = m1[1::]

rho1 = a5[0::5]						#density
rho1 = rho1[::-1]
rho1 = rho1 * 1E3
rho1 = rho1[1::]

P1 = a4[0::5]						#Pressure
P1 = P1[::-1]
P1 = P1*(0.1)
P1 = P1[1::]

Cp = a3[3::5]						#Specific heat capacity at constant pressure
Cp = Cp[::-1]
Cp = Cp[1::]


gamma2 = np.loadtxt('zgamma.txt', unpack=True, skiprows=0)
r2 = np.loadtxt('zradius.txt', unpack=True, skiprows = 0)
t2 = np.loadtxt('ztemp.txt', unpack=True, skiprows = 0)				#these are from the secondary data file.
rho2 = np.loadtxt('zrho.txt', unpack=True, skiprows = 0)
P2 = np.loadtxt('zpressure.txt', unpack=True, skiprows = 0)
m2 = np.loadtxt('zmass.txt', unpack=True, skiprows = 0)

gamm=[]										#allows gamma (from secondary set of data) to be used with first set of data
for i in arange(0, len(r1), 1):
	y = interp(r1[i], r2, gamma2)
	gamm.append(y)
gamma1 = array(gamm)

P2 = P2 * 0.1             #Converting Dynes to Pa
rho2 = rho2 * 1E3         #Converting g/cm^3 to kg/m^3

############################################## Constants ##############################################

radiussol = 6.96342E8	#Radius of Sun
kxc = 1.3806488E-23     #Boltzmanns constant
uu = 1.660538E-27       #Atomic mass constant
mu = 0.6                #average mass of sun
G = 6.67428E-11         #Gravitational Constant
gasconstant = 8.3144621	#Gas Constant
sbconst = 5.670400E-8	#Stefann-Boltzmann Constant
c = 3E8					# Speed of light
solarmass = 1.9891E30	#solar mass in kg
xx1 = r1 * radiussol	#makes array for radius of sun in metres as opposed to fraction (as with r1)
xx2 = r2 * radiussol
"""----------------------------------------------------------------------------------------------------------------------"""

############################################## Adjustments of constants/user inputs ##############################################
#These are used if you want to particularly specify a cell's conditions.

m_cell = float(raw_input("Mass (Mkg): "))            	#Enter mass of cell in mega-kilograms
r_initial1 = float(raw_input("Radial Fraction: "))    	#Radial Fraction of cell
t_initial = float(raw_input("Temperature (MK): "))   	#Temperature of cell

#m_cell = 1000
#t_initial = 2.5						
#r_initial1 = 0.9
#r_initial2 = 0.8

m_cell = m_cell * 1E6									#puts mass of cell into mega-kgs.
t_initial = t_initial *1E6								#temperature of cell into mega-kelvin
n = (m_cell * 1000)/2.0158								#number of moles in cell
R = 8.3144621											#Gas constant
rx1 = r_initial1 * 6.96342E8							#initial radius of cell in metres
"""----------------------------------------------------------------------------------------------------------------------"""

############################################## 	Various Functions  	##############################################

def sunvol(xx1):										#calculates radius of sun at a particular point 'xx1'.
	return (4./3)*pi*(xx1)**3
 
def idealtemp(P1):										#calculates the 'ideal' temperature of the Sun, given pressure and density by Model S.
	return ((P1 * mu * uu) / (rho1 * kxc))				#This is a measure of how much the Model S Sun acts like an ideal gas.
	
def idealdens(P1, t1):									#calculates the 'ideal' density of the Sun, given pressure and temperature by Model S.
	return ((P1 * mu * uu) / (t1 * kxc))				#This is also a measure of how much the Model S Sun acts like an ideal gas.

def dens(P, temp):                                 		#Calculating density of cell for pressures 'P' and cell temperatures via ideal gas law.
    return ((P * mu * uu) / (temp * kxc))        

def volbub(m_cell, P, temp):                        	#Calculating volume of cell for all using the density calculated above.
	aaa = m_cell
	bbb = (P * mu * uu)
	ccc = (temp * kxc)
	ddd = bbb/ccc
	eee = aaa/ddd
	return eee

def radiuscell(m_cell, P, temp):         				#Radius of the cell calculated from volume above. Assumes spherical cell.
	aaa = 3./(4. * pi)
	bbb = aaa * volbub(m_cell, P, temp)
	ccc = bbb ** (1./3.)
	return ccc

def areacell(m_cell, P, temp):							#Surface area of Cell. Calculated from cell radius above.
	return (4 * pi * radiuscell(m_cell, P, temp)**2)

def g(m, z):												#Function for gravity at all positions within the Sun.
    return (G *m * 1.9891E30/ (z)**2)

def buoy(m, z, m_cell, P, t_initial, rho):					#Calculating buoyancy of cell for all r. Function is broken down line by line as:
	aaa = (G *m * 1.9891E30/ (z)**2)						#aaa = Gravity (acceleration) at position
	bbb = (m_cell / ((P * mu * uu) / (t_initial * kxc)))	#bbb = volume at position
	ccc = rho												#ccc = ambient density at position
	ddd = aaa * bbb * ccc 									#combining above for buoyancy force         
	eee = (ddd)/m_cell										#eee = buoyancy acceleration (divide by mass)
	return eee    

def netacceleration(m, z, m_cell, P, t_initial, rho): 		#Buoyancy acceleration minus gravity. Broken down as follows:
	aaa = (G *m * 1.9891E30/ (z)**2)						#aaa = Gravity (acceleration) at position
	bbb = (m_cell / ((P * mu * uu) / (t_initial * kxc)))	#bbb = volume at position
	ccc = rho												#ccc = ambient density at position
	ddd = aaa * bbb * ccc									#combining above for buoyancy force
	eee = (ddd)/m_cell										#eee = buoyancy acceleration
	fff = eee - aaa											#fff = buoyancy acceleration -gravity acceleration = net acceleration
	return fff

def pheight(P, rho, m, x):									#Pressure Scale Height calculated for each point.
	return (P/(rho * g(m, x)))

############################################## Calculating Differences Between top/bottom of Cell ##############################################
""" This was used for plot number (). It works like so:

The cell radius is found. This is added and taken away from the cells position (at cell's centre) to get position at top and bottom of cell.
The conditions at these two points (temperature, pressure, density) are found. The difference between these conditions are found, then divided by the value of the condition at the centre of the cell.
This gives a rough way of measuring how well the spherical cell model works as the cell rises."""

def rminus(z, m_cell, P, t_initial):       					#position of cell minus radius, or position of bottom of cell
	return (z - (((3 * (m_cell / ((P * mu * uu) / (t_initial * kxc))))/4*pi)**(1./3.)))

def rplus(z, m_cell, P, t_initial):       					#position of cell plus radius, or position of top of cell
	return (z + (((3 * (m_cell / ((P * mu * uu) / (t_initial * kxc))))/4*pi)**(1./3.)))

def Pdifference(xx1, m_cell, P1, t_initial):      			#difference in pressure between these two points
	plus   = rplus(xx1, m_cell, P1, t_initial)				
	minus  = rminus(xx1, m_cell, P1, t_initial)				
	Pplus  = interp(plus, xx1, P1)							#pressure at top of cell
	Pminus = interp(minus, xx1, P1)							#pressure at bottom of cell
	return ((Pminus-Pplus)/P1)								#comparison to centre of cell

def rhodifference(xx1, m_cell, P1, t_initial):    			#This does for density what Pdifference(...) does for pressure.
	plus   = rplus(xx1, m_cell, P1, t_initial)
	minus  = rminus(xx1, m_cell, P1, t_initial)
	rhoplus  = interp(plus, xx1, rho1)
	rhominus = interp(minus, xx1, rho1)
	return ((rhominus-rhoplus)/rho1)

def tdifference(xx1, m_cell, P1, t_initial):      			#This does for temperature what Pdifference(...) does for pressure.
	plus   = rplus(xx1, m_cell, P1, t_initial)
	minus  = rminus(xx1, m_cell, P1, t_initial)
	tplus  = interp(plus, xx1, t1)
	tminus = interp(minus, xx1, t1)
	z = ((tminus-tplus)/t1)
	return z

############################################## Adiabatic Temperature Gradient  ##############################################
""" There were a number of different ways used to calculate the adiabatic temperature gradient of a cell. They were numbered 1-4, and found using the following functions."""

def adbtempgrad1(gamma1, t1, P1):						#first adiabatic method		
	a = (gamma1-1)/(gamma1)
	b = (t1/P1)
	c = -((g(m2, xx2)) * rho2)
	return a * b * c

def adbtempgrad2(m2, xx2, rho2, t2, P2):				#second adiabatic method
	a = -((g(m2, xx2)) * rho2)
	b = (t2/P2)
	c = Cp
	return a * b * c	

def adbtempgrad3(m1, xx1, Cp):							#third adiabatic method
	a = -g(m1, xx1)										
	b = Cp+1E-25										#The '1E-25' number was added, to prevent the code from thinking it was dividing by 0.
	return a/b

#The 4th adiabatic method required d(rho)/dR. This was found using a for loop with the geometric equation for slope.
slop=[]													#An empty dummy list 'slop' was made, which could be filled later.
for i in arange(0, len(r1), 1):							#for every point in the radius array r1,
	zzz = r1

	if i == 1:											#if it was the first point, the slope was set to 0.
		slo = 0
	if i == len(r1)-1:									#For the last point, the last point and the 3rd last point were used.
		x1 = (zzz[i-2] * radiussol)
		x2 = (zzz[i] * radiussol)
		y1 = interp(x1, xx1, rho1)
		y2 = interp(x2, xx1, rho1)
		slo = (y2-y1)/(x2-x1)
	else:												#For all other points, the points to the left and right of it were used (named x1, x2)
		x1 = (zzz[i-1] * radiussol)						
		x2 = (zzz[i+1] * radiussol)
		y1 = interp(x1, xx1, rho1)						#The corresponding density values were found (y1 and y2 respectively)
		y2 = interp(x2, xx1, rho1)
		slo = (y2-y1)/(x2-x1)							#These were plugged into the geometric slope equation 
	slop.append(slo)									#This value was put into the dummy list created earlier
rhoslope = array(slop)									#List was renamed and made an array so it could be manipulated.

# The following For loop finds the pressure gradient (dP/dR), much like the loop above.
#This isn't actually used in the code, however, and may best be commented out to speed it up.
slop=[]													
for i in arange(0, len(r1), 1):
	zzz = r1

	if i == 1:
		slo = 0
	if i == len(r1)-1:
		x1 = (zzz[i-2] * radiussol)
		x2 = (zzz[i] * radiussol)
		y1 = interp(x1, xx1, P1)
		y2 = interp(x2, xx1, P1)
		slo = (y2-y1)/(x2-x1)
	else:
		x1 = (zzz[i-1] * radiussol)
		x2 = (zzz[i+1] * radiussol)
		y1 = interp(x1, xx1, P1)
		y2 = interp(x2, xx1, P1)
		slo = (y2-y1)/(x2-x1)
	slop.append(slo)
Pslope = array(slop)

def adbtempgrad4(gamma1, t1, P1, rho1):				#4th adiabatic method, using the density gradient 'rhoslope'.
	a = (gamma1-1)
	b = (t1)/(P1)
	c = (P1)/(rho1)
	d = rhoslope
	return a * b * c * d
############################################## Actual Temperature Gradient  ##############################################

#This method calculates the actual temperature gradient.
slop=[]												#Dummy list, ready to be populated.
for i in arange(0, len(r1), 1):						#for every point in the radius array r1,
	zzz = r1

	if i == 1:										#if it was the first point, the slope was set to 0.
		slo = 0
	if i == len(r1)-1:								#For the last point, the last point and the 3rd last point were used.
		x1 = (zzz[i-2] * radiussol)
		x2 = (zzz[i] * radiussol)
		y1 = interp(x1, xx1, t1)
		y2 = interp(x2, xx1, t1)
		slo = (y2-y1)/(x2-x1)
	else:											#For all other points, the points to the left and right of it were used (named x1, x2)
		x1 = (zzz[i-1] * radiussol)				
		x2 = (zzz[i+1] * radiussol)
		y1 = interp(x1, xx1, t1)					##The corresponding temperature values were found (y1 and y2 respectively)
		y2 = interp(x2, xx1, t1)
		slo = (y2-y1)/(x2-x1)
	slop.append(slo)								#This value was put into the dummy list created earlier
tslope = array(slop)								#List was renamed and made an array so it could be manipulated.

############################################## Adiabatic Temperature Change  ##############################################
#Calculates the temperature of a cell or cells as they rise through the Sun.

rx1 = r_initial1 * radiussol				#Position of cell's in metres.
tempfina, distanc=[], []					# Two dummy lists, to be populated later.
asdf = interp(r_initial1, r1, t1)			#'asdf' is the temperature of the surroundings at the cell's position.  
temperatures=(2*asdf, 1.25*asdf, 1.1*asdf, 1.01*asdf)	#This gives a list of cell temperatures- 2, 1.25, 1.1 and 1.01 times the ambient temp.

for z in temperatures:						#For every temperature listed above:
	tempfina=[]	
	for i in (xx1):							#For every position along the solar radius:
		if i< (r_initial1*radiussol):		#If the position is less than the initial cell radius, the temperature is set negligibly small.
			temp2 = 0.01
		if i > (r_initial1*radiussol):		#If the position is more than the initial cell radius,
			temp1 = z						#The conditions at these points are found...
			PP1 = interp(rx1, xx1, P1)
			gammma1 = interp(rx1, xx1, gamma1)
			gammma2 = interp(i, xx1, gamma1)
			PP2 = interp(i, xx1, P1)

			aaa = (PP1)**(1-gammma1)		#And plugged into an equation for ideal gas temperature.
			bbb = (temp1)**(gammma1)
			ccc = (PP2)**(1-gammma2)
			ddd = aaa * bbb
			eee = ddd / ccc
			temp2 = (eee)**(1/gammma2)		#This is returned here.

		tempfina.append(temp2)				#And plugged into the list 'tempfina'

	if z == 2*asdf:							#These next 4 'ifs' put the temperatures into lists numbered 1-4.
		tempfinal1=array(tempfina)			#This is so the temperature loss of cells for different initial temperatures can be compared.
	if z == 1.25*asdf:
		tempfinal2 = array(tempfina)
	if z == 1.1*asdf:
		tempfinal3 = array(tempfina)	
	if z == 1.01*asdf:
		tempfinal4 = array(tempfina)
		

"""----------------------------------------------------------------------------------------------------------------------"""
"""----------------------------------------------------------------------------------------------------------------------"""
"""----------------------------------------------------------------------------------------------------------------------"""
"""=================================================================================================="""

############################################## BIG BAD FOR/WHILE LOOPS ##############################################

"""=================================================================================================="""

time = 0.1											#time step chosen- 0.1 seconds.
u, d, s = 0, 0, 0									#initial speed 'u'=0, counter 'd'=0, and displacement 's'=0

timelist, distlist, speedlist, acclist = [0], [0], [0], [0] #This will create lists for time, cell displacement, cell speed and cell acceleration.

rxx1 = (r_initial1 + 0.01) * 6.96342E8				#initial radius of cell in metres. 0.01 was added, otherwise cell begins erroneously sinking.
ghjk = 0
kjhg = 0

jj = [0.9]											#Cell positions to be tested. Works nicely if you use an array such as arange(0.7, 0.95, 0.01)
kk = (2, 1.25, 1.1, 1.01)							#Temperature factors to be used. These are multiples of solar temp at initial cell radius.

average=[]											#The following set empty lists for average cell speeds, cell positions, temperatures, and 
posit =[]											#final speeds.
temper=[]
finalspee =[]
for j in jj:										#For each cell position
	for kj in kk:									#For each temperature factor in 'kk'

		lop = interp(j, r1, t1)						#This finds the solar temp at initial cell radius.
		k = kj * lop								#Multiplies temperature factor by solar temperature.
		kjhg = kjhg + 1								#counter for telling how complete the loop is.

		tempdif = 0
		tempfina, distanc=[], []
		
		#The following finds the temperature loss (as above) for every cell used in the loop.
		rx1 = j * radiussol
		for i in (xx1):
			if i< (j*radiussol):
				temp2 = 0.01
			if i > (j*radiussol):
				temp1 = k
				PP1 = interp(rx1, xx1, P1)
				gammma1 = interp(rx1, xx1, gamma1)
				gammma2 = interp(i, xx1, gamma1)
				PP2 = interp(i, xx1, P1)

				aaa = (PP1)**(1-gammma1)
				bbb = (temp1)**(gammma1)
				ccc = (PP2)**(1-gammma2)
				ddd = aaa * bbb
				eee = ddd / ccc
				temp2 = (eee)**(1/gammma2)
			tempfina.append(temp2)
		tempfinal=array(tempfina)	
		#This temperature array 'tempfinal' is used in the array below.


		rxx1 = (j+ 0.001) * 6.96342E8				#rxx1 is the changing position of the cell used in the code.
		u, d, s = 0, 0, 0							#resets initial values for speed, counter and displacement.
		templist, kineticlist, timelist, distlist, speedlist, acclist = [0], [0], [0], [0], [0], [0]	#making the lists empty again.
		
		while (rxx1/6.96342E8)< 1: 					#while the cell is beneath the solar surface (where r/solarradius=1)
			#These find the the conditions at whatever position the cell is at
			mm1 = interp(rxx1, xx1, m1)				#mass of Sun
			PP1 = interp(rxx1, xx1, P1)				#pressure of Sun
			blaergh = interp(rxx1, xx1, tempfinal)	#Temperature of cell
			rhorho1 = interp(rxx1, xx1, rho1)		#density of Sun
	
			a = netacceleration(mm1, rxx1, m_cell, PP1, blaergh, rhorho1)	#This finds the net acceleration acting upon the cell.
			interp(rxx1, xx1, netacceleration(m1, xx1, m_cell, P1, blaergh, rho1))		

			if u<0:									#this condition stops the loop if the cell begins sinking (u= initial velocity)
				average.append(0)					#this adds a 0-value to the average cell velocity.
				posit.append(j)						#position
				temper.append(k)					#temperature
				finalspee.append(0)					#final cell speed
				break 

			v = a * time + u                        # First equation of motion, where v=final velocity
			s = u * time + 0.5 * a * time**2        # Second equation of motion, where s=displacement
			d = d + 1                               #counter
			u = v									#sets new initial velocity as the final velocity.
			v = 0                            		# sets final velocity after step as initial velocity
			rxx1 = rxx1 + s                   		# sums total distance travelled so far, and changes rxx1 accordingly
			zxc = (0.5)*(m_cell)*(u * u)		
			distlist.append(rxx1)					#populates the distance, speed, acceleration, and time lists.
			speedlist.append(u)
			acclist.append(a)
			timelist.append(d*time) 
			kineticlist.append(zxc)
			
		tyu = sum(speedlist)/(d/10.)				#finds average velocity for cell.

		if speedlist[-1] > 0:						#If the final cell velocity is NOT 0, the cell has reached the surface.
			posit.append(j)							#position, temperature, final speed and average speed are added to lists.
			temper.append(k)
			finalspee.append(speedlist[-1])	
			average.append(tyu)

		if kj==(2):
			kineticlist1, timelist1, acclist1, speedlist1, distlist1 = kineticlist, timelist, acclist, speedlist, distlist
			kineticarray1, distarray1, speedarray1, accarray1, timearray1 = array(kineticlist1), array(distlist1), array(speedlist1), array(acclist1), array(timelist1)
			templist, kineticlist, timelist, distlist, speedlist, acclist = [0], [0], [0], [0], [0], [0]
			
		if kj==(1.25):
			kineticlist2, timelist2, acclist2, speedlist2, distlist2 = kineticlist, timelist, acclist, speedlist, distlist
			kineticarray2, distarray2, speedarray2, accarray2, timearray2 = array(kineticlist2), array(distlist2), array(speedlist2), array(acclist2), array(timelist2)
			templist, kineticlist, timelist, distlist, speedlist, acclist = [0], [0], [0], [0], [0], [0]
			
		if kj==(1.1):
			kineticlist3, timelist3, acclist3, speedlist3, distlist3 = kineticlist, timelist, acclist, speedlist, distlist
			kineticarray3, distarray3, speedarray3, accarray3, timearray3 = array(kineticlist3), array(distlist3), array(speedlist3), array(acclist3), array(timelist3)

		if kj==(1.01):
			kineticlist4, timelist4, acclist4, speedlist4, distlist4 = kineticlist, timelist, acclist, speedlist, distlist
			kineticarray4, distarray4, speedarray4, accarray4, timearray4 = array(kineticlist4), array(distlist4), array(speedlist4), array(acclist4), array(timelist4)
	
print "Process Time: ", clock() -t0, "seconds"		#prints time taken for code.

"""=================================================================================================="""

############################################## PLOTS AND SHITE ##############################################

"""=================================================================================================="""
"""
#######################################################################################
#PLOT NUMBER 1
# TEMPERATURE PRESSURE DENSITY MASS OF SUN (Surroundings)

close()
fig = figure(2)

bleh=dict(wspace= 0.25, hspace=0.2)
subplots_adjust(**bleh)

ax1 = subplot(221)
#setp( ax1.get_xticklabels(), visible=False)

title('(A)', fontsize=15)
plot(r1, t1, 'r', lw = 3.5)
xlabel(r'$r/R_\odot$', fontsize=17)
xlim([0, 1])
ylabel(r'Temperature ($K$)', fontsize = 15)
axvline(0.713, color='black', lw=2, linestyle='steps--')
axvline(0.2, color='black', lw=2, linestyle='steps--')
grid(True)

plt.show()

#####2


ax2 = subplot(222)
#setp( ax2.get_xticklabels(), visible=False)

title('(B)', fontsize=15)
plot(r1, P1, 'b', lw=3.5)
xlabel(r'$r/R_\odot$', fontsize=17)
xlim([0, 1])
ylabel(r'Pressure ($Pa$)', fontsize = 15)
axvline(0.713, color='black', lw=2, linestyle='steps--')
axvline(0.2, color='black', lw=2, linestyle='steps--')
grid(True)

####3
ax3 = subplot(223)
plot(r1, rho1, 'g', lw=3.5)

title('(C)', fontsize=15)
xlim([0, 1])
xlabel(r'$r/R_\odot$', fontsize=17)
ylabel(r'Density ($kgm^-^3$)', fontsize = 15)
axvline(0.713, color='black', lw=2, linestyle='steps--')
axvline(0.2, color='black', lw=2, linestyle='steps--')
ax = gca() 
ax.ticklabel_format(style='sci', axis='y') 
ax.yaxis.major.formatter.set_powerlimits((0,0)) 

grid(True)

####4
ax4 = subplot(224, sharex=ax2)

title('(D)', fontsize=15)
plot(r1, m1, 'c', lw=3.5)
xlim([0, 1])
ylim([0, 1])
xlabel(r'$r/R_\odot$', fontsize=17)
axvline(0.713, color='black', lw=2, linestyle='steps--')
axvline(0.2, color='black', lw=2, linestyle='steps--')
ylabel(r'Mass of Sun ($m/M_\odot$)', fontsize = 15)
grid(True)

show()



#######################################################################################
# PLOT NUMBER 2
# TEMPERATURE PRESSURE DENSITY MASS OF SUN (Convective Zone)

close()
fig = figure(2)

bleh= dict(wspace= 0.3, hspace=0.2)
subplots_adjust(**bleh)

ax1 = subplot(221)
#setp( ax1.get_xticklabels(), visible=False)

title('(E)', fontsize=15)
plot(r1, t1, 'r', lw = 3.5)
xlim([0.7, 1])
ylim([0, 2.5E6])
xlabel(r'$r/R_\odot$', fontsize=17)
ax1.ticklabel_format(style='sci', axis='y') 
ax1.yaxis.major.formatter.set_powerlimits((0,0)) 
ylabel(r'Temperature ($K$)', fontsize = 15)

grid(True)

#####2


ax2 = subplot(222)
#setp( ax2.get_xticklabels(), visible=False)
title('(F)', fontsize=15)
xlabel(r'$r/R_\odot$', fontsize=17)
plot(r1, P1, 'b', lw=3.5)
xlim([0.7, 1])
ylim([0, 0.7E13])
ylabel(r'Pressure ($Pa$)', fontsize = 15)

grid(True)

####3
ax3 = subplot(223)
plot(r1, rho1, 'g', lw=3.5)
title('(G)', fontsize=15)

xlim([0.7, 1])
ylim([0, 2.5E2])
xlabel(r'$r/R_\odot$', fontsize=17)
ylabel(r'Density ($kgm^-^3$)', fontsize = 15)

ax = gca() 
ax.ticklabel_format(style='sci', axis='y') 
ax.yaxis.major.formatter.set_powerlimits((0,0)) 

grid(True)

####4
ax4 = subplot(224, sharex=ax2)
title('(H)', fontsize=15)
plot(r1, m1, 'c', lw=3.5)
xlim([0.7, 1])
ylim([0.97, 1])
xlabel(r'$r/R_\odot$', fontsize=17)

ylabel(r'Mass of Sun ($m/M_\odot$)', fontsize = 15)
grid(True)

show()

#######################################################################################
#PLOT NUMBER 3
# OPACITY OF SUN
close()

semilogy(r1, opacity, 'k', lw=2)
ylim([1E-2, 1E6])
xlim([0, 1.02])
grid()


axvline(0.713, color='black', lw=2, linestyle='steps--')
axvline(0.2, color='black', lw=2, linestyle='steps--')
xlabel(r'$r/R_\odot$', fontsize=18)
ylabel(r'Opacity ($m^2/kg$)', fontsize = 18)
from pylab import *

text(0.05, 316228, 'Core', rotation=00, fontsize=19)
text(0.339, 316228, 'Radiative Zone', rotation=00, fontsize=19)
text(0.77, 316228, 'Convective', rotation=00, fontsize=19)
text(0.82, 90851, 'Zone', rotation=00, fontsize=19)
text(0.65, 1100, 'Tachocline', rotation=90, fontsize=19)

show()



#######################################################################################
#PLOT NUMBER 4
# IDEAL GAS SUN
close()
ax1 = subplot(211)
plot(r1, t1, 'r')
plot(r1, idealtemp(P1), 'k')
ylabel(r'Temperature ($K$)', fontsize = 18)
legend((r'$T_{\odot}$', r'$T_{ideal}$'), loc='upper right')
setp( ax1.get_xticklabels(), visible=False)
xlim([0, 1])
axvline(0.713, color='black', lw=2, linestyle='steps--')
axvline(0.2, color='black', lw=2, linestyle='steps--')
grid()

ax2 = subplot(212)
plot(r1, (t1-idealtemp(P1))/idealtemp(P1)*100)
xlabel(r'$r/R_\odot$', fontsize=18)
ylabel(r'$\left(\frac{T_{\odot}-T_{ideal}}{T_{ideal}}\right)\%$', fontsize=20)
xlim([0, 1])
axvline(0.713, color='black', lw=2, linestyle='steps--')
axvline(0.2, color='black', lw=2, linestyle='steps--')
grid()

show()

 


#######################################################################################
# PLOT NUMBER 5
# PRESSURE SCALE HEIGHT WITH 'ERROR' BARS

xheight = np.linspace(0, 1, 56)
er =[]
for i in arange(0, len(xheight), 1):
	err = interp(xheight[i], r1, pheight(P1, rho1, m1, xx1))
	er.append(err)
yheight = array(er)

yheighter=[]
for i in arange(0, len(xheight), 1):
	if i < 38:
		yter = 0
	else:
		yter = yheight[i]
	yheighter.append(yter)

yheighterr = array(yheighter)

zeroes = zeros(56)
close()
figure(1)

plot(r1, pheight(P1, rho1, m1, xx1,), 'b', lw=2)
plot(r1, pheight(P1, rho1, m1, xx1), 'k', lw=2)
plot(r1, pheight(P1, rho1, m1, xx1), 'r', lw=2)

errorbar(xheight, yheight, xerr=[zeroes, 3*(yheighterr/radiussol)], fmt='b.', lw=2)
errorbar(xheight, yheight, xerr=[zeroes, 2*(yheighterr/radiussol)], fmt='k.', lw=2)
errorbar(xheight, yheight, xerr=[zeroes, (yheighterr/radiussol)], fmt='r.', lw=2)

xlim([0.7, 1.02])
ylim([0, 0.6E8])
grid(True)
xlabel(r'$r/R_\odot$', fontsize=18)
ylabel(r'Pressure Scale Height ($m$)', fontsize=18)

legend((r'$\alpha$ = 3', r'$\alpha$ = 2', r'$\alpha$ = 1'), loc = 'bottom left')


show()


#######################################################################################
# PLOT NUMBER 6
# ADIABATIC/ACTUAL TEMPERATURE GRADIENTS

close()

bleh=dict(wspace= 0.4, hspace=0.2)
plot(r1, abs(tslope), 'k', lw=3)
plot(r1, abs(adbtempgrad1(gamma1, t1, P1)), 'r', lw=2)
#plot(r1, abs(adbtempgrad4(gamma1, t1, P1, rho1)), 'g')

xlim([0.5, 1])
ylim([0, 0.02])

axvline(0.713, color='black', lw=2, linestyle='steps--')
ylabel(r'$|\frac{dT}{dr}|$', fontsize = 22)
xlabel(r'$r/R_\odot$', fontsize=18)
legend((r'$|\frac{dT}{dr}|_{act}$', r'$|\frac{dT}{dr}|_{adb}=(\frac{\gamma-1}{\gamma})(\frac{T}{P})(-g\rho)$',), loc='lower left')
grid(True)

 
show()


#######################################################################################
# PLOT NUMBER 7
# TEMPERATURE DIFFERENCE OF CELL

close()

plot(r1, t1, 'r', lw=2)
plot( r1, tempfinal1, 'k', r1, tempfinal2, 'b')

legend(('Solar Temperature', r'$T_{initial} = 2 \times T_{0.9R_{\odot}}$', r'$T_{initial} = 1.25 \times T_{0.9R_{\odot}}$'), loc='upper center')

xlim([0.9, 1])
ylim([0, 1.2E6])

xlabel(r'$r/R_\odot$', fontsize=18)
ylabel(r'Temperature ($K$)', fontsize = 18)
ylabel
grid()
show()



#######################################################################################
# PLOT NUMBER 8
# CELL RADIUS/VOLUME

#NOTE: Mass of cell needs to be at least 1E14kg to be seen with the xlimits below.
close()
fig = figure(1)
m_cell=1E14
ax1 = subplot(211)
setp( ax1.get_xticklabels(), visible=False)
#title('Volume and Radius of Cell', fontsize = 22)
semilogy(r1, volbub(m_cell, P1, tempfinal1), 'k', r1, volbub(m_cell, P1, tempfinal2), 'b')
xlim([0.9, 1.01])
ylim([1E12, 1E19])
ylabel('Volume of Cell ($m^{3}$)', fontsize = 17)
legend((r'$T_{initial}=2\times T_{0.9R_{\odot}}$', r'$T_{initial}=1.25\times T_{0.9R_{\odot}}$'), loc='upper left')

#legend((r'Cell= $3\times10^6 K$', r'Cell= $2.6\times10^6 K$'), loc='upper center')
grid(True)

ax2 = subplot(212)
semilogy(r1, radiuscell(m_cell, P1, tempfinal1), 'k', r1, radiuscell(m_cell, P1, tempfinal2), 'b')
xlim([0.9, 1.01])
ylim([1E4, 1E6]) 
xlabel(r'$r/R_\odot$', fontsize=18)
ylabel('Radius of Cell ($m$)', fontsize = 17)
legend((r'$T_{initial}=2\times T_{0.9R_{\odot}}$', r'$T_{initial}=1.25\times T_{0.9R_{\odot}}$'), loc='upper left')
#legend((r'Cell= $3\times10^6 K$', r'Cell= $2.6\times10^6 K$'), loc='upper center')
grid(True)


show()

#######################################################################################
# PLOT NUMBER 9
# APPLICABILITY OF 1-D SIM.
#Note: This also depends on the cell mass chosen.

close()
figure(1)
semilogy(r1, Pdifference(xx1, m_cell, P1, tempfinal1), 'r') 
semilogy(r1, rhodifference(xx1, m_cell, P1, tempfinal1), 'k')
semilogy(r1, tdifference(xx1, m_cell, P1, tempfinal1), 'b')
xlim([0.9, 1.001])
ylim([1E-6, 1E0])
xlabel(r'$r/R_\odot$', fontsize=18)
ylabel(r'$\frac{\bigtriangleup x}{x}$', fontsize = 30)
legend(('$x $ = Pressure', '$x $  = Density', '$x $  = Temperature'), loc='upper center')

#axvline(0.9994407736403459, color='black')
grid(True)

show()


###########################
############################################################
# PLOT NUMBER 10
# KINEMATICS OF CELL

close()
F = pylab.gcf()


bleh=dict(wspace= 0.3, hspace=0.3)
subplots_adjust(**bleh)#left=0.0, right=1.0, bottom=0.0, top=1.0)

####1
ax1 = subplot(211)
plot(distarray1/radiussol, accarray1, 'k', distarray2/radiussol, accarray2, 'b', distarray3/radiussol, accarray3, 'g', distarray4/radiussol, accarray4, 'r')
setp( ax1.get_xticklabels(), visible=False)
#title(r'$T_{cell} = 2\times T_{ext}$', fontsize=18)
xlim([0.9, 1.002])
ylabel(r'Acceleration ($ms^{-2}$)', fontsize = 17)

grid(True)

####2
ax2 = subplot(212)
plot(distarray1/radiussol, speedarray1/1000, 'k', distarray2/radiussol, speedarray2/1000, 'b', distarray3/radiussol, speedarray3/1000, 'g', distarray4/radiussol, speedarray4/1000, 'r')
xlim([0.9, 1.002])
ylim([0, 200])
ylabel(r'Velocity ($kms^{-1}$)', fontsize=17)
xlabel(r'$r/R_\odot$', fontsize=18)

legend((r'$T_{initial}=2\times T_{0.9R_{\odot}}$', r'$T_{initial}=1.25\times T_{0.9R_{\odot}}$', r'$T_{initial}=1.1\times T_{0.9R_{\odot}}$', r'$T_{initial}=1.01\times T_{0.9R_{\odot}}$'), loc='upper left')
grid()

show()


#######################################################################################

# PLOT NUMBER 11
# GRAVITY/BUOYANCY ACTING UPON CELL INNIT

close()

figure(1)
plot(r1, netacceleration(m1, xx1, m_cell, P1, tempfinal1, rho1), 'g', r1, g(m1, xx1), 'b--', r1, buoy(m1, xx1, m_cell, P1, tempfinal1, rho1), 'r--') 
grid(True)
xlim([0.9, 1])
ylim([-250, 700])

axhline(0, color='black', lw=2)
xlabel(r'$r/R_\odot$', fontsize=18)
ylabel(r'Acceleration  ($ms^{-2}$)', fontsize = 18) 
legend(('Net Acceleration', 'Gravity', ' Buoyancy'), loc = 'bottom left')

show()

#######################################################################################
# PLOT NUMBER 12
# SCATTER PLOT OF AVERAGE/FINAL VELOCITIES
#This requires a bit of fiddling around. See note on PDF

posits = array(posit)					#This makes the lists arrays.
tempers = array(temper)					# so they can be manipulated
finalspees = array(finalspee)
averages = array(average)

averages=averages/10000					#This changes the average velocity into km/s. The extra factor of 10 is to correct a mistake I had 
tempers=tempers/1E6						# in calculating the average velocity above.
finalspees=finalspees/1000				# Puts final speed in km/s
close()

plt.scatter(posits, tempers, c=averages, marker='s', s=400, cmap = 'hot', edgecolors='none')

xlim([0.695, 0.9495])
ylim([0.5, 2.9])

ylabel(r'Initial Cell Temperature ($10^{6}K$)', fontsize = 19)
xlabel(r'Initial Cell Position $r/R_\odot$', fontsize=19)

cbar = colorbar(ticks=[0, 50, 100, 150, 200, 250, 300, 350], orientation='vertical')

plot(r1, t1/1E6, 'w', lw=3)						#These add multiples of the solar temperature.
plot(r1, (1.25*t1)/1E6, 'w')#, lw=3)
plot(r1, (1.5*t1)/1E6, 'w')#, lw=3)
plot(r1, (1.75*t1)/1E6, 'w')#, lw=3)

show()

"""
