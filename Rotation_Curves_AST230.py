import numpy as np
import matplotlib.pyplot as plt
from astropy.io import ascii
from astropy import units as u
from astropy.constants import G

plt.style.use('ggplot')

def v_circ(m, r):
    return np.sqrt(G*m/r)
def V_sphere(r):
    return (4/3)*np.pi*r**3
def m_radius(density, V_sphere):
    return density*V_sphere
def g(c):
    return 1/(np.log(1+c)-c/(1+c))
def m_NFW(c, s, m_virial):
    return g(c)*m_virial*(np.log(1+c*s)-c*s/(1+c*s))
def v_hern(m_hern, r, a):
    return np.sqrt(r*G*m_hern/((r+a)**2))
def v_MN(m, r, a):
    return np.sqrt(G*m*r**2/(r**2+a**2)**(3/2))
def gra_dev(m, r, a):
    return G*m*r/(r**2+a**2)**(3/2)
def v_triMN(r ,dev1, dev2, dev3):
    return np.sqrt(r*(dev1+dev2+dev3))
def v_total(Mvir,Mhern,c,a,M1,M2,M3,a1,a2,a3,r):
    crit_density= (3*(67.4*u.km/u.s/u.Mpc)**2/(8*np.pi*G)).to(u.kg/(u.m**3))
    r_virial= np.cbrt(3*Mvir/(800*np.pi*crit_density)).to(u.kpc)
    return np.sqrt(v_triMN(r, gra_dev(M1, r, a1), gra_dev(M2, r, a2), gra_dev(M3, r, a3))**2+v_hern(Mhern,r,a)**2+v_circ(m_NFW(c, r/r_virial, Mvir), r)**2)


#2. a
m=10**(12)*u.Msun
r= range(1, 101)*u.kpc
plt.plot(r, v_circ(m, r).to(u.km/u.s), c='#444444')
plt.title('Radius-velocity relation')
plt.xlabel('Radio [kpc]')
plt.ylabel('Velocidad circular [km/s]')
plt.savefig('2.a.png', dpi=300)
plt.show()

#2. b
volume= V_sphere(100*u.kpc)
density= m/volume
r= range(1, 100001)*u.pc
plt.plot(r.to(u.kpc), v_circ(m_radius(density, V_sphere(r)).to(u.Msun), r).to(u.km/u.s), c='#444444')
plt.title('Radius-velocity relation with uniform density')
plt.xlabel('Radio [kpc]')
plt.ylabel('Velocidad circular [km/s]')
plt.savefig('2.b.png', dpi=300)
plt.show()

#2. c
crit_density= (3*(67.4*u.km/u.s/u.Mpc)**2/(8*np.pi*G)).to(u.kg/(u.m**3))
print('Current critical density: ', crit_density)

m_virial= 10**12*u.Msun
r_virial= np.cbrt(3*m_virial/(800*np.pi*crit_density)).to(u.kpc)
print('Radio de Virial: ', r_virial)

r= range(1, 101)*u.kpc
plt.plot(r, v_circ(m_NFW(5, r/r_virial, m_virial), r).to(u.km/u.s), c='#444444', label='c=5')
plt.plot(r, v_circ(m_NFW(15, r/r_virial, m_virial), r).to(u.km/u.s), c='#B06576', label='c=15')
plt.title('Radius-velocity relation of NFW halo')
plt.xlabel('Radio [kpc]')
plt.ylabel('Velocidad circular [km/s]')
plt.legend(loc='best')
plt.savefig('2.c.png', dpi=300)
plt.show()

#3. b
r= range(1, 100001)*u.pc
m_hern= 10**10*u.Msun
plt.plot(r.to(u.kpc), v_hern(m_hern, r, 0.5*u.kpc).to(u.km/u.s), c='#444444')
plt.title('Radius-velocity relation with Hernquits distribution')
plt.xlabel('Radio [kpc]')
plt.ylabel('Velocidad circular [km/s]')
plt.savefig('3.b.png', dpi=300)
plt.show()

#4. b
r= range(1, 100001)*u.pc
m= 3*10**10*u.Msun
plt.plot(r.to(u.kpc), v_MN(m, r, 0.5*u.kpc).to(u.km/u.s), c='#444444')
plt.title('Radius-velocity relation with Miyamoto-Nagai distribution')
plt.xlabel('Radio [kpc]')
plt.ylabel('Velocidad circular [km/s]')
plt.savefig('4.b.png', dpi=300)
plt.show()

#4. c
r= range(1, 100001)*u.pc
m1= 1.9487*3*10**10*u.Msun
m2= -1.3077*3*10**10*u.Msun
m3= 0.2242*3*10**10*u.Msun
a1= 2.0074*3*u.kpc
a2= 4.4441*3*u.kpc
a3= 0.7333*3*u.kpc
plt.plot(r.to(u.kpc), v_triMN(r, gra_dev(m1, r, a1), gra_dev(m2, r, a2), gra_dev(m3, r, a3)).to(u.km/u.s), c='#444444')
plt.title('Radius-velocity relation using a triple-MN disk')
plt.xlabel('Radio [kpc]')
plt.ylabel('Velocidad circular [km/s]')
plt.savefig('4.c.png', dpi=300)
plt.show()

#5.
data= ascii.read('Data.dat')

x=data['R[kpc]']
y=data['Vcirc[km/s]']
error=data['Err_Vcirc[km/s]']

r=range(1, 200001)*u.pc
Mvir=5.5*10**11*u.Msun
c=25
Mhern=2.3*10**11*u.Msun
a=9*u.kpc
M1=8*10**10*u.Msun
M2=-5*10**10*u.Msun
M3=4*10**10*u.Msun
a1=10*u.kpc
a2=15*u.kpc
a3=30*u.kpc
crit_density= (3*(67.4*u.km/u.s/u.Mpc)**2/(8*np.pi*G)).to(u.kg/(u.m**3))
r_virial= np.cbrt(3*Mvir/(800*np.pi*crit_density)).to(u.kpc)

plt.errorbar(x, y, yerr=error, fmt='o', c='#444444', capsize=3, markersize=2, label='Data')
plt.plot(r.to(u.kpc),v_total(Mvir,Mhern,c,a,M1,M2,M3,a1,a2,a3,r).to(u.km/u.s),c='#DF5047', label='V$_{circ,tot}$')
plt.plot(r.to(u.kpc), v_circ(m_NFW(c, r/r_virial, Mvir), r).to(u.km/u.s), c='#65AFB0', label='NFW')
plt.plot(r.to(u.kpc), v_hern(Mhern, r, a).to(u.km/u.s), c='#5933B9', label='Hernquist')
plt.plot(r.to(u.kpc), v_triMN(r, gra_dev(M1, r, a1), gra_dev(M2, r, a2), gra_dev(M3, r, a3)).to(u.km/u.s), c='#5FB321',label='Triple-MN')
plt.title('Radius-velocity relation of data with different fits')
plt.xlabel('Radio [kpc]')
plt.ylabel('Velocidad circular [km/s]')
plt.legend(loc='best')
plt.savefig('5.png', dpi=300)
plt.show()
