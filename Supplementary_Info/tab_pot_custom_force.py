import numpy as np
import math

def WCA_energy(b, r):
# Calculate WCA energy
  E_pot = 0
  val1 = math.pow((b / r), 12)
  val2 = -math.pow((b / r), 6)
  E_pot = 4 * epsilon * (val1 + val2 + 0.25)
  return E_pot

def WCA_forces(b, r):
# Calculate WCA forces
  Force = 0
  val1 = 24 * math.pow(b, 6) / math.pow(r, 7)
  val2 = -48 * math.pow(b, 12) / math.pow(r, 13)
  Force = -epsilon*(val1 + val2)
  return Force

def cus_energy(x, epsilon=20, sigma=1, n=12, m=6, bump_height=3, bump_center=1.2, bump_width=0.1):
    # Lennard-Jones-like potential with an additional Gaussian bump
    lj_potential = epsilon * ((sigma / x)**n - (sigma / x)**m)
    gaussian_bump = bump_height * np.exp(-((x - bump_center)**2 / (2 * bump_width**2)))
    potential = lj_potential + gaussian_bump
    return potential

def cus_force(x):
    A = 3  # Example values, replace with actual values
    B = 1.2
    C = 0.1
    E = np.exp(1)  # Base of the natural logarithm
    sigma = 1  # Replace with actual value
    epsilon = 20  # Replace with actual value
    term1 = -((A * (-B + x)) / (C**2 * E**((-B + x)**2 / (2 * C**2))))
    term2 = ((6 * sigma**6) / x**7 - (12 * sigma**12) / x**13) * epsilon# * (- (sigma**6 / x**6) + sigma**12 / x**12)
    return -(term1 + term2)

def Tail_energy(b, r):
# Calculate extra Tail energy
  E_pot = 0
  if (r < r_cut):
    E_pot = -1 * epsilon
  else:
    val1 = math.cos((math.pi * (r - r_cut)) / (2 * wca))
    E_pot = -1 * epsilon * math.pow(val1, 2)
  return E_pot

def Tail_forces(b, r):
  Force = 0
  if (r < r_cut):
    Force = 0;
  else:
    val1 = math.sin((math.pi * (r - r_cut)) / wca)
    Force = - epsilon * math.pi * val1 / (2 * wca)
  return Force

def tabs(w_cut,eps):
    global sigma
    global epsilon
    global r_cut
    global wca
    wca = 1.5
    epsilon = 1
    sigma = 1
    r_space = 0.01
    r_max = 2.5 # wca + 2**(1/6)+ r_space
    r_cut = 2**(1/6)
    rs = np.arange(r_space,r_max,r_space)
    with open( 'tabulated_potential.dat', 'w' ) as g:
        g.write('# Tabulated potential for Cooke 3-bead lipid model, Wc = 1.5\n\n')
        g.write('HEAD_HEAD\n')
        g.write('N {0:.0f} R 0.000001 {1:.6f}\n\n'.format(len(rs)+1,r_max))
        g.write('1 0.000001 0.000000 0.000000\n')
        for i in range(len(rs)):
            if rs[i] <= r_cut*0.95:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], WCA_energy(0.95,rs[i]), WCA_forces(0.95,rs[i]) ) )
            else:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], 0, 0 ) )
        g.write('\n')
        g.write('HEAD_TAIL\n')
        g.write('N {0:.0f} R 0.000001 {1:.6f}\n\n'.format(len(rs)+1,r_max))
        g.write('1 0.000001 0.000000 0.000000\n')
        for i in range(len(rs)):
            if rs[i] <= r_cut*0.95:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], WCA_energy(0.95,rs[i]), WCA_forces(0.95,rs[i]) ) )
            else:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], 0, 0 ) )
        g.write('\n')
        g.write('TAIL_TAIL\n')
        g.write('N {0:.0f} R 0.000001 {1:.6f}\n\n'.format(len(rs)+1,r_max))
        g.write('1 0.000001 0.000000 0.000000\n')
        for i in range(len(rs)):
            if rs[i] < r_cut:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], WCA_energy(sigma,rs[i])+Tail_energy(sigma, rs[i]), WCA_forces(sigma,rs[i])+Tail_forces(sigma, rs[i]) ) )
            elif rs[i] <= r_cut+wca:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], Tail_energy(sigma, rs[i]), Tail_forces(sigma, rs[i]) ) )
            else:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], 0, 0 ) )
        g.write('\n')
        g.write('TEST_HEAD\n')
        g.write('N {0:.0f} R 0.000001 {1:.6f}\n\n'.format(len(rs)+1,r_max))
        g.write('1 0.000001 0.000000 0.000000\n')
        for i in range(len(rs)):
            if rs[i] <= r_cut*0.95:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], WCA_energy(1.1,rs[i]), WCA_forces(1.1,rs[i]) ) )
            else:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], 0, 0 ) )
        g.write('\n')
        g.write('PART_HEAD\n')
        g.write('N {0:.0f} R 0.000001 {1:.6f}\n\n'.format(len(rs)+1,r_max))
        g.write('1 0.000001 0.000000 0.000000\n')
        wca *= 0.3
        for i in range(len(rs)):
            if rs[i] < 2.5:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], cus_energy( rs[i]), cus_force(rs[i]) ) )
            else:
                g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], 0, 0 ) )
        g.write('\n')
        g.write('ZERO_ZERO\n')
        g.write('N {0:.0f} R 0.000001 {1:.6f}\n\n'.format(len(rs)+1,r_max))
        g.write('1 0.000001 0.000000 0.000000\n')
        for i in range(len(rs)):
            g.write("{0:.0f} {1:.6f} {2:.6f} {3:.6f}\n".format( i+2, rs[i], 0, 0 ) )
        g.write('\n')
        g.close()
