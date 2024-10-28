
import numpy as np
import argparse
import filtools


parser = argparse.ArgumentParser(description='Calculate max. acceleration from porb, pulsar, companion mass assuming circular orbit')
parser.add_argument('-P', '--porb', help='Porb in seconds', default=3600, type=float)
parser.add_argument('-c', '--comp_mass', help='Companion mass in solar mass units',
                    type=float, default=1.0)
parser.add_argument('-p', '--pulsar_mass', help='Pulsar mass in solar mass units',
                    type=float, default=1.4)
args = parser.parse_args()

# input = filtools.FilterbankIO()
# input.read_header(args.fil)

# T0 = input.header['tstart']
angular_velocity  = 2 * np.pi/args.porb

companion_mass = args.comp_mass
pulsar_mass = args.pulsar_mass
mass_function = companion_mass**3/(pulsar_mass + pulsar_mass)**2
T0 =  4.925490947e-06 # Sun Time constant
speed_of_light = 2.99792458e+08
a_max = angular_velocity**(4/3) * (T0 * mass_function)**(1/3) * speed_of_light



print(f'{a_max}') 