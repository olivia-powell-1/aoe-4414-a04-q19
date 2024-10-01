# ecef_to_eci.py
#
# Usage: python3 eci_to_ecef.py year month day hour minute second eci_x_km eci_y_km eci_z_km
#  converts from the eci frame to the ecef frame
# Parameters:
#  year
#  month
#  day
#  hour
#  minute
#  second
#  ecef_x_km
#  ecef_y_km
#  ecef_z_km
#  ...
# Output:
#  eci_x_km
#  eci_y_km
#  eci_z_km
#
# Written by Olivia Powell
# Other contributors: None
#
# Optional license statement, e.g., See the LICENSE file for the license.

# import Python modules
# e.g., import math # math module
import sys # argv
import math as m
import numpy as np

# "constants"
r_e_km = 6378.137 # redius of earth
e_e = 0.081819221456 # eccentricity of earth
w = 7.292115*10**(-5) # omega in rad/sec

# helper functions

## function description
## calculated denominator
def calc_denom(ecc, lat_rad):
  return m.sqrt(1.0-(ecc**2)*(m.sin(lat_rad)**2))

## calculate jd_frac
def solve_jd_frac(year, month, day, hour, minute, sec):
  jd = day-32075+\
    int(1461*(year+4800+int((month-14)/12))/4)+\
    int(367*(month-2-int((month-14)/12)*12)/12) \
    -int(3*int((year+4900+int((month-14)/12))/100)/4)
  jd_midnight = jd-0.5
  d_frac = (sec+60*(minute+60*hour))/86400
  jd_frac = jd_midnight+d_frac
  return jd_frac

# initialize script arguments
year = int(0)
month = int(0)
day = int(0)
hour = int(0)
minute = int(0)
sec = float('nan')
ecef_x_km = float('nan')
ecef_y_km = float('nan')
ecef_z_km = float('nan')

# parse script arguments
if len(sys.argv)==10:
  year = int(sys.argv[1])
  month = int(sys.argv[2])
  day = int(sys.argv[3])
  hour = int(sys.argv[4])
  minute = int(sys.argv[5])
  sec = float(sys.argv[6])
  ecef_x_km = float(sys.argv[7])
  ecef_y_km = float(sys.argv[8])
  ecef_z_km = float(sys.argv[9])
  ...
else:
  print(\
   'Usage: '\
   'python3 ecef_to_eci.py year month day hour minute second ecef_x_km ecef_y_km ecef_z_km'\
  )
  exit()

# write script below this line
# given the inputs solve for jd_frac
jd_frac = solve_jd_frac(year, month, day, hour, minute, sec)

# convert jd_frac to gdst in sec
tuti = (jd_frac-2451545.0)/36525
gmst = 67310.54841+(876600*60*60+8640184.812866)*tuti+0.093104*tuti**2-(6.2e-6)*tuti**3

# convert gdst in sec to gdst in rad
gmst_rad = (gmst%86400*w+2*m.pi)%(2*m.pi)
                
# rotation matrix on eci vector using angle
rz = np.array([[m.cos(-gmst_rad), m.sin(-gmst_rad), 0], [-m.sin(-gmst_rad), m.cos(-gmst_rad), 0], [0, 0, 1]])
ecef_vec = np.array([[ecef_x_km], [ecef_y_km], [ecef_z_km]])
eci_vec = np.matmul(rz,ecef_vec)
eci_x_km = eci_vec[0]
eci_y_km = eci_vec[1]
eci_z_km = eci_vec[2]

# print
print(eci_x_km)
print(eci_y_km)
print(eci_z_km)