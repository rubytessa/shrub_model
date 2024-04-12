## OTC calculations 
library(tidyverse)

# inputs - set the large dimensions
x = 120 # long width
z = 60 # height
phi = 60*2*pi/360 # panel angle 


# derived dimensions
a = x/2
y = sqrt(3)*a
theta = atan(y/(a*cos(phi)))
d = z/sin(phi)
c = d/sin(theta)
e = a*z/y * (cos(theta)/sin(theta))
b = a - 2*e


# print
x
y
z
phi

a
b
c
d
theta*360/(2*pi)

## dimensions to cut
long = a+4
short = b+4
height = d

long/2.54
short/2.54
height/2.54

#3 area
print(str_c((long + short + e)*60, " cm2"))

# 7080 cm2 per 

## cm to inches: 1 inch = 2.54 cm

### Set Panel Dimensions ----- 

a = 55
phi = 60*2*pi/360 # panel angle 
d = 50

x = 2*a
y = sqrt(3)/2*x
z = d*sin(theta)
theta


e = d/tan(theta)
b = a-2*e
c = d/sin(theta)

(a_in = a/2.54)
(b_in = b/2.54)
(d_in = d/2.54)




