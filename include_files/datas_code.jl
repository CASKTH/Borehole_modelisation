##### COMMAND OF HEAT PUMP ______________________________________________
### TEMP Sources T1, T2 and ground 


#Time span for the simulation
tspan = 86400*4;

# Here the material is water 
# density in kg/m3
ro_fluid = 1000 #pure water density 
c_fluid = 4182

L =1.0 # lenght of the BOREHOLE in meters
## DIMENSIONS OF THE BOREHOLE ###########################
#borehole radius 
r_b = 0.08;
#D value distance between center of pipe and center of BHE 
D = 0.025;
#r p is the outside pipe diameter
# x_c shank spacing (m)
x_c = 0.5*r_b
#internal Pipe radius 
r_pipe = 0.27*r_b

#mass flow rate ########################################
mass_density_fluid=ro_fluid
flow_velocity = 0.01# m/s
volume_flow_rate = round(flow_velocity*(pi*r_pipe^2), digits = 7)
specific_heat_fluid = c_fluid 

m = mass_density_fluid*volume_flow_rate

#EQUATION (13) of [1]
Const_ = m*c_fluid

# ### HEAT TRANSFER COEF COMPUTATION ______________________________________________
# Temp_fluid = source.v
# q = source.i
# Temp_pipe = source_g.v
#convection coefficient (W m−2 K−1 
# h = q/(Temp_fluid- Temp_pipe) # q_fluid/(Tfluid - Tgrout_solid)
h = 100
# #In the
# TRCM, the fluid convection resistance (Rf) and pipe conduction
# resistance (Rp) are lumped together within Rfp ¼ Rf þ Rp
k_f= 0.4 #thermal conductivity of fluid 


#Nu = (h*L)/k_f

### RESISTOR COMPUTATION ________________________________________________________

##THERMAL conductivity ###########################

#Rb borehole resistance 
#thermal conductivity grout
k_g = 2.9; 
#thermal conductivity pipe (W m−1 K−1 )
k_pipe = 1.0;
#thermal conductivity soil 
k_soil = 3 
## useful value used 
teta = (k_g-k_soil)/(k_g+k_soil)



#outside pipe radius (m) same r_pire+ 5%
r_po = round(r_pipe*0.05 +r_pipe, digits = 5)

# for the soil 
r_soil = 10; # 20 metres de diameter


## RESISTOR ###########################
## 
# EQUATION (7) of [1]
#R'b =
Rb = (1/(4*pi*k_g))*(log(r_b/r_pipe) +log(r_b/(2*x_c)) + (teta*log(((r_b/x_c)^4)/(((r_b/x_c)^4)-1))))
Rb = round(Rb, digits = 5)
# conductive EQUATION (5) of [1]
Rcond = round((log(r_po/r_pipe))/(4*pi*k_pipe),digits = 5)

# convective EQUATION (4) of [1]
R_f = 1/(4*pi*r_pipe*h)

# EQUATION (3) for pipe of [1]
R_pipe = R_f + Rcond
#  #The third contribution is a 2D thermal resistance of the grout, and
# is the only one considered in this section. We will assume that R′pipe = 0, but it can be added to the grout resistance in a real com-
# putation.

# internal resistancefor equivalent symetric case 

# EQUATION (28) of [1] 
# Using variables at (9)
R_a = (1/(pi*k_g))*log((2*(r_b/r_pipe)*(((r_b/x_c)^2+1))^teta)/((r_b/x_c)*((r_b/x_c)^2-1)^teta))
R_a = round(R_a, digits = 5)

# EQUATION (9) of [2] 
x_cst = log((sqrt(r_b^2+2*r_po^2)/(2*r_po)))/log(r_b/(r_po*sqrt(2)))

# EQUATION (3) of [1]
#R1 = R2 = 2*Rbt
#Rb = Rbt = Rpipe + R'b 
# Rbt= Rb + R_pipe
R_g = (2*Rb-R_pipe)
R_GB = round((1-x_cst)*R_g, digits = 5)

# EQUATION (18) of [2] 
#  direct coupling resistance
# R_12 = round((2*R_GB*((R_a-R_pipe*2)-2*x_cst*R_g)/(2*R_GB-(R_a-2*R_pipe)+2*x_cst*R_g)), digits = 5)
R_12 = round((2*R_g*(R_a-2*(R_f+Rcond))/(2*R_g-(R_a-2*(R_f+Rcond)))))

# EQUATION inspired of (5) of [1] 
# r_b BHE radius

r_soil_1 = 0.5;
r_soil_2 = 2.5;
r_soil_3 = 10;

Rsoil_1 = round((log(r_soil_1/r_po))/(4*pi*k_soil), digits = 5)
Rsoil_2 = round((log(r_soil_2/r_soil_1))/(4*pi*k_soil), digits = 5)
Rsoil_3 = round((log(r_soil_3/r_soil_2))/(4*pi*k_soil), digits = 5)

### VOLUME COMPUTATION ________________________________________________________
# Here the material is Soil 
# density of the Silty Sand and Gravel  put to 1442 kg/m3
ro_soil = 1442
c_soil = 800 #soil dry # FOR 


# Here the material is PIPE 
# density of the PIPE put to 1600 kg/m3
ro_pipe = 1600
c_pipe = 2185

# Here the material is GROUT 
# density of the grout put to 1600 kg/m3
rog = 1600
#specific heat 
cg = 800#2000 

#L = 1 meter deep BHE and soil 
# #Vg ¼ pr2
# pi*rb² - 2*Vf - 2Vp - Vgg
# Vgg = 4 *r_po*D - pi*r_po²
# EQUATIONS (13) (14) of [2]
# according to the exemple of the MTRCM with nf = np = 2 and ng = ngg ) ns = 3
np = 2; # pour Rcond
nf = np;# pour R_f

ng = 3
ngg = ng
ns = ng


Vf = round(π*r_pipe^2*L/nf, digits = 5)
Vp = round((π*r_po^2*L - Vf)/np, digits = 5)
Vgg = round(L*(4*r_po*D - pi*r_po^2)/ngg , digits = 5)
Vg = round(L*(pi*r_b^2 - 2*Vf - 2Vp - Vgg)/ng, digits = 5)

Vsoil_1 = round(((pi*r_soil_1^2)*L - (pi*r_po^2)*L), digits = 5)
Vsoil_2 = round(((pi*r_soil_2^2)*L - (pi*r_soil_1^2)*L), digits = 5)
Vsoil_3 = round(((pi*r_soil_3^2)*L - (pi*r_soil_2^2)*L), digits = 5)
#C= rog cg Vg/2
# rog, cg  and vg are respectively the density, the specific heat and
# the volume of the material corresponding to the n j resistances
# surrounding each node

### CAPACITOR COMPUTATION ________________________________________________________
#EQUATION (11) of [2]
C_grout = round(rog*cg*Vg/2, digits = 3)
C_gg = rog*cg*Vgg
C_pipe = ro_pipe*c_pipe*Vp/2
C_fluid = ro_fluid*c_fluid*Vf/2

C_soil_1 = round(ro_soil*c_soil*Vsoil_1/2, digits = 5)
C_soil_2 = round(ro_soil*c_soil*Vsoil_2/2, digits = 5)
C_soil_3 = round(ro_soil*c_soil*Vsoil_3/2, digits = 5)



### LAYOUT CONSTRUCTION ___________________________________________________________
#EQUATION (12) and the explaination gived at 2.2 of [2]

Rf = round(R_f/nf, digits = 5)
Rp = Rcond/np

Rgg = R_12/ngg
Rg = R_g/ng
# Rs = Rsoil/ns
#capacitor divided


#joint capacitance 
C_g_p_gg = C_gg/2 + C_grout/2 + C_pipe/2 
# C_g_p_gg = C_grout
C_g_s = C_grout + C_soil_1
# C_g_s = C_soil_1 + C_grout/2