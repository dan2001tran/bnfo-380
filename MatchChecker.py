import pandas
import numpy
import sympy
import math



Lambda_D = 0.114
Lambda_R = 1.34
mu_D = 0.08
mu_R = 0.836
delta_D = 1
delta_R = 66.36
omega_D = float(1)
omega_R = 0
nu_D = 0.09
nu_R = 0.7
beta_RR = 33.25
beta_RD = 1
beta_DD = 1.58 * (10**-7)
sigma_D = 1/6
sigma_R = 0.13
gamma_D_1 = 1 / 3
gamma_D_2 = 1 / 3
gamma_D_3 = 1 / 3
gamma_R_1 = 1 / 3
gamma_R_2 = 1/3
gamma_R_3 = 1 / 3
C_V = 1
C_D = 1

d_er = mu_R + nu_R + sigma_R * gamma_R_1 + sigma_R*gamma_R_3

d_ed = mu_D + sigma_D*(gamma_R_1 + gamma_R_2 + gamma_R_3)

SpecialR_D = (beta_DD * sigma_D * gamma_D_3) / (d_ed * (mu_D + delta_D))

SpecialR_R = (beta_RR * sigma_R * gamma_R_3) / (d_er * (mu_R + delta_R))

squareRoot = []
#equation to solve a quadratic equation

def equationroots( a, b, c): 
  
    # calculating discriminant using formula
    dis = b * b - 4 * a * c 
    sqrt_val = math.sqrt(abs(dis)) 
      
    # checking condition for discriminant
    if dis > 0: 
        print(" real and different roots ") 
        print((-b + sqrt_val)/(2 * a)) 
        print((-b - sqrt_val)/(2 * a)) 

        root1 = (-b + sqrt_val)/(2 * a)
        root2 = (-b - sqrt_val)/(2 * a)
        if root1 > 0:
            squareRoot.append((-b + sqrt_val)/(2 * a))
        else:
            squareRoot.append((-b - sqrt_val)/(2 * a))
    elif dis == 0: 
        print(" real and same roots") 
        print(-b / (2 * a)) 

        squareRoot.append(-b / (2 * a))
    # when discriminant is less than 0
    else:
        print("Complex Roots") 
        print(- b / (2 * a), " + i", sqrt_val) 
        print(- b / (2 * a), " - i", sqrt_val) 
  

# disease free eq
Sr = (Lambda_R/mu_R) * ((omega_R + mu_R)/(omega_R+mu_R+nu_R))

Rr = (Lambda_R/mu_R) * ((nu_R)/(omega_R+mu_R+nu_R))

Er = 0

Ir = 0

Nr = Sr + Rr + Er + Ir

Sd = (Lambda_D/mu_D) * ((omega_D + mu_D)/(omega_D + mu_D + nu_D))

Rd = (Lambda_D/mu_D) * ((nu_D)/(omega_D + mu_D + nu_D))

Ed = 0

Id = 0

#
Nd = Sd + Rd + Ed + Id

print('Solving the disease free equilibrium')


diff_SR0 = Lambda_R + omega_R*gamma_R_1*Er + omega_R*Rr - (mu_R + nu_R + beta_RR * (Ir/Nr))*Sr
print('The value of dSR/dt is')
print(diff_SR0)

diff_ER0 = beta_RR * (Ir/Nr) * Sr - (nu_R + mu_R + sigma_R *
                                    gamma_R_3 + sigma_R * gamma_R_1) * Er
print('The value of dER/dt is')
print(diff_ER0)

diff_IR0 = sigma_R*gamma_R_3*Er - (mu_R + delta_R)*Ir
print('The value of dIR/dt is')
print(diff_IR0)

diff_RR0 = Sr * nu_R + Er * nu_R - (omega_R + mu_R) * Rr
print('The value of dRR/dt is')
print(diff_RR0)

diff_SD0 = Lambda_D + omega_D * Rd + sigma_D*gamma_D_1*Ed - (mu_D + beta_RD*(Ir/Nr) + beta_DD*(Id/Nd) + nu_D) * Sd
print('The value of dSD/dt is')
print(diff_SD0)

diff_ED0 = (beta_RD*(Ir/Nr) + beta_DD*(Id/Nd))*Sd - (mu_D + sigma_D * + sigma_D*gamma_D_2 + sigma_D*gamma_D_1) * Ed
print('The value of dED/dt is')
print(diff_ED0)

diff_ID0 = sigma_D*gamma_D_3*Ed - (mu_D + sigma_D) * Id
print('The value of dID/dt is')
print(diff_ID0)

diff_RD0 = sigma_D*gamma_D_2*Ed + nu_D*Sd - (omega_D + mu_D)*Rd
print('The value of dRD/dt is')
print(diff_RD0)


#Endemic Solution checker for where it is endemic in only dogs
def EndemicSoutionChecker(Id):


    Er = 0

    Ir = 0

    Sr_1 = (Lambda_R/mu_R) * ((omega_R + mu_R)/(omega_R+mu_R+nu_R))

    Vr_1 = (Lambda_R/mu_R) * ((nu_R)/(omega_R+mu_R+nu_R))

    Lambda_D = mu_D * Nd + sigma_D*Id

    Nd = (Lambda_D - sigma_D*Id) / mu_D

    Ed = ((mu_D + sigma_D) * Id) / sigma_D* gamma_D_3

    Sd = (d_ed * (mu_D + sigma_D) * Nd) / (beta_DD * sigma_D * gamma_D_3)

    Rd = (sigma_D * gamma_D_2 * Ed + nu_D * Sd) / (omega_D + mu_D)

    print('solving the endemic dog population equations')

    diff_SR0 = Lambda_R + omega_R*gamma_R_1*Er + omega_R*Rr - (mu_R + nu_R + beta_RR * (Ir/Nr))*Sr
    print('The value of dSR/dt is')
    print(diff_SR0)


    diff_RR0 = Sr * nu_R + Er * nu_R - (omega_R + mu_R) * Rr
    print('The value of dRR/dt is')
    print(diff_RR0)

    diff_SD0 = Lambda_D + omega_D * Rd + sigma_D*gamma_D_1*Ed - (mu_D + beta_RD*(Ir/Nr) + beta_DD*(Id/Nd) + nu_D) * Sd
    print('The value of dSD/dt is')
    print(diff_SD0)

    diff_ED0 = (beta_RD*(Ir/Nr) + beta_DD*(Id/Nd))*Sd - (mu_D + sigma_D * + sigma_D*gamma_D_2 + sigma_D*gamma_D_1) * Ed
    print('The value of dED/dt is')
    print(diff_ED0)

    diff_ID0 = sigma_D*gamma_D_3*Ed - (mu_D + sigma_D) * Id
    print('The value of dID/dt is')
    print(diff_ID0)

    diff_RD0 = sigma_D*gamma_D_3*Ed + nu_D*Sd - (omega_D + mu_D)*Rd
    print('The value of dRD/dt is')
    print(diff_RD0)


#solving id in order to plug it into our equations

a = (nu_D * d_ed * sigma_D * (mu_D + sigma_D))/ (beta_DD * sigma_D * gamma_D_3 * mu_D * (omega_D + mu_D))

b = -1 - (sigma_D/mu_D) - ((mu_D+sigma_D)/sigma_D*gamma_D_3) - (d_ed*sigma_D*(mu_D+sigma_D))/(mu_D*beta_DD*sigma_D*gamma_D_3) - (nu_D*d_ed*Lambda_D*(mu_D+sigma_D))/(beta_DD*sigma_D*gamma_D_3*mu_D*(omega_D+mu_D)) - (gamma_D_2*(mu_D+sigma_D))/mu_D*(omega_D+mu_D)

c = (Lambda_D*(beta_DD*sigma_D*gamma_D_3 - d_ed*mu_D - d_ed*sigma_D))/(beta_DD*sigma_D*gamma_D_3*mu_D)

equationroots(a, b, c)

print('id is', end=' ')
print(squareRoot[0])

EndemicSoutionChecker(squareRoot[0])

