import pandas
import numpy

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
gamma_R_3 = 2 / 3
C_V = 1
C_D = 1
Nr = 1
Nd = 1

# disease free eq
Sr = (Lambda_R/mu_R) * ((omega_R + mu_R)/(omega_R+mu_R+nu_R))

Vr = (Lambda_R/mu_R) * ((nu_R)/(omega_R+mu_R+nu_R))

Rr = (nu_R/(mu_R + omega_R)) * Sr

Er = 0

Ir = 0

Sd = (Lambda_D/mu_D) * ((omega_D + mu_D)/(omega_D + mu_D + nu_D))

Vd = (Lambda_D/mu_D) * ((nu_D)/(omega_D + mu_D + nu_D))

Rd = (nu_D/(mu_D + omega_D)) * Sd

Ed = 0

Id = 0

diff_SR = Lambda_R + omega_R*gamma_R_1*Er + \
    omega_R*Rr - (mu_R + nu_R + beta_RR * (Ir/Nr))*Sr
print('The value of dSR/dt is')
print(diff_SR)

diff_ER = beta_RR * (Ir/Nr) * Sr - (nu_R + mu_R + sigma_R *
                                    gamma_R_3 + sigma_R * gamma_R_1) * Er
print('The value of dER/dt is')
print(diff_ER)

diff_IR = sigma_R*gamma_R_3*Er - (mu_R + delta_R)*Ir
print('The value of dIR/dt is')
print(diff_IR)

diff_RR = Sr * nu_R + Er * nu_R - (omega_R + mu_R) * Rr
print('The value of dRR/dt is')
print(diff_RR)

diff_SD = Lambda_D + omega_D * Rd + sigma_D*gamma_D_1*Ed - (mu_D + beta_RD*(Ir/Nr) + beta_DD*(Id/Nd) + nu_D) * Sd
print('The value of dSD/dt is')
print(diff_SD)

diff_ED = (beta_RD*(Ir/Nr) + beta_DD*(Id/Nd))*Sd - (mu_D + sigma_D * + sigma_D*gamma_D_2 + sigma_D*gamma_D_1) * Ed
print('The value of dED/dt is')
print(diff_ED)

diff_ID = sigma_D*gamma_D_3 - (mu_D + sigma_D) * Id
print('The value of dID/dt is')
print(diff_ID)

diff_RD = sigma_D*gamma_D_3 + nu_D*Sd - (omega_D + mu_D)*Rd
print('The value of dRD/dt is')
print(diff_RD)
