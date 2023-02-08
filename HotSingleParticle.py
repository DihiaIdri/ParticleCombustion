import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt


class HotSingleParticle:

    def __init__(self, constM, constG):
        self.C0 = constM.C0  # Oxygen concentration in gas
        self.k0 = constM.k0  # Arrhenius reaction rate prefactor
        self.Ta = constM.Ta  # Activation temperature = Activation energy/universal gas constant
        self.D = constM.D  # Diffusivity m/s
        self.Sh = constM.Sh  # Sherwood non dimensional number
        self.q = constM.q
        self.stoich = constM.stoich
        self.Cond = constM.Cond  # Conductivity W/m-K
        self.Nu = constM.Nu
        self.gas = constG
        self.Tg = constG.Tg
        print(self.Tg)

    def reaction_rate(self):
        rMin = 1  # micro meter, 10**-6
        rMax = 25
        r, Tg = HotSingleParticle.finding_particle_temperature(self, rMin, rMax, self.D)
        plt.plot(r, Tg)
        plt.show()



    # def finding_heat_of_reaction(self):

    def finding_particle_temperature(self, rMin, rMax, D):
        r = []
        Temp = []
        for radius in np.arange(rMin, rMax, 2):
            radius = radius*10**(-6)
            # D = SingleParticle.diffusivity_O2(self)
            b = HotSingleParticle.mass_transfer_coefficient(self, D, radius)  # m/s
            h = HotSingleParticle.convective_HT_coefficient(self, radius)  # W/m2-K
            T0 = np.array([1000])  # guess
            # Add the [0] at the end else Tp comes in array form
            Ts = fsolve(HotSingleParticle.finding_Ts, T0, args=(self, b, h), xtol=1e-6)
            Temp.append(Ts)
            r.append(radius)
            print(Ts)
            print(radius)
        return r, Temp

    @staticmethod
    def finding_Ts(T, self, b, h):
        # Note that the way we write the equation to solve influences the final answer.
        const = self.q*self.stoich*b*self.C0/h
        a = T - const*(self.k0*np.exp(-self.Ta/T))/(self.k0*np.exp(-self.Ta/T) + b)
        c = self.Tg
        return a - c

    def heat_release_rate(self, radius):
        Qr = []
        Temp = []
        for T in range(500, 7000, 250):  # increment by 250
            b = self.Sh*self.D/(2*radius)
            k = self.k0*np.exp(-self.Ta/T)
            Da = k/(k + b)
            Qr.append(self.q*self.stoich*b*Da*self.C0)
            Temp.append(T)
        return Qr, Temp

    def diffusivity_O2(self):  # DO2
        VA = self.gas.VI  # diffusion volumes
        MA = self.gas.MI
        VB = self.gas.VO
        MB = self.gas.MO
        P = self.gas.P
        T = self.gas.Tg

        # Fuller-Schettler-Giddins
        numerator = pow(10, -3)*pow(T, 7/4)*(1/MA + 1/MB)**(1/2)
        denominator = (P/101325)*(pow(VA, 1/3) + pow(VB, 1/3))**2
        # Gilliland method
        # numerator = 435.7*pow(T, 3/2)*(1/MA + 1/MB)**(1/2)
        # denominator = (P*(pow(VA, 1/3) + pow(VB, 1/3))**2)
        return (numerator/denominator)/pow(100, 2)  # to convert cm2/s into m2/s

    def mass_transfer_coefficient(self, D, radius):
        return self.Sh*D/(2*radius)  # m/s

    def convective_HT_coefficient(self, radius):
        return self.Nu*self.Cond/(2*radius)  # W/m2-K

    def Arrhenius_rate_law(self, T):
        # Where T is the particle surface temperature
        return self.k0*np.exp(-self.Ta/T)

    def surface_concentration(self, k, b):  # Cs
        return self.C0*(b/(k+b))

    def normalized_Damkoler_number(self, k, b):
        return k/(k+b)

    # Moving particles
    # def sherwood_number(self):
    # def Nusselt_number(self):
