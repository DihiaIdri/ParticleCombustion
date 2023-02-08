import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

class SingleParticle:

    def __init__(self, constM, constG, globalReact):
        self.R = 8.314510  #Universal gas constant J/mol-K
        self.gas = constG
        self.GR = globalReact
        self.Metal = constM

        self.C0, self.CI = SingleParticle.oxygen_concentration(self)  # kg/m3, Oxygen concentration in gas
        # Diffusivity and concentration  depend on tempertaure and pressure.
        # Here I assume standard state, which is a mistake. because I am considering different gas tempertaure.
        self.D = SingleParticle.diffusivity_O2_FSG(self, self.gas.Tg)  # Diffusivity m/s
        self.stoich = SingleParticle.stoichiometric_coefficient(self)
        kMixture = SingleParticle.conductivity_gas_mixture(self)
        print(kMixture)

        self.k0 = constM.k0  # Arrhenius reaction rate prefactor
        self.Ta = constM.Ta  # Activation temperature = Activation energy/universal gas constant
        self.Sh = constM.Sh  # Sherwood non dimensional number
        self.Nu = constM.Nu
        self.q = constM.q
        self.Cond = constM.Cond  # Conductivity W/m-K

    def reaction_rate(self):
        radius1 = 30*10**(-6)
        radius2 = 20*10**(-6)
        radius3 = 10*10**(-6)
        Qr1, T = SingleParticle.heat_release_rate(self, radius1)
        Qr2, T = SingleParticle.heat_release_rate(self, radius2)
        Qr3, T = SingleParticle.heat_release_rate(self, radius3)

        plt.plot(T, Qr1)
        plt.plot(T, Qr2)
        plt.plot(T, Qr3)
        plt.show()

        rMin = 2  #micro meter, 10**-6
        rMax = 15
        # Diffusivity is a function of the fgas temperature
        r, Tg = SingleParticle.finding_ignition_temperature(self, rMin, rMax, self.D)
        plt.plot(r, Tg)
        plt.show()

        # Study of diffusivity
        Diffusivity = []
        T = []
        for Tgas in range(300, 2000, 200):
            Dchange = SingleParticle.diffusivity_O2_FSG(self, Tgas)
            Diffusivity.append(Dchange)
            T.append(Tgas)
        plt.plot(T, Diffusivity)
        plt.show()

    # def finding_heat_of_reaction(self):

    def finding_ignition_temperature(self, rMin, rMax, D):
        r = []
        Tg = []
        for radius in np.arange(rMin, rMax, 0.1):
            radius = radius*10**(-6)
            # D = SingleParticle.diffusivity_O2(self)
            #A FUNCTION OF Tgas!!!!!!!
            b = SingleParticle.mass_transfer_coefficient(self, D, radius)  # m/s
            h = SingleParticle.convective_HT_coefficient(self, radius)  # W/m2-K
            T0 = np.array([2000])  # guess
            # Add the [0] at the end else Tp comes in array form
            Tp = fsolve(SingleParticle.finding_T, T0, args=(self, b, h), xtol=1e-6)[0]

            k = SingleParticle.Arrhenius_rate_law(self, Tp)
            Da = SingleParticle.normalized_Damkoler_number(self, k, b)
            Qr = self.q*self.stoich*b*Da*self.C0
            Tgas = Tp - (Qr/h)
            Tg.append(Tgas)
            r.append(radius)

        return r, Tg


    @staticmethod
    def finding_T(T, self, b, h):
        # Note that the way we write the equation to solve influences the final answer.
        a = (T*(self.k0*np.exp(-self.Ta/T) + b))**2/(self.k0*np.exp(-self.Ta/T))
        c = self.q*self.stoich*self.C0*self.Ta*(b**2)/h
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

    def conductivity_gas_mixture(self):
        a = self.GR.FO*self.gas.kO + self.GR.FI*self.gas.kI
        b = self.GR.FO/self.gas.kO + self.GR.FI/self.gas.kI
        print(a, 1/b)
        return (1/2)*(a + 1/b)

    def diffusivity_O2_FSG(self, Tgas):  # DO2
        VA = self.gas.VI  # diffusion volumes
        MA = self.gas.MI
        VB = self.gas.VO
        MB = self.gas.MO
        P = self.gas.P

        # Fuller-Schettler-Giddins
        numerator = pow(10, -3)*pow(Tgas, 7/4)*(1/MA + 1/MB)**(1/2)
        denominator = (P/101325)*(pow(VA, 1/3) + pow(VB, 1/3))**2
        return (numerator/denominator)/pow(100, 2)  # to convert cm2/s into m2/s

    def oxygen_concentration(self):
        # Need Global Reaction sheet and the gas constant data.
        FO = self.GR.FO  # mole fraction of oxidizer
        FI = self.GR.FI  # mole fraction of inert
        # number of moles
        nO = 1
        nI = nO*FI/FO
        # finding concentration of each component
        V = (nO + nI)*self.R*self.gas.Tg/self.gas.P  # m3
        cO = nO/(V*self.gas.MO)  #kg/m3
        cI = nI/(V*self.gas.MI)
        # print(nO, nI, V, cO, cI)
        return cO, cI

    def stoichiometric_coefficient(self):
        # Ratio of the reactant mass to oxidizer mass
        # Multiply the Global Reaction (GR) number of moles.
        massO = self.GR.gO*self.gas.MO
        massR = self.GR.gR*self.Metal.MR
        return massR/massO

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