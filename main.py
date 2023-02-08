import numpy as np
import pandas as pa
from SingleParticle import SingleParticle
from HotSingleParticle import HotSingleParticle


def particle_combustion():
    # Extracting data from excel Sheet
    metal_data = pa.read_excel('/Users/dihiaidrici/Desktop/ParticleCombustion.xlsx', sheet_name='MetalData', header=1, index_col=0)
    gas_data = pa.read_excel('/Users/dihiaidrici/Desktop/ParticleCombustion.xlsx', sheet_name='GasData', header=1, index_col=0)
    global_reaction = pa.read_excel('/Users/dihiaidrici/Desktop/ParticleCombustion.xlsx', sheet_name='GlobalReaction', header=0, index_col=0)

    # radius = pa.read_excel('/Users/dihiaidrici/Desktop/ParticleCombustion.xlsx', sheet_name='radiusRange', header=0, index_col=0)
    # remember that indexing begins at 0, Rows: Al[0], Zr[1]
    Al = metal_data.iloc[0]
    Zr = metal_data.iloc[1]
    # Gas data
    ArO2 = gas_data.iloc[0]
    AirO2 = gas_data.iloc[1]
    AirCO2 = gas_data.iloc[2]
    #Pick a global reaction
    Al_21_O2 = global_reaction.iloc[0]
    Al_30_O2 = global_reaction.iloc[1]
    Al_40_O2 = global_reaction.iloc[2]

    Al_combustion = SingleParticle(Al, AirO2, Al_21_O2)
    Al_combustion.reaction_rate()

    # Al_combustion = HotSingleParticle(Al, AirO2)
    # Al_combustion.reaction_rate()

    # Zr_combustion = SingleParticle(Zr, ArO2, radius)
    # Zr_combustion.reaction_rate()


if __name__ == '__main__':
    particle_combustion()
