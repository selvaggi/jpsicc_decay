##################################################################################################
#      HBoostDecay: Weighted events for H > c c~ mu- mu+                                         #
#      Yang Ma @ INFN Bologna        11/27/2022                                                  #
##################################################################################################
import os
import numpy as np
import pandas as pd


class Decay:
    def __init__(self, weight_dir, SDCfactor):
        # We use the following LDMEs
        # Color singlet (3s11), 1904.11542
        # Color octet (3s18), 1403.3612
        self.LDME_List = {"3s11": 1.0952, "3s18": 0.011}
        self.Charm_Mass = 1.5
        self.SDCfactor = SDCfactor

        # Load the momentums
        # Moment1,2,3,4 are for charm, anti-charm, muon, anti-muon, respectively
        self.Momentum_data_list = []
        for i in range(1, 5):
            self.Momentum_data_list.append(
                pd.read_csv(
                    "{}/Moment{}.dat".format(weight_dir, i),
                    sep="\s+",
                    dtype=np.float64,
                    names=["p0", "px", "py", "pz"],
                )
            )

        self.Momentum_data = [1, 2, 3, 4]
        for i in range(0, 4):
            self.Momentum_data[i] = self.Momentum_data_list[i].values

        # Load the Weights
        # "3s11","3s18" are the contributions for color-singlet and color-octet, respectively
        self.Weight_data_List = []
        self.Fork_States_jpsi = ["3s11", "3s18"]
        for state in self.Fork_States_jpsi:
            self.Weight_data_List.append(
                pd.read_csv(
                    "{}/weight{}.dat".format(weight_dir, state),
                    sep="\s+",
                    dtype=np.float64,
                    names=["Weight"],
                )
            )
        self.Weight_data = [1, 2]
        for i in range(0, 2):
            self.Weight_data[i] = self.Weight_data_List[i].values

    # Collect the prefactors that include the LDME
    def NRQCD_Parameter(self, WaveFunction_origin, LDME3s18, MC):
        prefactor_3s11 = WaveFunction_origin / (96 * MC * np.pi)
        prefactor_3s18 = LDME3s18 / (96 * MC)
        return (prefactor_3s11, prefactor_3s18)

    # Privde the event for a Mother with momentum Mother_Momentum
    def Decay_Events(self, Mother_Momentum):

        #   Generate a random number as the event index
        def Event_Index():
            return np.random.randint(0, len(self.Weight_data[1]) - 1)

        def Boost_Matrix(p):
            def Get_Velocity(p):
                s = slice(1, 4)
                return p[s] / p[0]

            def Get_LT_Gamma(v):
                return 1 / np.sqrt(1 - np.sum(np.square(v)))

            def Get_LT_Gamma2(v):
                if np.sum(np.square(v)) < 1e-15:
                    return 0.0
                else:
                    return (Get_LT_Gamma(v) - 1) / np.sum(np.square(v))

            v = Get_Velocity(p)
            v2 = np.sum(np.square(v))
            gam = Get_LT_Gamma(v)
            gam2 = Get_LT_Gamma2(v)

            matrix = np.array(
                [
                    [gam, v[0] * gam, v[1] * gam, v[2] * gam],
                    [v[0] * gam, 1 + (v[0] ** 2) * gam2, v[0] * v[1] * gam2, v[0] * v[2] * gam2],
                    [v[1] * gam, v[0] * v[1] * gam2, 1 + (v[1] ** 2) * gam2, v[1] * v[2] * gam2],
                    [v[2] * gam, v[0] * v[2] * gam2, v[1] * v[2] * gam2, 1 + (v[2] ** 2) * gam2],
                ]
            )
            return matrix

        def LT_Boost(p1, p2):
            return np.dot(Boost_Matrix(p1), p2)

        event = Event_Index()
        #    event=100

        #   Final state momentums in the mother rest frame
        pc = self.Momentum_data[0][event]
        pcb = self.Momentum_data[1][event]
        pmu = self.Momentum_data[2][event]
        pmub = self.Momentum_data[3][event]

        #   Boost to the mother moving frame
        pc2 = LT_Boost(Mother_Momentum, pc)
        pcb2 = LT_Boost(Mother_Momentum, pcb)
        pmu2 = LT_Boost(Mother_Momentum, pmu)
        pmub2 = LT_Boost(Mother_Momentum, pmub)

        prefactors = self.NRQCD_Parameter(
            self.LDME_List["3s11"], self.LDME_List["3s18"], self.Charm_Mass
        )
        weight_cs = self.SDCfactor * self.Weight_data[0][event][0] * prefactors[0]
        weight_co = self.SDCfactor * self.Weight_data[1][event][0] * prefactors[1]
        weight = weight_cs + weight_co

        return (event + 1, pc2, pcb2, pmu2, pmub2, weight, weight_cs, weight_co)
