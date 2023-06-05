#####################################################
# Calculates the thermal conductivity (Slack)       #
# Version: 1.2                                      #
# Author: Kaja Bilińska ( k.bilinska@intibs.pl )    #
# Data: 05.06.2023                                  #
#####################################################

import numpy as np
from scipy import constants
import re
import pandas as pd
from bs4 import BeautifulSoup


class ThermalConductivity():

    def __init__(self):
        """
        Calculates thermal conductivity of the particular compound (half-Heusler)
        Analytical formula by Slack, calculations numerical
        """

        # get csv of stable hH (name of compound --> index)
        self.df_34_stable = pd.read_csv('stable_34_time.csv', index_col=0)

        # get values of atomic masses for the particular elements (csv provided with the use of ChatGPT3)
        self._elements_preprocessing()

    def _elements_preprocessing(self):
        """
        Get the csv with values of atomic mass (given in g and Da) for each of elements constituent below

        !!!
        For streamlit --> provide full period for Heusler + inforamtion about supported grouns&periods (csv on github)
        !!!
        """

        self.df_elements = pd.read_csv('stable_34_elements.csv', index_col=0)

        # list of elements constituent
        self.elements = ['Sc', 'Pd', 'Bi', 'Hf', 'Ni', 'Sn', 'Zr', 'Nb', 'Ru', 'Sb', 'Ti', 'V', 'Ir', 'Fe', 'Pt', 'Pb', 'Co', 'Ta', 'Rh', 'As', 'Te']

    def _calculate_thermal_conductivity(self, compound, XCF, T=300):
        """
        calculates thermal conductivity for compound given
        Args:
            compound (str): name of compound
            XCF (str): GGA / MBJ
            T (float, int > 0): value of temperature in K
        """

        # correct values of themperature
        if float(T) <= .0:
            print('Too low temperature given. Use Kelvins. Fin')
            exit()

        # get the and summary_mass of elements constituent in the particular alloy
        elements = re.findall('[A-Z][a-z]*', compound)
        l_summary_mass = []
        for element in elements:
            l_summary_mass.append(self.df_elements.loc[element, 'mass_g_1']) # mass_g_1 <-- mass of the single atom
        summary_mass = sum(l_summary_mass)

        number_of_atoms = len(element)

        # get value of volume of cell --> vasprun.xml
        with open(f'EIGENVAL/{compound}/{XCF}/vasprun.xml', 'r') as file:
            lines = file.readlines()

        for line in lines:
            if '<i name="volume">' in line:

                # parse the HTML-like file
                soup = BeautifulSoup(line, 'html.parser')

                # find the element with name="volume"
                element = soup.find('i', attrs={'name': 'volume'})

                try:
                    # extract the value after the <i> tag --> volume of cell
                    volume = float(element.text.strip())
                    self.df_34_stable.loc[compound, 'volume'] = volume
                except:
                    print(f'Something went wrong with {compound} ({XCF}). Fin')
                    exit()

                break

        # calculate cubic a**3, Angst**3 --> cm**-3
        # f_cell_volume_cm = d_volume[s_compound] * 1.8897259886 ** 3 * 10 ** (-24)  # cm**3
        f_cell_volume_cm = volume * 10 ** (-24)  # cm**3

        # density, g / cm**3
        rho = summary_mass / f_cell_volume_cm

        # density, g / cm**3 --> kg / m**3
        f_unit = 10**6
        rho_SI = rho / f_unit

        # elastic constants, calculated from VASP
        C11 = self.df_34_stable.loc[compound, 'C11']
        C12 = self.df_34_stable.loc[compound, 'C12']
        C44 = self.df_34_stable.loc[compound, 'C44']

        # some in-between calculations for bulk and Shear
        B_V = ( C11 + 2*C12 ) / 3
        B_R = ( ( C11 + C12 ) * C11 - 2 * C12**2 ) / ( 3 * ( C11 - C12 ) )

        G_V = ( ( C11 - C12 ) + 3*C44 ) / 5
        G_R_1 = ( 18 * B_V ) / ( (C11 + C12) * C11 - 2 * C12**2 )
        G_R_2 = 6 / ( C11 - C12 )
        G_R_3 = 9 / C44
        G_R = 15 *  ( G_R_1 + G_R_2 + G_R_3 )**(-1)

        # print(f'{compound}: C11={C11}, C12={C12}, C44={C44}')
        # print(f'{compound}: B_R={np.round(B_R,2)}, B_V={np.round(B_V,2)}')
        # exit()

        # bulk modulus
        B = ( B_V + B_R ) / 2
        self.df_34_stable.loc[compound, f'B_{XCF}'] = np.round(B,2)

        # Shear modulus
        G = ( G_V + G_R ) / 2
        self.df_34_stable.loc[compound, f'G_{XCF}'] = np.round(G,2)

        # Shear anisotropy factor
        A = ( 2 * C44 ) / ( C11 - C12 )
        self.df_34_stable.loc[compound, f'A_{XCF}'] = np.round(A,2)

        # Poisson's ratio
        v = ( 3 * B - 2 * G ) / ( 6 * B + 2 * G )
        self.df_34_stable.loc[compound, f'v_{XCF}'] = np.round(v,2)

        # Young's modulus
        Y = ( 9 * B * G ) / ( 3 * B + G )
        self.df_34_stable.loc[compound, f'Y_{XCF}'] = np.round(Y,2)

        # transverse sound velocity, m/s
        v_t = np.sqrt(G / rho_SI)
        self.df_34_stable.loc[compound, f'v_t_{XCF}'] = np.round(v_t,2)

        # longitudinal sound velocity, m/s
        v_l = np.sqrt((B + 4 * G / 3) / rho_SI)
        self.df_34_stable.loc[compound, f'v_l_{XCF}'] = np.round(v_l,2)

        # average sound velocity, m/s
        v_a = (1 / 3 * (2 / v_t ** 3 + 1 / v_l ** 3)) ** (-1 / 3)
        self.df_34_stable.loc[compound, f'v_a_{XCF}'] = np.round(v_a,2)

        # Gruneisen parameter
        gamma = (9 - 12 * (v_t / v_l) ** 2) / (2 + 4 * (v_t / v_l) ** 2)
        self.df_34_stable.loc[compound, f'gamma_{XCF}'] = np.round(gamma,2)

        # Boltzmann const, eV/K
        k_B = constants.physical_constants['Boltzmann constant in eV/K'][0]

        # volume of the unit cell, m**3 ???
        Omega = f_cell_volume_cm * 10 ** (-6)

        # Planck reduced (eV * s) / shall not be reduced?...
        hbar = 6.58211957 * 10 ** (-16)

        # theta_D - Debye temperature
        theta_D = (hbar * 2 * constants.pi / k_B) * (3 * number_of_atoms / (4 * constants.pi * Omega)) ** (1 / 3) * v_a
        self.df_34_stable.loc[compound, f'theta_D_{XCF}'] = np.round(theta_D,2)

        # theta_a - acoustic-mode Debye temperature
        theta_a = theta_D * number_of_atoms ** (-1 / 3)
        self.df_34_stable.loc[compound, f'theta_a_{XCF}'] = np.round(theta_a,2)

        # phys const
        A_const = 3.08 * 10 ** (-8) # -6

        # volume per atom in Angstroem**3
        delta = volume / number_of_atoms / 10 ** (-30)

        M_avg = summary_mass / number_of_atoms

        K_l = A_const * M_avg * theta_D * delta * gamma ** (-2) * T ** (-1)
        self.df_34_stable.loc[compound, f'K_l_{XCF}'] = np.round(K_l,2)

    def calculate_thermal_condictivity(self):

        XCFs = ['GGA', 'MBJ']

        for XCF in XCFs:

            for index, row in self.df_34_stable.iterrows():
                self._calculate_thermal_conductivity(index, XCF)

        print(self.df_34_stable)
        self.df_34_stable.to_csv('stable_34_Slack.csv')


TC = ThermalConductivity()
TC.calculate_thermal_condictivity()
