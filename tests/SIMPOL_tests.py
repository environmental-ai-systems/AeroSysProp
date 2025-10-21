'''programme for testing the SIMPOL implementation in AeroSysProp'''

# list of changes to SIMPOL.smarts
# group 2 (carbon number on acid side of amide) left as is
# groups 23, 24, 25 (primary, secondary and tertiary amides had
# the requirement for the carbon in an amide group to be bonded
# to another carbon removed, as it is possible for that
# carbon to be bonded to other atoms and still be an amide,
# e.g. formamide)
# groups 18, and 19 for primary and secondary amines were
# extended so that as well as flagging for amines as
# defined conventionally, they also allowed inclusion
# of hydroxylamines (e.g.dimethyl-hydroxylamine, 
# sometimes called (N,N-Dimethylhydroxylamine)), 
# which Pankow and Asher (2008) classify
# as amines in their Table 1c, and therefore Simon O'Meara
# interprets as meaning for SIMPOL.1, hydroxylamines should
# activate the amine SMARTS. Other changes, see notes in groups.py

# Issues with Bezaatpour et al. (2025), to gain some of the manually
# estimated vapour pressures in their Table 2, specifically dimethyl-hydroxylamine,
# they must activate the alkyl hydroxyl group, however footnote h of Table 2 of
# Pankow and Asher (2008) says this group should only be activated when attached
# to a non-aromatic carbon, whereas in dimethyl-hydroxylamine it is attached to 
# a Nitrogen. Also, the manually estimated and VaPOrS estimated vapour pressure
# in their Table 2 for diethyl-ethanedicarboxylate does not agree with the
# estimate for that SMILES string when using T=333.15 in the paper-cited
# VaPOrS.py python programme. However, when I search this SMILES string in
# ChemSpider it calls the molecule diethyl 2-butenedioate. diethyl-ethanedicarboxylate
# may be representing diethyl ethanedioate (diethyl oxalate more commonly) which is 
# listed as an ester in Table 1d of
# Pankow and Asher (2008) and the internet says is more commonly called 
# diethyl oxalate, and has the SMILES CCOC(=O)C(=O)OCC. This SMILES also does not
# return the provide manually estimated vapour pressure, but at least it gives
# consistency with a molecule named in the original paper. The SMILES that
# does reproduce the manually estimated vapour pressure is CCOC(=O)CCC(=O)OCC, but
# as described above this molecule does not seem to correspond to any in the
# Pankow and Asher (2008) paper.
# Note that the Bezaatpour et al. 2025 paper and corresponding VaPOrS.py
# programme wrongly asigns a ketone group to structures like C(=O)O[O], such
# as in [O]OC(=O)CC(=O)C=O, which is incorrect as the carbonyl C in the ketone
# must bond to two other carbons, not an oxygen because of the definition of
# a ketone, this is true for all compounds in their Table 4 and some compounds in
# their Table 5, whilst other compounds in their table 5 relate to issues with the
# old AeroSysProp version, including: OC1COC(=O)C1=O 

import sys
import os
import openbabel.pybel as pybel
import numpy as np
import math
import scipy.constants as si
	
def SIMPOL_test(): # function to test SIMPOL implementation in AeroSysProp
    # include aerosysprop directory in search path
    parent_directory = os.path.abspath('../')
    sys.path.append(parent_directory)

    from aerosysprop.methods import vapour_pressures
    from aerosysprop import data

    # get the table of coefficient values for each contributing group
    # from a AeroSysProp, which is based on
    # Table 5 of https://doi.org/10.5194/acp-8-2773-2008
    b_onevals = np.array(list(data.SIMPOL_1.values()))
    b_twovals = np.array(list(data.SIMPOL_2.values()))
    b_thrvals = np.array(list(data.SIMPOL_3.values()))
    b_fouvals = np.array(list(data.SIMPOL_4.values()))
    
    # set temperatures to estimate vapour pressures at (K), the
    # 333.15 K referenced in Section 3.1 of 
    # Bezaatpour et al. 2025
    T = 333.15

    # calculate pure component liquid vapour pressures log10(atm) based
    # on the SIMPOL.1 equation in Eq. 6 of 
    # https://doi.org/10.5194/acp-8-2773-2008
    # and using the coefficients in Table
    bk = (b_onevals/T+b_twovals+b_thrvals*T+b_fouvals*math.log(T)).reshape(-1, 1)

    # manually set number of occurrences of each contributing group
    # per component, and get SMILES strings of each component
    [vki, SMILES] = manual_v_set()
    # get manually estimated vapour pressures (log10(atm))
    pLzero_ref_man = np.sum(vki*bk.reshape(1, -1), axis= 1)

    # prepare to hold AeroSysProp estimates (log10(atm))
    pLzero_aero = []

    for SMILESi in range(len(SMILES)):
        # convert SMILES string to pybel object
        Pybel_object = pybel.readstring('smi', SMILES[SMILESi])
        
        # estimate pure component liquid vapour pressure log10(atm)
        # from AeroSysProp implementation of SIMPOL.1
        pLzero_aero.append(vapour_pressures.simpol(Pybel_object, T))

    return()

# Function containing Simon O'Meara's manual counting of functional
# groups per component
def manual_v_set():
    
    # prepare results array with components in rows and functional 
    # groups in columns
    vki = np.zeros((63, 31))
    # the zeroeth group always has a count of 1 (see text around Eq. 6
    # of Pankow and Asher 2008)
    vki[:, 0] = 1
    # for each component in Tables 2, 3, 4 and 5 of
    # Bezaatpour et al. 2025, manually supply count on contributing groups
    vki[0, [1, 2, 22]] = [1, 1, 1] # start of Table 2 of Bezaatpour et al. 2025
    vki[1, [1, 2, 23]] = [3, 1, 1]
    vki[2, [1, 2, 23]] = [2, 1, 1]
    vki[3, [1, 2, 24]] = [5, 1, 1]
    vki[4, [1, 2, 24]] = [3, 1, 1]
    vki[5, [1, 7, 19]] = [2, 0, 1]
    vki[6, [1, 3, 11]] = [11, 1, 1]
    vki[7, [1, 3, 11]] = [11, 1, 1]
    vki[8, [1, 3, 11]] = [10, 1, 1]
    vki[9, [1, 3, 11]] = [9, 1, 1]
    vki[10, [1, 3, 11]] = [8, 1, 1]
    vki[11, [1, 3, 11]] = [8, 1, 1]
    vki[12, [1, 3, 11]] = [10, 1, 2]
    vki[13, [1, 3, 11]] = [10, 1, 2]
    vki[14, [1, 3, 11]] = [10, 1, 2]
    vki[15, [1, 11]] = [10, 2]
    vki[16, [1, 11]] = [6, 2]
    vki[17, [1, 16, 30]] = [5, 1, 1]
    vki[18, [1, 16, 30]] = [4, 1, 1] # end of Table 2 of Bezaatpour et al. 2025
    vki[19, [1, 8]] = [1, 1] # start of Table 3 of Bezaatpour et al. 2025
    vki[20, [1, 10]] = [1, 1]
    vki[21, [1, 11]] = [2, 1]
    vki[22, [1, 4, 5]] = [15, 2, 2] # end of Table 3 of Bezaatpour et al. 2025
    vki[23, [1, 8, 9]] = [4, 1, 1] # start of Table 4 of Bezaatpour et al. 2025
    vki[24, [1, 4, 10]] = [8, 1, 1]
    vki[25, [1, 4, 10]] = [9, 1, 1]
    vki[26, [1, 4, 9]] = [9, 1, 1]
    vki[27, [1, 4, 8]] = [9, 1, 1]
    vki[28, [1, 4, 8]] = [9, 1, 1]
    vki[29, [1, 4, 7, 9]] = [10, 1, 1, 1]
    vki[30, [1, 4, 9]] = [10, 1, 1]
    vki[31, [1]] = [2]
    vki[32, [1, 8, 9]] = [5, 1, 1]
    vki[33, [1, 9]] = [7, 3]
    vki[34, [1, 9]] = [5, 2]
    vki[35, [1, 7, 9]] = [4, 1, 1]
    vki[36, [1, 7, 9]] = [5, 1, 2]
    vki[37, [1, 7, 9]] = [7, 1, 2]
    vki[38, [1, 7, 9]] = [5, 1, 1]
    vki[39, [1, 8]] = [3, 1]
    vki[40, [1, 8]] = [2, 1]
    vki[41, [1, 7]] = [2, 1]
    vki[42, [1, 8, 9, 15]] = [7, 1, 2, 1] # Table 4 ends
    vki[43, [1, 4, 5, 8, 13]] = [6, 1, 1, 1, 1]# Table 5 of Bezaatpour et al. 2025 begins
    vki[44, [1, 4, 7, 9, 11]] = [4, 1, 1, 1, 1]
    vki[45, [1, 4, 5, 6, 11]] = [4, 1, 1, 1, 1]
    vki[46, [1, 5, 8, 9]] = [5, 1, 1, 1]
    vki[47, [1, 5, 7, 8, 9]] = [6, 1, 1, 1, 1]
    vki[48, [1, 8, 11]] = [4, 2, 1]
    vki[49, [1, 8, 10, 11]] = [4, 1, 1, 1]
    vki[50, [1, 4, 8, 13]] = [4, 1, 1, 1]
    vki[51, [1, 8]] = [2, 1]
    vki[52, [1, 7, 8]] = [3, 1, 1]
    vki[53, [1, 7]] = [2, 1]
    vki[54, [1, 5, 8]] = [4, 1, 1]
    vki[55, [1, 5, 8]] = [4, 1, 1]
    vki[56, [1, 5, 8, 16]] = [4, 1, 1, 1]
    vki[57, [1, 4, 5, 7, 16, 26]] = [6, 2, 1, 3, 1, 1]
    vki[58, [1, 4, 5, 7, 16, 26]] = [6, 2, 1, 3, 1, 1]
    vki[59, [1, 4, 5, 7, 16, 26, 27]] = [6, 2, 1, 3, 1, 1, 1]
    vki[60, [1, 4, 5, 7, 15, 16, 26]] = [6, 2, 1, 2, 1, 1, 1]
    vki[61, [1, 4, 5, 7, 15, 16, 26]] = [6, 2, 1, 2, 1, 1, 1]
    vki[62, [1, 4, 5, 7, 15, 16, 26, 27]] = [6, 2, 1, 2, 1, 1, 1, 1]


    # SMILES strings for considered components, taken from 
    # Bezaatpour et al. 2025; 130 components in total
    SMILES = ['C(=O)N', # formamide (0) # start of Table 2 of Bezaatpour et al. 2025
              'CCNC=O', # ethyl-formamide
              'CNC=O ', # methyl-formamide
              'CCN(CC)C=O', # diethyl-formamide
              'CN(C)C=O', # dimethyl-formamide
              'CN(C)O', # dimethyl-hydroxylamine (N,N-Dimethylhydroxylamine)
              'CCCCOC(=O)c1ccccc1', # n-butyl-benzoate (n-butyl benzoate in Pankow and Asher 2008)
              'CC(C)COC(=O)c1ccccc1', # 2-methyl-propyl-benzoate
              'CCCOC(=O)c1ccccc1', # n-propyl-benzoate
              'CCOC(=O)c1ccccc1', # ethyl-benzoate
              'COC(=O)c1ccccc1', # methyl-benzoate (10)
              'CC(=O)Oc1ccccc1', # acetic-acid,-phenyl-ester (phenyl acetate in Pankow and Asher 2008)
              'COC(=O)c1ccccc1C(=O)OC', # dimethyl-1,2-benzenedicarboxylate
              'COC(=O)c1cccc(c1)C(=O)OC', # dimethyl-benzene-1,3-dicarboxylate
              'COC(=O)c1ccc(cc1)C(=O)OC', # dimethyl-benzene-1,4-dicarboxylate
              'CCCCOC(=O)C(=O)OCCCC ', # di-n-butyl-ethanedicarboxylate
              'CCOC(=O)C(=O)OCC', # diethyl-ethanedicarboxylate (diethyl ethanedioate in Pankow and Asher 2008 and diethyl oxalate more commonly, note that SMILES string different to Table 2 of Bezaatpour et al. 2025 as they mistakenly included a carbon-carbon double bond, whereas neither the double bond nor the carbons in it are present in diethyl oxalate)
              'CCOC(=O)C(C)[N+](=O)[O-]', # ethyl-2-nitropropionate 
              'CC(C(=O)OC)[N+](=O)[O-]', # methyl-2-nitro-propionate 
              'C=O', # formaldehyde 
              'OC=O', # formic acid
              'COC=O', # methyl ester
              'C/C1=C/CCC(=C)C2CC(C)(C)C2CC\\1', # beta-caryophyllene
              '[O]OC(=O)CC(=O)C=O', # start of Table 4 of Bezaatpour et al. 2025
              '[O]OC(=O)C1CC(C(=O)O)C1(C)C',
              '[O]OC(=O)CC1CC(C(=O)O)C1(C)C',
              '[O]OC(=O)C1CC(C(=O)C)C1(C)C',
              'O=CCC1CC(C(=O)[O])C1(C)C',
              'O=CCC1CC(C(=O)O[O])C1(C)C',
              '[O]OC(=O)CC1CC(C(=O)CO)C1(C)C',
              '[O]OC(=O)CC1CC(C(=O)C)C1(C)C',
              'CC(=O)O[O]',
              '[O]OC(=O)CC(=O)CC=O',
              '[O]OC(=O)CC(=O)CC(=O)C(=O)C',
              '[O]OC(=O)CC(=O)C(=O)C',
              'CC(=O)C(O)C(=O)O[O]',
              '[O]OC(=O)CC(=O)C(=O)CO',
              '[O]OC(=O)CC(=O)CC(O)C(=O)C',
              'CC(=O)C(O)CC(=O)O[O]',
              '[O]OC(=O)CC=O',
              '[O]OC(=O)C=O',
              'OCC(=O)O[O]',
              'O=CC(=O)CC(=O)C(C)(ON(=O)=O)C(=O)O[O]',
              '[O]OC(=O)C1OC1C=CC=O', # start of Table 5 of Bezaatpour et al. 2025
              'OC1COC(=O)C1=O',
              'C1OC(=O)C=C1',
              '[O]OC(=O)C=CC(=O)C=O',
              'OC(C=O)C(=O)C=CC(=O)O[O]',
              'O=CCOC(=O)C=O',
              'O=CCOC(=O)C(=O)O',
              '[O]OC(=O)C1OC1C=O',
              '[O]OC(=O)C=O',
              'OC(C=O)C(=O)O[O]',
              'OCC(=O)O[O]',
              'O=CC=CC(=O)[O]',
              '[O]OC(=O)C=CC=O',
              'O=CC(=CC(=O)[O])N(=O)=O',
              'O=N(=O)C12OOC(C2O)C(O)([O])C(=C1)O',
              '[O]OC1(O)C(=CC2(OOC1C2O)N(=O)=O)O',
              'OOC1(O)C(=CC2(OOC1C2O)N(=O)=O)O',
              'O=N(=O)OC1C2OOC1(C=C(O)C2([O])O)N(=O)=O',
              '[O]OC1(O)C(=CC2(OOC1C2ON(=O)=O)N(=O)=O)O',
              'OOC1(O)C(=CC2(OOC1C2ON(=O)=O)N(=O)=O)O', # end of Table 5 of Bezaatpour et al. 2025
              'CC(=O)C(=O)C(O)=CC(=O)C(O)C=O', # start of Table 6 of Bezaatpour et al. 2025
            'CC(=O)C(=O)C(O)=CC(=O)C(OO)C=O', 
            'CC(=O)C(=O)C(O)=CC(OO)C(=O)C=O', 
            'CC(=O)C(=O)C(O)=CC(O[O])C(O)C=O', 
            'CC(=O)C(=O)C(O)=CC(ON(=O)=O)C(O)C=O', 
            'CC(=O)C(=O)C(O)=CC(OO)C(O)C=O', 
            'CC(=O)C(=O)C(O)=CC(O)C(OO)C=O', 
            'CC(=O)C(=O)C(O)=CC(=O)C(O)C(=O)OO', 
            'CC(=O)C(=O)C(O)=CC(O[O])C(OO)C=O', 
            'CC(=O)C(=O)C(O)C(OO)C=CC(O[O])=O', 
            'CC(=O)C(=O)C(O)=CC(ON(=O)=O)C(OO)C=O', 
            'CC(=O)C(=O)C(O)=CC(O)C(O)C(=O)OO', 
            'CC(=O)C(=O)C(O)=CC(=O)C(OO)C(=O)OO', 
            'CC(=O)C(=O)C(O)=CC(O[O])C(O)C(=O)OO', 
            'CC(=O)C(=O)C(O)=CC(ON(=O)=O)C(O)C(=O)OO', 
            'CC(=O)C(=O)C(O)=CC(OO)C(OO)C(=O)O[O]', 
            'CC(=O)C(=O)C(O)=CC(O[O])C(OO)C(=O)OO', 
            'CC(=O)C(=O)C(O)=CC(OO)C(OO)C(=O)ON(=O)=O', 
            'CC(=O)C(=O)C(O)=CC(ON(=O)=O)C(OO)C(=O)OO', 
            'CC(=O)C(=O)C(O)=CC(OO)C(OO)C(=O)OO', 
            '[O]OCC(=O)C(=O)C(O)=CC(OO)C(O)C(OO)=O', 
            'O=N(=O)OCC(=O)C(=O)C(O)=CC(OO)C(O)C(OO)=O', 
            'O=C(C=O)C(O)=CC(=O)C(O)C=O', 
            'O=C(C=O)C(O)=CC(=O)C(OO)C=O', 
            'O=C(C=O)C(O)=CC(OO)C(=O)C=O', 
            'O=C(C=O)C(O)=CC(O[O])C(O)C=O', 
            'O=C(C=O)C(O)=CC(O)C(OO)C=O', 
            'O=C(C=O)C(O)=CC(=O)C(O)OO', 
            'O=C(C=O)C(O)=CC(=O)C(O)C(=O)OO', 
            'O=C(C=O)C(O)=CC(O[O])C(OO)C=O', 
            'O=C(C=O)C(O)=CC(OO)C(OO)C=O', 
            'O=C(C=O)C(O)=CC(=O)C(OO)C(=O)OO', 
            'O=C(C=O)C(O)=CC(O[O])C(O)C(=O)OO', 
            'O=C(C=O)C(O)=CC(O[O])C(OO)C(=O)OO', 
            'O=C(C=O)C(O)=CC(OO)C(OO)C(=O)OO', 
            'O=CCC(=O)C(O)=CC(=O)C(O)C=O', 
            'O=CCC(=O)C(O)=CC(=O)C(OO)C=O', 
            'O=CCC(=O)C(O)=CC(OO)C(=O)C=O', 
            'O=CCC(=O)C(O)=CC(O[O])C(O)C=O', 
            'O=CCC(=O)C(O)=CC(ON(=O)=O)C(O)C=O', 
            'O=CCC(=O)C(O)=CC(O)C(OO)C=O', 
            'O=CCC(=O)C(O)=CC(=O)C(O)C(=O)OO', 
            'O=CCC(=O)C(O)=CC(O[O])C(OO)C=O',
            'O=CCC(=O)C(O)=CC(ON(=O)=O)C(OO)C=O',
            'O=CCC(=O)C(O)=CC(OO)C(OO)C=O',
            'O=CCC(=O)C(O)=CC(O[O])C(O)C(=O)OO',
            'O=CCC(=O)C(O)=CC(ON(=O)=O)C(O)C(=O)OO',
            'O=CC(=O)C(=O)C(O)=CC(OO)C(O)C(=O)OO',
            'O=CCC(=O)C(O)=CC(OO)C(OO)C(=O)O[O]',
            '[O]OC(C=O)C(=O)C(O)=CC(OO)C(OO)C=O',
            'O=CCC(=O)C(O)=CC(O[O])C(OO)C(OO)=O',
            'O=CCC(=O)C(O)=CC(OO)C(OO)C(=O)ON(=O)=O',
            'O=N(=O)OC(C=O)C(=O)C(O)=CC(OO)C(OO)C=O',
            'O=CCC(=O)C(O)=CC(ON(=O)=O)C(OO)C(OO)=O',
            'O=CCC(=O)C(O)=CC(OO)C(OO)C(=O)OO',
            'O=CC(OO)C(=O)C(O)=CC(OO)C(OO)C=O',
            'O=CC(O[O])C(=O)C(O)=CC(OO)C(O)C(=O)OO',
            'O=CCC(=O)C(=O)C(OO)C(OO)C(O)C(=O)O[O]',
            'O=CC(ON(=O)=O)C(=O)C(O)=CC(OO)C(O)C(=O)OO',
            '[O]OC(=O)C(OO)C(=O)C(O)=CC(OO)C(OO)C=O',
            'O=CCC(=O)C(=O)C(OO)C(OO)C(OO)C(=O)O[O]',
            'O=N(=O)OC(=O)C(OO)C(=O)C(O)=CC(OO)C(OO)C=O',
            'OOC(=O)C(OO)C(=O)C(O)=CC(OO)C(OO)C=O',
            'O=CCC(=O)C(O)=CC(C(O)C=O)OOOC(C=O)C=CC(O)C(=O)CC=O',
            'O=CCC(=O)C(O)=CC(C(O)C=O)OOOC1C=CC2OOC1C2(O)CC=O',
            'O=CCC(=O)C(O)=CC(C(O)C=O)OOC(C(O)C=O)C=C(O)C(=O)CC=O',
            'O=CCC(=O)C(O)=CC(C(OO)C=O)OOC(C(OO)C=O)C=C(O)C(=O)CC=O',
    ]

    return(vki, SMILES)

# Table 2 of Pankow and Asher 2008 supplement 
# https://doi.org/10.5194/acp-8-2773-2008
# gives the
# Antoine coefficients for estimating pure component
# liquid saturation vapour pressures (atm), however these
# values are for the reference data used to train and 
# test SIMPOL.1, so do not have much use when testing
# whether AeroSysProp estimates match SIMPOL.1 estimates
def Antoine_2008_res(T):

    # prepare to hold results
    pLzero_ref_ac = []

    # formamide (uses eq. 0.1 and 0.2 of supplement)
    pLzero_ref_ac.append(0.9869*10**(7.5852-3881.305/(T+27.655)))
    # ethyl-formamide (uses eq. 0.7 of supplement)
    pLzero_ref_ac.append(np.exp((471.2/(si.R*58.44))*(1.8*(1-(58.44/T))+0.8*math.log(58.44/T))))
    
    # convert to log10(atm)
    pLzero_ref_ac = np.log10(pLzero_ref_ac)

    return(pLzero_ref_ac)

SIMPOL_test() # call function to test SIMPOL implementation in AeroSysProp