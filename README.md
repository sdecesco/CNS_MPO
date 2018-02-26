#CNS_MPO
            Python code that takes a SDF file as input to calculate the central nervous system multiparameter optimization (CNS MPO) score

This code calculate the score as detailed in this publication : 
Moving beyond Rules: The Development of a Central Nervous System Multiparameter Optimization (CNS MPO) Approach To Enable Alignment of Druglike Properties 
Travis T. Wager, Xinjun Hou, Patrick R. Verhoest, and Anabella Villalobos
ACS Chemical Neuroscience 2010 1 (6), 435-449
DOI: 10.1021/cn100008c

======================================================================

            CNS MPO Score and Solubility forecaster index (3.2)
            
=========================== Compatibility ============================

                        Compatible with Python 2.7 & >3.5
                        
===================== Prerequiste installations ======================

            Chemaxon Marvin suite WITH license
            cxcalc module required to work 	(used for pKa pred)
 					(make sure to add cxcalc to $PATH)
            by default : C:\Program Files (x86)\ChemAxon\MarvinBeans\bin
            How to use :
                        - input the sdf file name when requested (.sdf included)
                        - the script output an sdf with name_out.sdf
                                    This sdf contains the fields with :
                                            - CNS MPO Score
                                            - Solubility forecaster index (SFI)[optional]
                                            - bpKa,logD(7.4),logP,MW,HBD,Ar,TPSA
                                            - All individual components of MPO Score

======================================================================
