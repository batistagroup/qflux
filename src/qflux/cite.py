import os


def cite(save_file=False, save_dir=None):
    """
    Function to generate bibtex citations for preprints of all QFlux parts. 

    
    """
    citations = """ 
                @article{qflux.2026_1, 
                         year    = {2026}, 
                         title   = {{QFlux}: Classical Foundations for Quantum Dynamics Simulation. Part I - Building Intuition and Computational Workflows}, 
                         author  = {Allen, Brandon C and Dan, Xiaohan and Cabral, Delmar G A and Vu, Nam P and Cianci, Cameron and Soudackov, Alexander V and Dutta, Rishab and Kais, Sabre and Geva, Eitan and Batista, Victor S}, 
                         journal = {{ChemRxiv}}, 
                         doi     = {10.26434/chemrxiv.10001765/v1}
                        }       

                @article{qflux.2026_2, 
                         year    = {2026}, 
                         title   = {{QFlux}: Quantum Circuit Implementations of Molecular Dynamics. Part {II} - Closed Quantum Systems}, 
                         author  = {Cabral, Delmar G A and Allen, Brandon C and Cianci, Cameron and Soudackov, Alexander V and Dan, Xiaohan and Vu, Nam P and Dutta, Rishab and Kais, Sabre and Geva, Eitan and Batista, Victor S}, 
                         journal = {{ChemRxiv}}, 
                         doi     = {10.26434/chemrxiv.10001766/v1}
                        }

                @article{qflux.2026_3, 
                         year    = {2026}, 
                         title   = {{QFlux}: Quantum Circuit Implementations of Molecular Dynamics. Part {III} - State Initialization and Unitary Decomposition}, 
                         author  = {Soudackov, Alexander V and Cabral, Delmar G A and Allen, Brandon C and Dan, Xiaohan and Vu, Nam P and Cianci, Cameron and Dutta, Rishab and Kais, Sabre and Geva, Eitan and Batista, Victor S}, 
                         journal = {{ChemRxiv}}, 
                         doi     = {10.26434/chemrxiv.10001767/v1}
                        }

                @article{qflux.2026_4, 
                         year    = {2026}, 
                         title   = {{QFlux}: An Open-Source Toolkit for Quantum Dynamics Simulations on Quantum Computers. Part {IV} - Dilation Method for Open Quantum Systems}, 
                         author  = {Dan, Xiaohan and Shivpuje, Saurabh and Wang, Yuchen and Cabral, Delmar G A and Allen, Brandon C and Khazaei, Pouya and Soudackov, Alexander V and Hu, Zixuan and Lyu, Ningyi and Geva, Eitan and Kais, Sabre and Batista, Victor S}, 
                         journal = {{ChemRxiv}}, 
                         doi     = {10.26434/chemrxiv.10001768/v1}
                        }

                @article{qflux.2026_5, 
                         year    = {2026}, 
                         title   = {{QFlux}: An Open-Source Toolkit for Quantum Dynamics Simulations on Quantum Computers. Part V - Adaptive Variational Quantum Algorithms for Open Quantum Systems}, 
                         author  = {Shivpuje, Saurabh and Soudackov, Alexander V and Dan, Xiaohan and Wang, Yuchen and Allen, Brandon C and Cabral, Delmar G A and Hu, Zixuan and Lyu, Ningyi and Geva, Eitan and Batista, Victor S and Kais, Sabre}, 
                         journal = {{ChemRxiv}}, 
                         doi     = {10.26434/chemrxiv.10001769/v2}
                        }

                @article{qflux.2026_6, 
                         year    = {2026}, 
                         title   = {{QFlux}: An Open-Source Toolkit for Quantum Dynamics Simulations on Quantum Computers. Part {VI} - The Generalized Quantum Master Equation}, 
                         author  = {Dan, Xiaohan and Khazaei, Pouya and Allen, Brandon C and Lyu, Ningyi and Wilson, Callie and Mulvihill, Ellen and Wang, Yuchen and Shivpuje, Saurabh and Kais, Sabre and Batista, Victor S and Geva, Eitan}, 
                         journal = {{ChemRxiv}}, 
                         doi     = {10.26434/chemrxiv.10001770/v1}
                        }"""
    print("You can cite the QFlux pre-prints using the following bibtex entries:")
    print(citations)
    
    if save_file:
        if not save_dir:
            save_dir = os.getcwd()
        fname = "QFlux_citation.bib"
        with open(os.path.join(save_dir, fname), 'w+') as file: 
            file.write("\n".join(citations))

    return 
