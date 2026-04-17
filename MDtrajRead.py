# -*- coding: utf-8 -*-
"""
PdbRead reads and saves .pdb trajectories into numpy arrays
Outputs:
1. atmGrp: Indicates if it's atom or hetatm
2. atmId: Saves atom ID
3. atmNm: Saves atom Name
4. resNm: Saves the name of the amino acids
5. resId: Saves the ID of the amino acids
6. segId: Saves chain ID
"""
from __future__ import division
from MDAnalysis import Universe as unvs
import numpy as np


def save_atm_res_seg(trajFileName):
    from MDAnalysis import Universe as unvs
    import numpy as np
    """
    Makes lists of atoms and residues name and ID 
    by using "Universe" code of "MDAnalysis"
    N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. 
    MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
    J. Comput. Chem. 32 (2011), 2319-2327, doi:10.1002/jcc.21787. PMCID:PMC3144279
    
    Input: a traj file name
    Outputs: lists of atoms and residues name and ID 
    """
    u = unvs(trajFileName)
    # 
    atmsAndRes = list(u.atoms)
    
	atmGrp, atmId, atmNm, resNm, resId, segId = [],[],[],[],[],[]
	
    for a in xrange(0,len(atmsAndRes)):
	
        # chang them to strings
        atmsAndResStr = str(atmsAndRes[a])
		
        # remove extra characters at the begining and at the end of str
        atmsAndResStr = atmsAndResStr[1:-1]
		
        # split them to words
        atmsAndResStrWrds = atmsAndResStr.split()
		
        # creating different lists
        atmGrp.append(atmsAndResStrWrds[0])
		
        atmId.append(int(atmsAndResStrWrds[1][0:-1]))
		
        atmNm.append(atmsAndResStrWrds[2])
		
        resNm.append(atmsAndResStrWrds[8][0:-1])
		
        resId.append(int(atmsAndResStrWrds[10]))
		
        segId.append(atmsAndResStrWrds[-1])
    
    del atmsAndRes, atmsAndResStr, atmsAndResStrWrds
	
    return np.array(atmGrp), np.array(atmId), np.array(atmNm), \
	       np.array(resNm), np.array(resId), np.array(segId) 


def save_traj_pos(trajFileName):

    """
    Saves the trajectoy into a 3D array "traj"

    Input: a traj file name
    Output: array with [frames, atoms, coordinates] 
    
    By using "Universe" code of "MDAnalysis"
    N. Michaud-Agrawal, E. J. Denning, T. B. Woolf, and O. Beckstein. 
    MDAnalysis: A Toolkit for the Analysis of Molecular Dynamics Simulations.
    J. Comput. Chem. 32 (2011), 2319-2327, doi:10.1002/jcc.21787. PMCID:PMC3144279
    """
    from MDAnalysis import Universe as unvs
    import numpy as np
    
    u = unvs(trajFileName)
	
    # save the trajectory into an numpy array
    traj = np.zeros([len(u.trajectory),u.coord.n_atoms,3])
    
    for frm in xrange(0,len(u.trajectory)):
	
        u.trajectory[frm]
		
        print(frm)
		
        traj[frm] = u.coord.positions
		
    return traj

	
atmGrp, atmId, atmNm, resNm, resId, segId = save_atm_res_seg(trajFileName)

traj = save_traj_pos(trajFileName);

return 	atmGrp, atmId, atmNm, resNm, resId, segId, traj
	
