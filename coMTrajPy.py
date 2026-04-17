def traj_center_of_mass_res_bb_sc(traj, resId, atmNm):
    """
    Calculates the trajecotry of center of mass of residue, backbone and side chain
    
    Input: traj in [frames, atoms, coordinates] form 
    Output: array dimensions [residue, coordinates, frames]
    """
    import numpy as np
    # Atomic weights for Nitrogen, Carbon, Oxygen and Sulfur
    Nw, Cw, Ow, Sw = 14.01, 12.01, 16.00, 32.60
    weights = [Nw, Cw, Ow, Sw]
  
    
    nFrm = traj.shape[0]
	
    first_res = np.min(resId)
	
    last_res = np.max(resId)
    nRes = (last_res - first_res +1)
	
	# trajectory of residue center of mass
    CoM_res_traj = np.zeros([nRes,3,nFrm])
	
	# trajectory of backbone center of mass
    CoM_bb_traj = np.zeros([nRes,3,nFrm])
	
	# trajectory of backbone center of mass
    CoM_sc_traj = np.zeros([nRes,3,nFrm])
	
    bb_weight = np.sum(weights[0:3])+weights[1]
    
    # Change the order of dimensions into [atoms, coordinates, frames]
    traj = traj.swapaxes(1,2).swapaxes(2,0)
    
    for r in xrange(first_res-1, last_res):
        # Extracting trajectory of each residue
        res_traj = traj[resId==r+1,:,:]
        
        res_atm_name = atmNm[resId==r+1]
        
        # Weighting the traj of residue
        res_weight = 0
        res_weighted_traj = np.zeros([len(res_atm_name),3,nFrm])
		
        for a in xrange(len(res_atm_name)):
            if res_atm_name[a][0]=='N' :
                res_weighted_traj[a,:,:] = res_traj[a,:,:]*weights[0]
                res_weight += weights[0]
            elif res_atm_name[a][0]=='C' :
                res_weighted_traj[a,:,:] = res_traj[a,:,:]*weights[1]
                res_weight += weights[1]
            elif res_atm_name[a][0]=='O' :
                res_weighted_traj[a,:,:] = res_traj[a,:,:]*weights[2]
                res_weight += weights[2]
            elif res_atm_name[a][0]=='S' :
                res_weighted_traj[a,:,:] = res_traj[a,:,:]*weights[3]
                res_weight += weights[3]
            
        # Center of mass of residue
        CoM_res_traj[r,0,:] = res_weighted_traj[:,0,:].sum(0)/res_weight #x
        CoM_res_traj[r,1,:] = res_weighted_traj[:,1,:].sum(0)/res_weight #y
        CoM_res_traj[r,2,:] = res_weighted_traj[:,2,:].sum(0)/res_weight #z
        
        # Centre of mass of backbone
        # Location of backbone atoms
        bb_atms = ismember(res_atm_name, np.array(['N','CA','C','O']))
        
        CoM_bb_traj[r,0,:] = res_weighted_traj[bb_atms,0,:].sum(0)/bb_weight #x
        CoM_bb_traj[r,1,:] = res_weighted_traj[bb_atms,1,:].sum(0)/bb_weight #y
        CoM_bb_traj[r,2,:] = res_weighted_traj[bb_atms,2,:].sum(0)/bb_weight #z
        
        # Centre of mass of side chain
        if  len(res_atm_name) == 4:
            CoM_sc_traj[r,0,:] = np.nan
            CoM_sc_traj[r,1,:] = np.nan
            CoM_sc_traj[r,2,:] = np.nan
        else:
            CoM_sc_traj[r,0,:] = res_weighted_traj[~bb_atms,0,:].sum(0)/(res_weight - bb_weight) #x
            CoM_sc_traj[r,1,:] = res_weighted_traj[~bb_atms,1,:].sum(0)/(res_weight - bb_weight) #y
            CoM_sc_traj[r,2,:] = res_weighted_traj[~bb_atms,2,:].sum(0)/(res_weight - bb_weight) #z
        
        del r, res_traj, res_atm_name, res_weight, a, res_weighted_traj, bb_atms
        
    return CoM_res_traj, CoM_bb_traj, CoM_sc_traj
	
	
