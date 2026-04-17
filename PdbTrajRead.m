function [atmGrpRes, atmNam, resNam, nAtm, traj, time] = PdbTrajRead(trajFileName,nFrm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PdbTrajRead reads and saves the Gromacs pdb file(just heavy atoms) into matlab arrays
% Required functions: 
% PdbFirstFrame_reads the first frame of the trajectory
% PdbNumFrms_counts the number of frames of the trajectory if it's not known
% PdbTrajReadPos_reads and saves the cooridnates of the atoms into matlab arrays 
  
% outputs: 
% 1. 1st columnd: atom types(ATOM=1 or HETATM=2), 2nd column: residue number
%    for each atom in a array(atmGrpRes)
% 2. atom names in a cell(atmNam)
% 3. residue names in a cell(resNam)
% 4. number of atoms(nAtm) 
% 5. frames' timestamp(time)


    % 1. Get Topology
    [atmGrpRes, atmNam, resNam, nAtm] = PdbFirstFrame(trajFileName); 

    % 2. Determine Frame Count if not provided
    if nargin < 2 || isempty(nFrm)
        nFrm = PdbNumFrms(trajFileName);
    end

    % 3. Read Coordinates (Now outside the IF, so it always runs)
    [traj, time] = PdbTrajReadPos(trajFileName, nFrm, nAtm);
end
