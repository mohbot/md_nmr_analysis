function [traj, time] = PdbTrajReadPos(trajFileName, nFrm, nAtm)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% saves the coordinates and time stamps of the frames of a pdb file into matlab arrays
% traj dimensions: (frames, coordinates, atoms)


fid = fopen(trajFileName,'r');

traj = zeros(nFrm,3,nAtm);

time = zeros(nFrm,1);
  
    a = 0;
    f = 1;
    while ~feof(fid)
        firstWrd = fscanf(fid, '%s',1);
        switch firstWrd
            case 'TITLE'
                time(f) = fscanf(fid,'%*15c%d');
                
            case 'ATOM'
                a = a + 1;
                traj(f,:,a) = fscanf(fid, '%*26c%8f%8f%8f',3)';
                fgetl(fid);
                
            case 'HETATM'
                a = a + 1;
                traj(f,:,a) = fscanf(fid, '%*24c%8f%8f%8f',3)';
                fgetl(fid);
                
            case 'TER'
                fgetl(fid);
                
            case 'ANISOU'
                fgetl(fid);
            
                
            case 'ENDMDL'
                f = f + 1;
                a = 0;
                fgetl(fid);
                
            otherwise
        end
        
    end
   
fclose(fid);

clear fid a f firstWrd
