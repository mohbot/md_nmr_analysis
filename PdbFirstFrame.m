function [atmGrpRes, atmNam, resNam, nAtm] = PdbFirstFrame(trajFileName) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reads the first frame of the trajectory and saves the topological information
% input argument trajectory of heavy atoms without chain ID
% outputs: 
% 1. 1st columnd: atom types(ATOM=1 or HETATM=2), 2nd column: residue number
%    for each atom in a array(atmGrpRes)
% 2. atom names in a cell(atmNam)
% 3. residue names in a cell(resNam)
% 4. number of atoms(nAtm)


fid = fopen(trajFileName,'r');

firstWrd = fscanf(fid, '%s',1);

i = 0;

while ~strcmp(firstWrd,'ENDMDL')
    switch firstWrd
        case 'ATOM'
          i = i+1;
          atmGrpRes(i,1) = 1;
          atmNam{i} = fscanf(fid, '%*9c%3c',1);
          resNam{i} = fscanf(fid, '%*1c%3c',1);
          atmGrpRes(i,2) = fscanf(fid, '%*3c%3d',1);
          fgetl(fid);
          firstWrd = fscanf(fid, '%s',1);
        case 'HETATM'
          i = i+1;
          atmGrpRes(i,1) = 2;
          atmNam{i} = fscanf(fid, '%*7c%3c',1);
          fgetl(fid);
          firstWrd = fscanf(fid, '%s',1);
        case 'TER'
          fgetl(fid);
          firstWrd = fscanf(fid, '%s',1);
        otherwise
          fgetl(fid);
          firstWrd = fscanf(fid, '%s',1);
    end
    
end

fclose(fid);

atmNam = atmNam';
resNam = resNam';
nAtm = size(atmGrpRes,1);

clear fid firstWrd i ans