function nFrm = PdbNumFrms(trajFileName)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% counts the number of frames of a pdb file.
function nFrm = PdbNumFrms(trajFileName)

nFrm = 0;

fid = fopen(trajFileName,'r');

while ~feof(fid)
    firstWrd = fscanf(fid, '%s',1);
    if strcmp(firstWrd,'ENDMDL')
        nFrm = nFrm + 1;
        fgetl(fid);
    else
        fgetl(fid);
    end
end

clear fid firstWrd ans