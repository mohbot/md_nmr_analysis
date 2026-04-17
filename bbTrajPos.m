function bbTrajPos = bbTrajPos(traj,atmGrpRes,atmNam)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% save the position of backbone atoms
% Finding the four points coordinates to calculate Psi
% for the first residue and Phi for the seconed one
% Storing the residue number and BB coordinates in bbResPos


nFrm = size(traj,1);
nRes = max(atmGrpRes(:,2))- min(atmGrpRes(:,2))+1;

bbTrajPos = zeros(nRes,12,nFrm);

traj = permute(traj,[3 2 1]);

% if first residue is not res number 1
r = min(atmGrpRes(:,2));
% atms
a=1;

while a < size(atmGrpRes,1)-1
    
    if (isequal(atmNam{a}, 'N  ')) && (atmGrpRes(a,2) == r)
        % N coordinates
        bbTrajPos(r,1:3,:) = traj(a,:,:);
        
        % CA coordinates
        bbTrajPos(r,4:6,:) = traj(a+1,:,:);
        
        
        while ~(isequal(atmNam{a}, 'C  ')) && (atmGrpRes(a,2) == r)
            a=a+1;
        end
        
        % C coordinates
        bbTrajPos(r,7:9,:) = traj(a,:,:);
        
        % O coordinates
        bbTrajPos(r,10:12,:) = traj(a+1,:,:);
        
        a=a+2;
        % prints the residue number
        r = r+1
    end
   
end
clear a r