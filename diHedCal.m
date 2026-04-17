function [phi, psi, omega] = diHedCal(bbTrajPos,atmGrpRes)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% diHedCal calculates phi, psi and omega and for each reasidue in each frame
% inputs:
% trajecotry of backbone(bbTrajPos) 
% 1st columnd: atom types(ATOM=1 or HETATM=2), 2nd column: residue number(atmGrpRes) 

nFrm = size(bbTrajPos,3);
nRes = size(bbTrajPos,1);

[phi, psi, omega] = deal(zeros(nRes, nFrm));

for f=1:nFrm
    for r=min(atmGrpRes(:,2))+1:max(atmGrpRes(:,2))-1
            
            % for Phi calculation for the second residue
            % phi 2  C1 N2 CA2 C2
            coord(1,:) = bbTrajPos(r-1,7:9,f);
            coord(2,:) = bbTrajPos(r,1:3,f);
            coord(3,:) = bbTrajPos(r,4:6,f);
            coord(4,:) = bbTrajPos(r,7:9,f);
            
            phi(r,f) = localCalculateTorsionAngle(coord)*180/pi;
            clear coord
            
            % for Psi calculation for the first residue
            % psi 1  N1 CA1 C1 N2
            coord(1,:) = bbTrajPos(r,1:3,f);
            coord(2,:) = bbTrajPos(r,4:6,f);
            coord(3,:) = bbTrajPos(r,7:9,f);
            coord(4,:) = bbTrajPos(r+1,1:3,f);
            
            psi(r,f) = localCalculateTorsionAngle(coord)*180/pi;
            clear coord
            
            % for Omega calculation for the second residue
            % omega 2  CA1 C1 N2 CA2
            coord(1,:) = bbTrajPos(r,4:6,f);
            coord(2,:) = bbTrajPos(r,7:9,f);
            coord(3,:) = bbTrajPos(r+1,1:3,f);
            coord(4,:) = bbTrajPos(r+1,4:6,f);
            omega(r,f) = localCalculateTorsionAngle(coord)*180/pi;
            
    end
    f
       
end
        
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [torsionAngle] = localCalculateTorsionAngle(coord)
% Evaluate the torsion angles.

p1 = coord(1,:);
p2 = coord(2,:);
p3 = coord(3,:);
p4 = coord(4,:);

a = p2 - p1;
b = p3 - p2;
c = p4 - p3;

a = a/norm(a,2);
b = b/norm(b,2);
c = c/norm(c,2);

b_cross_c = [b(2).*c(3) - b(3).*c(2);
    b(3).*c(1) - b(1).*c(3);
    b(1).*c(2) - b(2).*c(1)];

x = -sum(conj(a).*c) + ((sum(conj(a).*b))*(sum(conj(b).*c)));
y = sum(conj(a).*b_cross_c');
torsionAngle = localNewAngle(x,y);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%           
function ang = localNewAngle(x,y)
% Calculate the angle represented by (x,y). The angle lies between -pi and
% pi.

ang = 0; % This is the default value. In case y == 0 and x ~= 0, ang = 0.
if (x ~= 0) && (y~= 0)
    c = x./sqrt(x.^2 + y.^2);
    ang = sign(y)*acos(c);
elseif (x == 0)
    if (y > 0)
        ang = pi/2;
    elseif (y < 0)
        ang = -pi/2;
    end
end
