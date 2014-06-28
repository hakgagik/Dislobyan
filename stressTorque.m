function [S,T]=stressTorque(origMols,mols,latticeParams)

[oCenters, oDirections]=mol2VectorField(origMols,0);
[centers,directions]=mol2VectorField(mols,0);

%---------Everything in this section is explained in perfectHelix---------%
firstPt=2;
secondPt=3;

while 1
    planeNormal=cross(oCenters(firstPt,:)-oCenters(1,:),oCenters(secondPt,:)-oCenters(1,:));
    if (norm(planeNormal) >0.1)
        break
    else
        secondPt=secondPt+1;
    end
    if secondPt>length(oCenters)
        error('Couldn''t find a plane.')
    end
end

planeNormal=planeNormal/norm(planeNormal);
theta=acos(planeNormal(3));
phi=atan2(planeNormal(2),planeNormal(1));

Rz=[ cos(phi) sin(phi) 0;
    -sin(phi) cos(phi) 0;
        0        0     1;];
Ry=[cos(theta) 0 -sin(theta);
        0      1      0     ;
    sin(theta) 0  cos(theta);];

oCenters=(Ry*Rz*oCenters')';
oCenters=oCenters(:,1:2);

K=convhull(oCenters(:,1),oCenters(:,2));

aVect=oCenters(K(2),:)-oCenters(K(1),:);
bVect=oCenters(K(end-1),:)-oCenters(K(1),:);

aVect=aVect/norm(aVect); bVect=bVect/norm(bVect);

alpha=latticeParams(2,1); beta=latticeParams(2,2);

y0=(bVect(1)*cosd(beta)-aVect(1)*cosd(alpha))/(aVect(2)*bVect(1)-bVect(2)*aVect(1));
x0=(bVect(2)*cosd(beta)-aVect(2)*cosd(alpha))/(aVect(1)*bVect(2)-bVect(1)*aVect(2));

cVect=[x0 y0 sqrt(1-x0^2-y0^2)];

aVect=aVect*latticeParams(1,1);
bVect=bVect*latticeParams(1,2);
cVect=cVect*latticeParams(1,3);

%-------------------------------------------------------------------------%

% Recalculate the actual a, b, and c vectors of the lattice. These will be
% used later.
realaVect=(Ry*Rz)\[aVect 0]';
realbVect=(Ry*Rz)\[bVect 0]';
realcVect=(Ry*Rz)\cVect';

%Find each molecule's neighbors
tri=delaunay(oCenters); %Calculate the delaunay triangulation

%Find the area of each triangle and throw out the degenerate ones.
areas=cross([oCenters(tri(:,1),1:2)-oCenters(tri(:,2),1:2) zeros(length(tri),1)],[oCenters(tri(:,1),1:2)-oCenters(tri(:,3),1:2) zeros(length(tri),1)]);
areas=sqrt(sum(areas.^2,2));
tri=tri(areas>.1,:);



end