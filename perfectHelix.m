function mols = perfectHelix(infile,order,latticeParams,burgers)
% Raises the selected molecules by a fraction of the burgers vector based
% on the number of molecules selected. Assumes that the input INFILE
% contains an MxNx1 supercell. Also assumes that the crystal consists of
% one molecule with only one orientation. BURGERS should be a real number
% indicating how many burgers vectors the helix needs to rise. Right click
% or hit enter to stop selecting points.

% Load in molecules
mols=readMol2(infile,order);
numMols=size(mols,1);

% Represent each molecule by the mean of its atom positions
molCenters=zeros(numMols,3);
molCenters(:,:)=mean(mols,2);

% Figure out a normal to the plane
% Begin by testing to see if the first three points are colinear
firstPt=2;
secondPt=3;

% Keep testing poitns until we find three non-colinear points. Plane normal
% is then the cross product of any two unique vectors formed by the three
% poitns.
while 1
    planeNormal=cross(molCenters(firstPt,:)-molCenters(1,:),molCenters(secondPt,:)-molCenters(1,:));
    if norm(planeNormal) > 0.001
        break
    else
        secondPt=secondPt+1;
    end
    
    if secondPt>length(molCenters)
        error('Couldn''t find a plane.')
    end
end

%Normalize planeNormal and generate rotation matricies to rotate the
%molecules such that the supercell lies in the xy plane. This will make it
%easier to display to the user. (See 6/3/14 lab notes for how this is done.)
planeNormal=planeNormal/norm(planeNormal);
theta=acos(planeNormal(3));
phi=atan2(planeNormal(2),planeNormal(1));

Rz=[ cos(phi) sin(phi) 0;
    -sin(phi) cos(phi) 0;
        0        0     1;];
Ry=[cos(theta) 0 -sin(theta);
        0      1      0     ;
    sin(theta) 0  cos(theta);];

molCenters=(Ry*Rz*molCenters')';
molCenters=molCenters(:,1:2); %Get rid of the z coordinate since it's 0 for all molecules.

% Figure out the directions of the a and b vectors so we can generate the c
% vector (which is then the scale and direction of the burgers vector).
% First use a convex hull to find the outermost molecules of the supercell.
K=convhull(molCenters(:,1),molCenters(:,2));

% Then the first and second molecules in the convex hull should form the a
% direction and the first and second to last molecules should form the b
% direction
aVect=molCenters(K(2),:)-molCenters(K(1),:);
bVect=molCenters(K(end-1),:)-molCenters(K(1),:);

aVect=aVect/norm(aVect); bVect=bVect/norm(bVect);

% Calculate the c direction based on a, b, and alpha and beta. (Again, see
% 6/3/13 lab notes for how this is done.)
alpha=latticeParams(2,1); beta=latticeParams(2,2);

y0=(bVect(1)*cosd(beta)-aVect(1)*cosd(alpha))/(aVect(2)*bVect(1)-bVect(2)*aVect(1));
x0=(bVect(2)*cosd(beta)-aVect(2)*cosd(alpha))/(aVect(1)*bVect(2)-bVect(1)*aVect(2));

cVect=[x0 y0 sqrt(1-x0^2-y0^2)];

aVect=aVect*latticeParams(1,1);
bVect=bVect*latticeParams(1,2);
cVect=cVect*latticeParams(1,3);

% Plot ALL the molecules!
h=figure(1);
plot(molCenters(:,1), molCenters(:,2),'o','MarkerFaceColor','y')
hold on

plot([0, aVect(1)*5],[0, aVect(2)*5],[0, bVect(1)*5],[0,bVect(2)]*5,'LineWidth',5)
titles={'','a','b'};
legend(titles)

pointsStack=zeros(1);
counter=1;

% Make sure the figure is open and begin polling for inputs.
while ishandle(h)
    [x,y,button]=ginput(1);
    % If the user left clicked, find the closest point and enter it into
    % the stack, making sure it's not in there already.
    if ~isempty(button) && button==1
        [~,clickedPoint]=min((molCenters(:,1)-x).^2+(molCenters(:,2)-y).^2);
        if ~sum(find(pointsStack==clickedPoint))
            pointsStack(counter)=clickedPoint;
            plot(molCenters(pointsStack(counter),1),molCenters(pointsStack(counter),2),'o','MarkerFaceColor','g');
            
            counter=counter+1;
            drawnow
        end
    else
        break
    end
end
hold off
close(1)

if pointsStack(1)
    numSelected = length(pointsStack);
    % Find the direction of the burgers vector. This is done by rotating
    % the cVector we calculated by the inverse of the rotation that placed
    % the molecule plane on the xy plane.
    realcVect=(Ry*Rz)\cVect';
    dispIncrement = realcVect*burgers/numSelected;
    
    % Raise all the atoms in each molecule by the appropriate amount.
    for k=1:numSelected
        mols(pointsStack(k),:,1)=mols(pointsStack(k),:,1)+k*dispIncrement(1);
        mols(pointsStack(k),:,2)=mols(pointsStack(k),:,2)+k*dispIncrement(2);
        mols(pointsStack(k),:,3)=mols(pointsStack(k),:,3)+k*dispIncrement(3);
    end
    
    tempMols=mols;
    
    % Output as a mols array
    mols=cell(numMols,1);
    for mol=1:numMols
        mols{mol}=zeros(size(tempMols,2),3);
        mols{mol}(:,:)=tempMols(mol,:,:);
    end
else
    mols={};
end
