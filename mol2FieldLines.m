function [pointsStack]=mol2FieldLines(mols,startPoint,plot)

forceMultiplier=1;
maxPoints=15000;

[centers, directions]=mol2VectorField(mols,0);



bubbleRadius = 8;
sigmaSqr=(bubbleRadius/2)^2;
a=1/sqrt(2*pi*sigmaSqr);
currentPoint=startPoint;
pointsStack=zeros(maxPoints,3);
pointsStack(1,:)=currentPoint;
counter=1;

for k=2:maxPoints
    pointDisps=bsxfun(@minus,centers,currentPoint);
    dist2CurrentPoint=(pointDisps(:,1).^2+pointDisps(:,2).^2+pointDisps(:,3).^2).^.5;
    nearbyPoints=find(dist2CurrentPoint<bubbleRadius);
    
    if isempty(nearbyPoints)
        break
    end
    
    pointWeights=a*exp(-dist2CurrentPoint(nearbyPoints).^2/(2*sigmaSqr));
    force=sum(bsxfun(@times, directions(nearbyPoints,:), pointWeights))*forceMultiplier;
    currentPoint=currentPoint+force;
    pointsStack(k,:)=currentPoint;
    if counter>=maxPoints
        break
    end
    counter=counter+1;
end

pointsStack=pointsStack(1:counter,:);

currentPoint=startPoint;
backwardPointsStack=zeros(maxPoints,3);

counter=0;

for k=1:maxPoints
    pointDisps=bsxfun(@minus,centers,currentPoint);
    dist2CurrentPoint=(pointDisps(:,1).^2+pointDisps(:,2).^2+pointDisps(:,3).^2).^.5;
    nearbyPoints=find(dist2CurrentPoint<bubbleRadius);
    
    if isempty(nearbyPoints)
        break
    end
    
    pointWeights=a*exp(-dist2CurrentPoint(nearbyPoints).^2/(2*sigmaSqr));
    
    force=sum(bsxfun(@times, directions(nearbyPoints,:), pointWeights))/sum(pointWeights)*forceMultiplier;
    currentPoint=currentPoint-force;
    backwardPointsStack(k,:)=currentPoint;
    if counter>=maxPoints
        break
    end
    counter=counter+1;
end
backwardPointsStack=backwardPointsStack(1:counter,:);

pointsStack=[flip(backwardPointsStack); pointsStack];

if plot
    whatToPlot=mod(1:size(pointsStack,1),10)==1;
    plot3(pointsStack(whatToPlot,1), pointsStack(whatToPlot,2),pointsStack(whatToPlot,3))
end