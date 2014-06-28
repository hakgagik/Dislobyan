function moveThroughField(mols,startPoint)

[centers, directions]=mol2VectorField(mols,0);

currentPoint=startPoint;
h=plot3(currentPoint(1),currentPoint(2),currentPoint(3),'o','MarkerFaceColor','g');
drawnow
forceMultiplier=1000;
pause

bubbleRadius=8;
sigmaSqr=(bubbleRadius/2)^2;
a=1/sqrt(2*pi*sigmaSqr);
counter=0;
while 1
    pointDisps=bsxfun(@minus,centers,currentPoint);
    dist2CurrentPoint=(pointDisps(:,1).^2+pointDisps(:,2).^2+pointDisps(:,3).^2).^.5;
    nearbyPoints=find(dist2CurrentPoint<bubbleRadius);
    
    if isempty(nearbyPoints)
        break
    end
    pointWeights=a*exp(-dist2CurrentPoint(nearbyPoints).^2/(2*sigmaSqr));
    force=sum(bsxfun(@times, directions(nearbyPoints,:), pointWeights))*forceMultiplier;
    currentPoint=currentPoint+force;
    if ~mod(counter,1000)
        delete(h)
        h=plot3(currentPoint(1),currentPoint(2),currentPoint(3),'o','MarkerFaceColor','g');
        drawnow
    end
    counter=counter+1;
end