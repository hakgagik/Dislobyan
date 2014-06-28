function [a,b,c]=latticeFT(centers, limits, Depth, maxDepth)

xSubDiv=10;
ySubDiv=10;
zSubDiv=10;

lMin=limits(1);
lMax=limits(2);
mMin=limits(3);
mMax=limits(4);
nMin=limits(5);
nMax=limits(6);

lSteps=(lMax-lMin)/(xSubDiv-1);
mSteps=(mMax-mMin)/(ySubDiv-1);
nSteps=(nMax-nMin)/(zSubDiv-1);

x=centers(:,1);
y=centers(:,2);
z=centers(:,3);

M=zeros(xSubDiv,ySubDiv,zSubDiv);
counter=0;

for l=lMin:lSteps:lMax
    for m=mMin:mSteps:mMax
        for n=nMin:nSteps:nMax
            counter=counter+1;
            [Ml,Mm,Mn]=ind2sub([xSubDiv ySubDiv zSubDiv],counter);
            M(Ml,Mm,Mn)=abs(sum(exp(-2*pi*1i*(l*x+m*y+n*z))));
        end
    end
end

lVect=lMin:lSteps:lMax;
mVect=mMin:mSteps:mMax;
nVect=nMin:nSteps:nMax;


newLimits=[0 0 0 0 0 0]; %Placeholder

if Depth<maxDepth
    [a,b,c]=latticeFT(centers,newLimits,Depth+1,maxDepth);
end