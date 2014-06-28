function [centers, directions]=mol2VectorField(mols,plot)

if ~iscell(mols)
    numMols=size(mols,1);
    tempMols=mols;
    mols=cell(numMols,1);
    for mol=1:numMols
        mols{mol}=zeros(size(tempMols,2),3);
        mols{mol}(:,:)=tempMols(mol,:,:);
    end
end

numMols=length(mols);

centers=zeros(numMols,3);
directions=zeros(numMols,3);

% Find the mean direction of each molecule
for k=1:numMols
    tempMol=mols{k};
    centers(k,:)=mean(tempMol);
    [~,~,V]=svd(bsxfun(@minus,tempMol,centers(k,:)));
    directions(k,:)=V(:,1);
   if dot(directions(k,:),(mols{k}(end,:)-mols{k}(1,:))/norm(mols{k}(end,:)-mols{k}(1,:)))<0;
       directions(k,:)=-directions(k,:);
   end
end

% ----------TESTING! COMMENT OUT THIS LINE WHEN NOT NECESSARY!----------- %
directions=bsxfun(@minus,directions,mean(directions));
directionNorms=sqrt(sum(directions.^2,2));
directions(abs(directionNorms-mean(directionNorms))>10*std(directionNorms),:)=zeros(size(directions(abs(directionNorms-mean(directionNorms))>10*std(directionNorms),:))); % Supress outliers
% ----------------------------------------------------------------------- %


% Plot said directions as vectors.
if plot
    centers=centers-directions/2;
    quiver3(centers(:,1),centers(:,2),centers(:,3),directions(:,1),directions(:,2),directions(:,3))
end
