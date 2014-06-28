% h=figure;
% [centers,directions]=mol2VectorField(mols,1);
% hold on
% axis equal
% drawnow
% 
% pause

% for k=0:10:110;
%     for j=0:10:150;
%         mol2FieldLines(mols,[k, j, 18])
%     end
%     drawnow
% end

skip=2;
for k=2:skip:length(centers)
    [~]=mol2FieldLines(mols, centers(k,:),1);
    drawnow
end

figure(h)

% directionNorms=sqrt(sum(directions.^2,2));
% 
% scatter3(centers(:,1),centers(:,2),centers(:,3),50,directionNorms,'filled');