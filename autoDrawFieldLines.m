function autoDrawFieldLines(infile,order)

mols=readMol2(infile,order);
figure('units','normalized','outerposition',[0 0 1 1])
[centers, ~]=mol2VectorField(mols,1);
hold on
axis equal
view(-90,90);
drawnow
skip=2;
for k=2:skip:length(centers)
    [~]=mol2FieldLines(mols, centers(k,:),1);
    fprintf('%d\n',k*100/length(centers));
end
drawnow
figure(1);

saveas(gcf,infile(1:end-5),'png');

close(1);

