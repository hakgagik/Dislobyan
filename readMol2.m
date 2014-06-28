function mols=readMol2(infile,order)

if ~strcmp(infile(end-4:end),'.mol2')
    error('Please select a .mol2 file.');
end

ifid = fopen(infile);

if ifid <2
    error('Input file does not exist in current directory')
end

rawData=textscan(ifid,'%d %s %f %f %f %s %s %s %s','HeaderLines',9);
fclose(ifid);
xin=rawData{3}(1:end);
yin=rawData{4}(1:end);
zin=rawData{5}(1:end);

molsPerCell=size(order,1);
atomsPerMol=size(order,2);

atomsPerCell=molsPerCell*atomsPerMol;

numAtoms=length(xin);
numCells=numAtoms/atomsPerCell;
numMols=numCells*molsPerCell;

if round(numMols)~=numMols
    error('Order matrix does not match the total number of atoms.')
end

mols = zeros(numMols,atomsPerMol,3);

for mol=1:numMols
    index=atomsPerMol*(mol-1)+1;
    endIndex=index+atomsPerMol-1;
    atomsBuffer=[xin(index:endIndex) yin(index:endIndex) zin(index:endIndex)];
    mols(mol,:,:)=atomsBuffer(order,:);
end

