function varargout = dislocateGUI(varargin)
% DISLOCATEGUI MATLAB code for dislocateGUI.fig
%      DISLOCATEGUI, by itself, creates a new DISLOCATEGUI or raises the existing
%      singleton*.
%
%      H = DISLOCATEGUI returns the handle to a new DISLOCATEGUI or the handle to
%      the existing singleton*.
%
%      DISLOCATEGUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DISLOCATEGUI.M with the given input arguments.
%
%      DISLOCATEGUI('Property','Value',...) creates a new DISLOCATEGUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before dislocateGUI_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to dislocateGUI_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help dislocateGUI

% Last Modified by GUIDE v2.5 19-Jul-2013 16:09:54

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @dislocateGUI_OpeningFcn, ...
    'gui_OutputFcn',  @dislocateGUI_OutputFcn, ...
    'gui_LayoutFcn',  [] , ...
    'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before dislocateGUI is made visible.
function dislocateGUI_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to dislocateGUI (see VARARGIN)

% Choose default command line output for dislocateGUI
handles.output = hObject;

plot1=0;
while ishandle(plot1)
    plot1=plot1+1;
end
plot2=plot1+1;
while ishandle(plot2)
    plot2=plot2+1;
end
handles.plot1=plot1;
handles.plot2=plot2;

handles.rotation=eye(3);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes dislocateGUI wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = dislocateGUI_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function infileEdit_Callback(hObject, eventdata, handles) %#ok<*INUSL,*DEFNU>
% hObject    handle to infileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of infileEdit as text
%        str2double(get(hObject,'String')) returns contents of infileEdit as a double
infile=get(hObject,'String');
%If the infile contains '.mat', enable variable inporting options.
if length(infile)>3
    if strcmp(infile(end-3:end),'.mat')
        set(handles.structureCheckbox,'Enable','on')
        set(handles.paramsCheckbox,'Enable','on')
        set(handles.formatCheckbox,'Enable','on')
        set(handles.molsCheckbox,'Enable','on')
    else
        set(handles.structureCheckbox,'Enable','off')
        set(handles.paramsCheckbox,'Enable','off')
        set(handles.formatCheckbox,'Enable','off')
        set(handles.molsCheckbox,'Enable','off')
    end
end


% --- Executes during object creation, after setting all properties.
function infileEdit_CreateFcn(hObject, eventdata, handles) %#ok<*INUSD>
% hObject    handle to infileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in loadButton.
function loadButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.errorBox,'String', '')
infile=get(handles.infileEdit,'String');
if length(infile)>4
    if strcmp(infile(end-3:end),'.mat') %Call the inportVariables function if the infile has a '.mat' in it.
        handles=importVariables(infile,handles);
    else
        %Import the supercell size if the infile contains it
        scStart=regexp(infile,'\d\d\d\d\d\d');
        if scStart
            handles.supercell=[str2double(infile(scStart:(scStart+1))), ...
                str2double(infile((scStart+2):scStart+3)), ...
                str2double(infile((scStart+4):(scStart+5)))];
            set(handles.supercellTable,'Data',handles.supercell);
        end
        orderPresent=1;
        if ~get(handles.covalentCheckBox,'Value')
            try
                %Make sure order has been defined before attempting to sort
                handles.order;
            catch %#ok<*CTCH>
                set(handles.errorBox,'String','Order array is missing, please import or create it.');
                orderPresent=0;
            end
        end
        if orderPresent
            if get(handles.storeSlot2Radio, 'Value')==1 %Check which slot the user wants mols saved in
                [handles, handles.mols2]=readFile(infile, handles);
            else
                [handles, handles.mols1]=readFile(infile, handles);
            end
        else
            set(handles.errorBox,'String','Check the command window!')
            atomsPerCell=input('How many atoms in a single cell?\n');
            [order,types,bonds]=readMols2(infile, atomsPerCell);
            handles.order=order;
            maxLength=length(order{1});
            orderLength=length(order);
            for n=2:orderLength
                if length(order{n})>maxLength
                    maxLength=length(order{n});
                end
            end
            paddedOrder=zeros(orderLength,maxLength);
            for n=1:orderLength
                for m=1:length(order{n})
                    paddedOrder(n,m)=order{n}(m);
                end
            end
            set(handles.orderTable,'Data',paddedOrder)
            set(handles.orderTable,'ColumnEditable',true(1,size(paddedOrder,2)))
            
            handles.types=types;
            set(handles.typesTable,'Data',handles.types)
            set(handles.typesTable,'ColumnEditable',true(1,length(handles.types)))
            
            handles.bonds=bonds;
            set(handles.bondsTable,'Data',handles.bonds)
            set(handles.bondsTable,'ColumnEditable',true(1,3))
        end
    end
    guidata(handles.figure1,handles)
else
    set(handles.errorBox,'String','Filename too short, don''t forget to include a file extention.')
end

function [order, types, bonds]=readMols2(infile,atomsPerCell)

ifid=fopen(infile,'r');
lineN=0;

while 1
    tline=fgetl(ifid);
    lineN=lineN+1;
    if strfind(tline, '@<TRIPOS>ATOM')
        atomsLine=lineN;
        break
    end
end

while 1
    tline=fgetl(ifid);
    lineN=lineN+1;
    if strfind(tline, '@<TRIPOS>BOND')
        bondsLine=lineN;
        break
    end
end

fclose(ifid);

ifid=fopen(infile,'r');
rawTypes=textscan(ifid, '%*d %s %*f %*f %*f %*s %*d %*s %*d','HeaderLines',atomsLine);
rawTypes=rawTypes{1};
rawTypes=rawTypes(1:atomsPerCell);
fclose(ifid);

ifid=fopen(infile,'r');
rawBonds=textscan(ifid, '%*d %d %d %d', 'HeaderLines',bondsLine);
fclose(ifid);
bonds=[rawBonds{1} rawBonds{2} rawBonds{3}];
[lastBond,~]=find(bonds(:,1:2)>atomsPerCell,1);

firstBonds=bonds(1:lastBond-1,1:2);

todo=java.util.LinkedList;
orderQueue=java.util.LinkedList;
while nnz(firstBonds)
    tempMol=java.util.LinkedList;
    [tempBond,~]=find(firstBonds>0,1);
    tempMol.add(firstBonds(tempBond,1));
    tempMol.add(firstBonds(tempBond,2));
    todo.add(firstBonds(tempBond,1));
    todo.add(firstBonds(tempBond,2));
    firstBonds(tempBond,:)=[0 0];
    while todo.size()
        currentMol=todo.remove();
        %Bonds Containing This Molecule
        [BCTM,col]=find(firstBonds==currentMol);
        for k=1:length(BCTM)
            todo.add(firstBonds(BCTM(k),-col(k)+3));
            if ~tempMol.contains(firstBonds(BCTM(k), -col(k)+3))
                tempMol.add(firstBonds(BCTM(k),-col(k)+3));
            end
            firstBonds(BCTM(k),:)=[0 0];
        end
    end
    orderQueue.add(tempMol);
end

order=cell(orderQueue.size(),1);
index1=0;
while orderQueue.size()
    index1=index1+1;
    tempMol=orderQueue.remove();
    order{index1}=zeros(tempMol.size(),1);
    index2=0;
    while tempMol.size()
        index2=index2+1;
        order{index1}(index2)=tempMol.remove;
    end
end

for k=1:length(rawTypes)
    rawTypes{k}=rawTypes{k}(rawTypes{k}>57);
end
types=rawTypes(order{1});
bonds=bonds(order{1});


% Reads a mat file for relevant information and writes it all into the
% handles structure based on user preferences.
function handles = importVariables(infile, handles)
fid=fopen(infile);
if infile>2
    fclose(fid);
    tempVars=load(infile);
    varNames=fieldnames(tempVars);
    numVars=length(varNames);
    % Run through the newly imported variables, look for keywords, and
    % import anything that matches them.
    if get(handles.structureCheckbox, 'Value') %Structure parameters
        for k=1:numVars
            if ~isempty(strfind(varNames{k},'param')) || ~isempty(strfind(varNames{k},'attice'))
                %Make sure the latticeParams array has the correct shape.
                if size(tempVars.(varNames{k}),1)==3
                    handles.latticeParams=tempVars.(varNames{k});
                    set(handles.paramsTable, 'Data', handles.latticeParams');
                elseif size(tempVars.(varNames{k}),1)==2
                    handles.latticeParams=(tempVars.(varNames{k}))';
                    set(handles.paramsTable, 'Data', handles.latticeParams');
                end
            end
            if ~isempty(strfind(varNames{k},'upercell')) || ~isempty(strfind(varNames{k}, 'uperCell'))
                supercell=[0 0 0];
                supercell(:)=tempVars.(varNames{k});
                set(handles.supercellTable, 'Data', supercell);
                handles.supercell=supercell;
            end
            if strfind(varNames{k},'otation')
                handles.rotation=tempVars.(varNames{k});
                set(handles.rotationTable, 'Data', handles.rotation);
            end
        end
    end
    if get(handles.paramsCheckbox, 'Value') %Dislocation parameters
        for k=1:numVars
            if strfind(varNames{k},'ocation')
                location=[0 0 0];
                location(:)=tempVars.(varNames{k});
                set(handles.locationTable,'Data',location)
                handles.location=location;
            elseif strfind(varNames{k}, 'ense')
                sense=[0 0 0];
                sense(:)=tempVars.(varNames{k});
                set(handles.senseTable,'Data',sense);
                handles.sense=sense;
            elseif strfind(varNames{k}, 'urgers')
                burgers=[0 0 0];
                burgers(:)=tempVars.(varNames{k});
                set(handles.burgersTable,'Data',burgers);
                handles.burgers=burgers;
            elseif strfind(varNames{k}, 'rientation')
                orientation=[0 0 0];
                orientation(:)=tempVars.(varNames{k});
                set(handles.orientationTable,'Data',orientation);
                handles.orientation=orientation;
            end
        end
    end
    if get(handles.formatCheckbox, 'Value') %Format parameters
        for k=1:numVars
            if strfind(varNames{k},'rder')
                order=tempVars.(varNames{k});
                % Make sure that order is a cell for storage and an array for display
                if iscell(order)
                    handles.order=order;
                    maxLength=length(order{1});
                    orderLength=length(order);
                    for n=2:orderLength
                        if length(order{n})>maxLength
                            maxLength=length(order{n});
                        end
                    end
                    paddedOrder=zeros(orderLength,maxLength);
                    for n=1:orderLength
                        for m=1:length(order{n})
                            paddedOrder(n,m)=order{n}(m);
                        end
                    end
                    set(handles.orderTable,'Data',paddedOrder)
                    set(handles.orderTable,'ColumnEditable',true(1,size(paddedOrder,2)))
                elseif isnumeric(order)
                    set(handles.orderTable,'Data',order)
                    set(handles.orderTable,'ColumnEditable',true(1,size(order,2)))
                    tempOrder=order;
                    order=cell(size(tempOrder,1),1);
                    for n=1:size(tempOrder,1)
                        order{n}=tempOrder(n,:);
                    end
                    handles.order=order;
                end
            elseif strfind(varNames{k}, 'ypes') %Make sure types is a cell
                types=tempVars.(varNames{k});
                if ischar(types)
                    tempTypes=types;
                    types=cell(size(tempTypes,1),1);
                    for n=1:size(types,1)
                        types{n}=tempTypes(n,:);
                    end
                end
                if size(types,1)>size(types,2)
                    types=types';
                end
                handles.types=types;
                set(handles.typesTable,'Data',handles.types)
                set(handles.typesTable,'ColumnEditable',true(1,length(handles.types)))
            elseif strfind(varNames{k}, 'onds')
                handles.bonds=tempVars.(varNames{k});
                set(handles.bondsTable,'Data',handles.bonds)
                set(handles.bondsTable,'ColumnEditable',true(1,3))
            end
        end
    end
    if get(handles.molsCheckbox, 'Value') % Mols array, this will probably not be used much
        for k=1:numVars
            if strfind(varNames{k}, 'mol')
                if get(handles.storeSlot2Radio, 'Value')
                    handles.mols2=tempVars.(varNames{k});
                    plotMols(handles,2);
                else
                    handles.mols1=tempVars.(varNames{k});
                    plotMols(handles,1);
                end
            end
        end
    end
    set(handles.errorBox, 'String', 'Anything successfully imported will be reflected in the above fields')
    guidata(handles.figure1,handles)
else
    set(handles.errorBox,'String','No such file in directory.')
end


% TODO: Implement this function.
function plotMols(handles,plotFrom)
if plotFrom==1
    mols=handles.mols1;
    figure(handles.plot1)
elseif plotFrom==2
    mols=handles.mols2;
    figure(handles.plot2)
else
    error('You somehow referenced an unknown slot.')
end
numMols=length(mols);
positions=zeros(numMols,3);
for m=1:numMols
    positions(m,:)=mols{m}(1,:);
end
plot3(positions(:,1),positions(:,2),positions(:,3),'bs')
axis equal
try
    for m=1:numMols
        positions(m,:)=mols{m}(2,:);
    end
    hold on
    plot3(positions(:,1),positions(:,2),positions(:,3),'gs')
    hold off
catch
end




% Reads a mol2 or pdb file and returns an ordered "mols" array
function [handles, mols] = readFile(infile, handles)
% infile     The input file
% order      Order of the atoms in each cell. This should either be an
% array in the form [mol1; mol2; ...] or a cell in the form {mol1, mol2,
% ...], where mol1, mol2, ... are column vectors which designate the atoms
% in the first cell of the input file that belong to their respective
% molecules. I. e. [7 8 1 2 3; 9 10 4 5 6] means that the 7th, 8th,
% 1st, 2nd, and 3rd atoms in the infile are members of the first
% molecule, while the 9th, 10th, 4th, 5th, and 6th atoms are members of the
% second molecule. The ordering of each row, in most cases, is irrelevant.

ifid=fopen(infile);
goodFile=1;
if ifid<=2
    set(handles.errorBox, 'String', 'File does not exist in current directory.')
    mols={};
else %Just a bunch of reading through the infile...
    if strcmp(infile(end-2:end),'pdb')
        set(handles.errorBox,'String','Reading .pdb file...')
        drawnow
        header=textscan(ifid,'%s',50);
        header=header{1};
        fclose(ifid);
        if strcmp(header(16),'CRYST1') %see if we can import lattice params
            headerlines=3;
            try
                supercell=handles.supercell;
                latticeParams=zeros(3,2);
                for param=1:6
                    latticeParams(param)=str2double(header(16+param));
                end
                latticeParams(1:3)=latticeParams(1:3)./supercell;
                handles.latticeParams=latticeParams;
                set(handles.paramsTable,'Data',latticeParams')
                guidata(handles.figure1,handles)
                drawnow
            catch
                set(handles.errorBox,'String','Couldn''t import lattice parameters since no supercell array is present.')
            end
        else
            headerlines=2;
        end
        if strcmp(header{20},'C') || strcmp(header{28},'C')
            weirdFile=1;
        else
            weirdFile=0;
        end
        ifid=fopen(infile);
        if weirdFile
            rawData=textscan(ifid,'%s %f %s %s %*s %f %f %f %f %f %f %s','HeaderLines',headerlines);
        else
            rawData=textscan(ifid,'%s %f %s %s %f %f %f %f %f %f %s','HeaderLines',headerlines);
        end
        fclose(ifid);
        xin=rawData{6}(1:end-1);
        yin=rawData{7}(1:end-1);
        zin=rawData{8}(1:end-1);
        
    elseif strcmp(infile(end-3:end),'mol2') %Same comments as above
        set(handles.errorBox,'String','Reading .mol2 file...')
        drawnow
        header=textscan(ifid,'%s');
        header=header{1};
        fclose(ifid);
        if strcmp(header(end-8),'@<TRIPOS>CRYSIN')
            try
                supercell=handles.supercell;
                latticeParams=zeros(3,2);
                for param=1:6
                    latticeParams(param)=str2double(header(end-8+param));
                end
                latticeParams(1:3)=latticeParams(1:3)./supercell;
                handles.latticeParams=latticeParams;
                set(handles.paramsTable,'Data',latticeParams')
                guidata(handles.figure1,handles)
            catch
                set(handles.errorBox,'String','Couldn''t import lattice parameters since no supercell array is present.')
            end
        end
        headerlines=9;
        ifid=fopen(infile);
        rawData=textscan(ifid,'%d %s %f %f %f %s %s %s %s','HeaderLines',headerlines);
        fclose(ifid);
        xin=rawData{3}(1:end);
        yin=rawData{4}(1:end);
        zin=rawData{5}(1:end);
    else
        set(handles.errorBox, 'String', 'Unrecognized or unsupported file format.')
        goodFile=0;
    end
    %Begin sorting the raw data
    mols={};
    if goodFile
        order=handles.order;
        molsPerCell=length(order);
        atomsPerMol=zeros(1,molsPerCell);
        for k=1:molsPerCell
            atomsPerMol(k)=length(order{k});
        end
        
        atomsPerCell=sum(atomsPerMol);
        numAtoms=length(xin);
        handles.numAtoms=numAtoms;
        numCells=numAtoms/atomsPerCell;
        numMols=numCells*molsPerCell;
        
        if round(numMols)~=numMols
            set(handles.errorBox,'String','Order matrix does not match the total number of atoms.')
            mols={};
        else
            mols=cell(numMols,1);
            for cel=1:numCells
                index=atomsPerCell*(cel-1)+1;
                endIndex=index+atomsPerCell-1;
                atomsBuffer=[xin(index:endIndex) yin(index:endIndex) zin(index:endIndex)];
                for m=1:molsPerCell
                    mols{molsPerCell*(cel-1)+m}=atomsBuffer(order{m},:);
                end
            end
            set(handles.errorBox,'String','Successfully read and sorted molecules.')
        end
    end
end


% --- Executes on button press in loadStructureButton.
function loadStructureButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadStructureButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.errorBox,'String','')
handles.latticeParams=get(handles.paramsTable,'Data');
handles.supercell=get(handles.supercellTable,'Data');
set(handles.errorBox,'String','Successfully updated supercell structure parameters.');
guidata(handles.figure1,handles);



% --- Executes on button press in structureCheckbox.
function structureCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to structureCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of structureCheckbox


% --- Executes on button press in paramsCheckbox.
function paramsCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to paramsCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of paramsCheckbox


% --- Executes on button press in formatCheckbox.
function formatCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to formatCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of formatCheckbox


% --- Executes on button press in molsCheckbox.
function molsCheckbox_Callback(hObject, eventdata, handles)
% hObject    handle to molsCheckbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of molsCheckbox


% --- Executes on button press in exportButton.
function exportButton_Callback(hObject, eventdata, handles)
% hObject    handle to exportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.errorBox,'String','')
outfile=get(handles.outfileEdit,'String');

if get(handles.xsdRadio, 'Value')
    extension='.xsd';
elseif get(handles.molsRadio, 'Value')
    extension='.mat';
else
    extention='.pdb';
end

if isempty(strfind(outfile,extension))
    outfile=[outfile extension];
end

if get(handles.xsdRadio, 'Value')
    writeXSD(handles,outfile)
elseif get(handles.molsRadio, 'Value')
    writeMols(handles,outfile)
else
    writePDB(handles,outfile)
end


function writeXSD(handles, outfile)

varsLoaded=1;
try
    bonds=handles.bonds;
    order=handles.order;
    types=handles.types;
    if get(handles.exportSlot1Radio,'Value')
        mols=handles.mols1;
    else
        mols=handles.mols2;
    end
catch
    varsLoaded=0;
    
end

if varsLoaded
    set(handles.errorBox,'String','Started writing...')
    drawnow
    if get(handles.covalentCheckBox,'Value')
        molsPerCell=length(mols);
    else
        molsPerCell=length(order);
    end
    if get(handles.covalentCheckBox,'Value')
        atomsPerMol=length(mols);
    else
        atomsPerMol=zeros(1,molsPerCell);
        for k=1:molsPerCell
            atomsPerMol(k)=size(order{k},2);
        end
    end
    
    atomsPerCell=sum(atomsPerMol);
    numMols=length(mols);
    numCells=numMols/molsPerCell;
    
    if get(handles.covalentCheckBox,'Value')
        dislocatedMols=zeros(numCells,atomsPerCell,3);
        for m=1:numMols
            dislocatedMols(1,m,:)=mols{m};
        end
    elseif length(order)==1 || length(handles.order{1})==length(handles.order{2})
        dislocatedMols=zeros(length(mols),length(mols{1}),3);
        for m=1:length(mols)
            dislocatedMols(m,:,:)=mols{m};
        end
    else
        dislocatedMols=zeros(numCells,atomsPerCell,3);
        for cel=1:numCells
            for m=1:molsPerCell
                dislocatedMols(cel,order{m},:)=mols{(cel-1)*molsPerCell+m};
            end
        end
    end
    
    mols=dislocatedMols;
    
    numAtoms=numel(mols)/3;
    if get(handles.covalentCheckBox,'Value')
        numMolecules=1;
    else
        numMolecules=size(mols,1);
    end
    atomsPerMol=size(mols,2);
    bondsPerMol=size(bonds,1);
    numBonds=bondsPerMol*numMolecules;
    atomIDs=5:(numAtoms+5);
    bondIDs=(numAtoms+6):(numAtoms+numBonds+5);
    
    connects=zeros(numBonds,3);
    if get(handles.covalentCheckBox,'Value')
        connects=[bonds(:,1:2)+4 bonds(:,3)];
    else
        for m=1:numMols
            milestone=0;
            connects(((m-1)*bondsPerMol+1):m*bondsPerMol,:)=[bonds(:,1:2)+(m-1)*atomsPerMol+4, bonds(:,3)];
            if m/numMols*100>milestone
                set(handles.errorBox,'String', ['Stage 1 of 4: ' num2str(round(m/numMols*100)) '% done...'])
                drawnow
                milestone=milestone+1;
            end
        end
    end
    
    single_connections=cell(atomsPerMol,1);
    if get(handles.covalentCheckBox,'Value')
        atomsPerMol=numAtoms;
    end
    for k=1:atomsPerMol
        [connections_buffer,~]=find(bonds(:,1:2)==k);
        single_connections{k}=connections_buffer;
    end
    
    connections=cell(numAtoms,1);
    numAtomsPlusFive=numAtoms+5;
    for m=1:numMolecules
        for at=1:atomsPerMol
            connections{(m-1)*atomsPerMol+at}=regexprep(num2str((single_connections{at}+bondsPerMol*(m-1)+numAtomsPlusFive)'),' +',',');
        end
        milestone=0;
        if m/numMols*100>milestone
            set(handles.errorBox,'String',['Stage 2 of 4: ' num2str(round(m/numMolecules*100)) '% done...'])
            drawnow
            milestone=milestone+1;
        end
    end
    
    
    names=cell(atomsPerMol,1);
    present_atoms=unique(types);
    num_unique=length(present_atoms);
    num_present=zeros(num_unique,1);
    if get(handles.covalentCheckBox,'Value')
        for k=1:numAtoms
            names{k}=[types{1} num2str(k)];
        end
    else
        for at=1:atomsPerMol
            for k=1:num_unique
                if strcmp(present_atoms{k},types{at})
                    num_present(k)=num_present(k)+1;
                    names{at}=[types{at} num2str(num_present(k))];
                    break
                end
            end
        end
    end
    if get(handles.covalentCheckBox,'Value')
        elemType=types{1};
        types=cell(1,numMols);
        for m=1:numMols
            types{m}=elemType;
        end
    end
    fid=fopen(outfile,'w+');
    fprintf(fid, '<?xml version="1.0" encoding="latin1"?>\r\n<!DOCTYPE XSD []>\r\n<XSD Version="4.2">\r\n');
    fprintf(fid, '\t<AtomisticTreeRoot ID="1" NumProperties="47" NumChildren="1">\r\n');
    fprintf(fid, '\t\t<Property Name="AngleEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="BendBendEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="BendTorsionBendEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="BondEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="EFGAsymmetry" DefinedOn="Atom" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="EFGQuadrupolarCoupling" DefinedOn="Atom" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="ElectrostaticEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="FaceMillerIndex" DefinedOn="GrowthFace" Type="MillerIndex"/>\r\n');
    fprintf(fid, '\t\t<Property Name="FacetTransparency" DefinedOn="GrowthFace" Type="Float"/>\r\n');
    fprintf(fid, '\t\t<Property Name="Force" DefinedOn="Bondable" Type="CoDirection"/>\r\n');
    fprintf(fid, '\t\t<Property Name="FrameFilter" DefinedOn="Trajectory" Type="String"/>\r\n');
    fprintf(fid, '\t\t<Property Name="HarmonicForceConstant" DefinedOn="HarmonicRestraint" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="HarmonicMinimum" DefinedOn="HarmonicRestraint" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="HydrogenBondEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="ImportOrder" DefinedOn="Bondable" Type="UnsignedInteger"/>\r\n');
    fprintf(fid, '\t\t<Property Name="InversionEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="IsBackboneAtom" DefinedOn="Atom" Type="Boolean"/>\r\n');
    fprintf(fid, '\t\t<Property Name="IsChiralCenter" DefinedOn="Atom" Type="Boolean"/>\r\n');
    fprintf(fid, '\t\t<Property Name="IsOutOfPlane" DefinedOn="Atom" Type="Boolean"/>\r\n');
    fprintf(fid, '\t\t<Property Name="LineExtentPadding" DefinedOn="BestFitLineMonitor" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="LinkageGroupName" DefinedOn="Linkage" Type="String"/>\r\n');
    fprintf(fid, '\t\t<Property Name="ListIdentifier" DefinedOn="PropertyList" Type="String"/>\r\n');
    fprintf(fid, '\t\t<Property Name="NMRShielding" DefinedOn="Atom" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="NonBondEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="NormalMode" DefinedOn="Bondable" Type="Direction"/>\r\n');
    fprintf(fid, '\t\t<Property Name="NormalModeFrequency" DefinedOn="Bondable" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="NumScanSteps" DefinedOn="LinearScan" Type="UnsignedInteger"/>\r\n');
    fprintf(fid, '\t\t<Property Name="OrbitalCutoffRadius" DefinedOn="Bondable" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="PlaneExtentPadding" DefinedOn="BestFitPlaneMonitor" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="PotentialEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="QuantizationValue" DefinedOn="ScalarFieldBase" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="RestraintEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="ScanEnd" DefinedOn="LinearScan" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="ScanStart" DefinedOn="LinearScan" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="SeparatedStretchStretchEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="SimulationStep" DefinedOn="Trajectory" Type="Integer"/>\r\n');
    fprintf(fid, '\t\t<Property Name="StretchBendStretchEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="StretchStretchEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="StretchTorsionStretchEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="TorsionBendBendEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="TorsionEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="TorsionStretchEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="ValenceCrossTermEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="ValenceDiagonalEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="VanDerWaalsEnergy" DefinedOn="ClassicalEnergyHolder" Type="Double"/>\r\n');
    fprintf(fid, '\t\t<Property Name="Velocity" DefinedOn="Bondable" Type="Direction"/>\r\n');
    fprintf(fid, '\t\t<Property Name="_Stress" DefinedOn="SymmetrySystem" Type="Matrix"/>\r\n');
    fprintf(fid, '\t\t<Molecule ID="2" NumChildren="1" Name="%s" XYZ="0,0,0">\r\n', outfile(1:end-4));
    fprintf(fid, '\t\t\t<Chain ID="3" NumChildren="1" XYZ="0,0,0">\r\n');
    fprintf(fid, '\t\t\t\t<SubUnit ID="4" NumChildren="%d" Name="MOL2" UserID="2" XYZ="0,0,0">\r\n',round(numBonds+numAtoms));
    
    for m=1:numMolecules
        for at=1:atomsPerMol
            k=(m-1)*atomsPerMol+at;
            fprintf(fid, '\t\t\t\t\t<Atom3d ID="%d" Name="%s" UserID="%d" XYZ="%.14f,%.14f,%.14f" Connections="%s" TemperatureType="Isotropic" Components="%s"/>\r\n', ...
                atomIDs(k), names{at}, k, mols(m,at,1), mols(m,at,2), mols(m,at,3), connections{k}, types{at});
        end
        milestone=0;
        if m/numMols*100>milestone
            set(handles.errorBox,'String',['Stage 3 of 4: ' num2str(round(m/numMolecules*100)) '% done...'])
            drawnow
            milestone=milestone+1; %#ok<*NASGU>
        end
    end
    
    for b=1:numBonds
        fprintf(fid, '\t\t\t\t\t<Bond ID="%d" UserID="%d" Connects="%d,%d"', ...
            bondIDs(b), b, connects(b,1), connects(b,2));
        if connects(b,3)==2
            fprintf(fid, ' Type="Double"/>\r\n');
        elseif connects(b,3)==3
            fprintf(fid, ' Type="Triple"/>\r\n');
        else
            fprintf(fid, '/>\r\n');
        end
        milestone=0;
        if m/numMols*100>milestone
            set(handles.errorBox,'String',['Stage 4 of 4: ' num2str(round(b/numBonds*100)) '% done...'])
            drawnow
            milestone=milestone+1;
        end
    end
    
    fprintf(fid, '\t\t\t\t</SubUnit>\r\n');
    fprintf(fid, '\t\t\t</Chain>\r\n');
    fprintf(fid, '\t\t</Molecule>\r\n');
    fprintf(fid, '\t</AtomisticTreeRoot>\r\n');
    fprintf(fid, '</XSD>\r\n');
    fclose(fid);
    set(handles.errorBox,'String',[outfile ' successfully exported.'])
else
    set(handles.errorBox,'String','One or more variables necessary for exporting are missing.')
end

function writePDB(handles,outfile)
varsLoaded=1;
try
    if get(handles.exportSlot1Radio,'Value')
        mols=handles.mols1;
    else
        mols=handles.mols2;
    end
catch
    varsLoaded=0;
end

if varsloaded
    numMols=length(mols);
    molLengths=zeros(1,numMols);
    for m=1:numMols
        molLengths(m)=length(mols{m});
    end
    numAtoms=sum(molLengths);
    positions=zeros(numAtoms,3);
    index=1;
    for m=1:numMols
        positions(index:index+molLengths(m)-1,:)=mols{m};
        index=index+molLengths(m);
    end
    
    S=date;
    C=clock;
    C=round(C(4:6));
    
    [~,dayname]=weekday(S);
    
    ofid=fopen(outfile,'w+');
    fprintf(ofid, 'REMARK   Materials Studio PDB file\r\n');
    fprintf(ofid, 'REMARK   Created:   %s %s %s %02.0f:%02.0f:%02.0f Eastern Standard Time %s\r\n', dayname, S(4:6), S(1:2), C(1), C(2), C(3), S(end-3:end));
    try
        latticeParams=handles.latticeParams;
        supercell=handles.supercell;
        fprintf(ofid, 'CRYST1 % 6.3f % 6.3f % 6.3f % 5.2f % 5.2f% 5.2f P1\r\n', latticeParamss(1)*supercell(1), latticeParams(2)*supercell(2), latticeParams(3)*supercell(3), latticeParameters(4), latticeParameters(5), latticeParameters(6));
    catch
    end
    for at=1:numAtoms
        fprintf(ofid,'ATOM% 7d % 3s  MOL     2    % 8.3f% 8.3f% 8.3f  1.00  0.00         % 3s\r\n', at,rawData{3}{at}, positions(at,1),positions(at,2),positions(at,3),rawData{11}{at});
    end
    fprintf(ofid,'TER\r\n');
    fclose(ofid);
else
    set(handles.errorBox,'String','One or more variables necessary for exporting are missing.')
end

function writeMols(handles,outfile)
varsLoaded=1;
try
    if get(handles.exportSlot1Radio,'Value')
        mols=handles.mols1;
    else
        mols=handles.mols2;
    end
catch
    varsLoaded=0;
end

if varsLoaded==1;
    save(outfile,'mols')
else
    set(handles.errorBox,'String','One or more variables necessary for exporting are missing.')
end

function outfileEdit_Callback(hObject, eventdata, handles)
% hObject    handle to outfileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of outfileEdit as text
%        str2double(get(hObject,'String')) returns contents of outfileEdit as a double


% --- Executes during object creation, after setting all properties.
function outfileEdit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to outfileEdit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in tempMapCheckBox.
function tempMapCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to tempMapCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of tempMapCheckBox
if get(hObject, 'Value')
    set(handles.plotSlot1Radio,'Enable','off')
    set(handles.plotSlot2Radio,'Enable','off')
else
    set(handles.plotSlot1Radio,'Enable','on')
    set(handles.plotSlot2Radio,'Enable','on')
end


% --- Executes on button press in plotButton.
function plotButton_Callback(hObject, eventdata, handles)
% hObject    handle to plotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.errorBox,'String','')
if ~get(handles.tempMapCheckBox,'Value')
    if get(handles.plotSlot1Radio,'Value')
        try
            handles.mols1;
            plotMols(handles,1)
        catch
            set(handles.errorBox,'String','No molecules found in Slot 1.')
        end
    else
        try
            handles.mols2;
            plotMols(handles,2)
        catch
            set(handles.errorBox,'String','No molecules found in Slot 2.')
        end
    end
else
    plotTempMap(handles)
end

function plotTempMap(handles)
mols1=handles.mols1;
mols2=handles.mols2;

numMols=length(mols1);
if numMols~=length(mols2)
    set(handles.errorBox,'String', 'The cell sizes are inconsistent.')
else
    try
        temperature=zeros(length(mols1),1);
        positions=zeros(length(mols1),3);
        for m=1:numMols
            temperature(m)=sum(sum((mols1{m}-mols2{m}).^2)).^.5;
            positions(m,:)=mean(mols1{m},1);
        end
      
        meantemp=mean(temperature);
        stdtemp=std(temperature);
        temperature(temperature<(meantemp-4*stdtemp))=meantemp-4*stdtemp;
        temperature(temperature>(meantemp+4*stdtemp))=meantemp+4*stdtemp;
        figure(handles.plot1)
        scatter3(positions(:,1),positions(:,2),positions(:,3),35,temperature,'filled')
    catch
        set(handles.errorBox,'String','The cell sizes are inconsistent.')
    end
end


% --- Executes on button press in loadFormatButton.
function loadFormatButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadFormatButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.errorBox,'String','')
order=get(handles.orderTable,'Data');
if ~iscell(order)
    tempOrder=order;
    molsPerCell=size(order,1);
    order=cell(molsPerCell,1);
    for k=1:molsPerCell
        order{k}=tempOrder(k,tempOrder(k,:)~=0);
    end
    handles.order=order;
end

bonds=get(handles.bondsTable,'Data');
if ~iscell(bonds)
    if nnz(bonds)==numel(bonds)
        handles.bonds=bonds;
    else
        set(handles.errorBox,'String','The Bonds array contains zeros, couldn''t update.')
    end
end

types=get(handles.typesTable,'Data');
empt=0;
for k=1:numel(types)
    empt=empt+double(isempty(types{k}));
end
if ~empt
    handles.types=types;
elseif ~isempty(types{1})
    set(handles.errorBox,'String','The Types array contains empty spaces, couldn''t update.')
end
set(handles.errorBox,'String','Successfully updated infile format parameters.')
guidata(handles.figure1,handles)




% --- Executes on button press in dislocateButton.
function dislocateButton_Callback(hObject, eventdata, handles)
% hObject    handle to dislocateButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.errorBox,'String','')
location=handles.location';
sense=handles.sense';
burgers=handles.burgers';
orientation=handles.orientation';
latticeParams=handles.latticeParams;
rotafier=handles.rotation;

rigid=get(handles.rigidCheckBox,'Value');

mols=handles.mols1;
numMols=length(mols);

r_gradient=0:.001:1;
r_gradLength=1001;
th_gradient=0:.001:1;
th_gradLength=1001;

a=latticeParams(1); b=latticeParams(2); c=latticeParams(3);
alpha=latticeParams(4); beta=latticeParams(5); gamma=latticeParams(6);
bVect=[0; b; 0];
aVect=[a*sind(gamma); a*cosd(gamma);0];
c2=c*cosd(alpha);
c1=(a*c*cosd(beta)-c*cosd(alpha)*cosd(gamma))/(a*sind(gamma));
c3=sqrt((c*sind(alpha))^2-c1^2);
cVect=[c1;c2;c3];

crystalizer=rotafier\[aVect,bVect,cVect];

rCrossE=cross(orientation,sense);
normE=norm(sense);
centralizer=[cross(sense,rCrossE)/norm(cross(sense,rCrossE)) rCrossE/norm(rCrossE) sense/normE]';

if rigid
    bPositions=zeros(numMols,3);
    for m=1:numMols
        bPositions(m,:)=mean(mols{m});
    end
else
    molLengths=zeros(1,length(mols));
    for m=1:numMols
        molLengths(m)=size(mols{m},1);
    end
    index=1;
    bPositions=zeros(sum(molLengths),3);
    for m=1:numMols
        bPositions(index:index+molLengths(m)-1,:)=mols{m};
        index=index+length(mols{m});
    end
end

numPoints=size(bPositions,1);
mover=diag(location)*ones(3,numPoints);
cPositions=(centralizer*(crystalizer\bPositions'-mover))';
aPositions=zeros(size(cPositions));
centralized_burgers=centralizer*burgers;

for m=1:numPoints
    x0=cPositions(m,2);
    y0=cPositions(m,1);
    if x0~=0
        newBurger=centralized_burgers*x0/abs(2*x0);
        angleRatio=((pi-mod(atan(abs(x0)/y0),pi))/pi);
        gradPosition=angleRatio*(th_gradLength-1)+1;
        lower=floor(gradPosition);
        upper=ceil(gradPosition);
        if lower==upper
            th_multiplier=th_gradient(upper);
        else
            th_multiplier=th_gradient(upper)*(gradPosition-lower)+th_gradient(lower)*(upper-gradPosition);
        end
        newBurger=newBurger*th_multiplier;
    else
        newBurger=burgers*0;
    end
    r0=sqrt(x0^2+y0^2);
    if r0>=normE
        aPositions(m,:)=cPositions(m,:)+newBurger';
    else
        distanceRatio=r0/normE;
        gradPosition=distanceRatio*(r_gradLength-1)+1;
        lower=floor(gradPosition);
        upper=ceil(gradPosition);
        if lower==upper
            r_multiplier=r_gradient(upper);
        else
            r_multiplier=r_gradient(upper)*(gradPosition-lower)+r_gradient(lower)*(upper-gradPosition);
        end
        aPositions(m,:)=cPositions(m,:)+r_multiplier*newBurger';
    end
end

aPositions=(crystalizer*(centralizer\aPositions'+mover))';

dislocatedMols=cell(numMols,1);

if rigid
    translator=aPositions-bPositions;
    for m=1:numMols
        dislocatedMols{m}=mols{m}+ones(size(mols{m}))*diag(translator(m,:));
    end
else
    index=1;
    for m=1:numMols
        dislocatedMols{m}=aPositions(index:index+molLengths(m)-1,:);
        index=index+length(mols{m});
    end
end

handles.mols2=handles.mols1;
handles.mols1=dislocatedMols;
set(handles.errorBox,'String','Dislocation successful!')
guidata(handles.figure1,handles);
if get(handles.plotCheckBox,'Value')
    plotMols(handles,1)
end


% --- Executes during object creation, after setting all properties.
function typesTable_CreateFcn(hObject, eventdata, handles)
% hObject    handle to typesTable (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
set(hObject,'ColumnFormat',{'char','char'})
set(hObject,'Data',cell(4,1));


% --- Executes on button press in rigidCheckBox.
function rigidCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to rigidCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of rigidCheckBox


% --- Executes on button press in plotCheckBox.
function plotCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to plotCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of plotCheckBox


% --- Executes during object creation, after setting all properties.
function plotButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to plotButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in loadRotationButton.
function loadRotationButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadRotationButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.errorBox,'String','')
handles.rotation=get(handles.rotationTable,'Data');
set(handles.errorBox,'String','Successfully updated the custom rotation matrix.')
guidata(handles.figure1,handles);


% --- Executes on button press in loadParamsButton.
function loadParamsButton_Callback(hObject, eventdata, handles)
% hObject    handle to loadParamsButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
set(handles.errorBox,'String','')
handles.location=get(handles.locationTable,'Data');
handles.sense=get(handles.senseTable,'Data');
handles.burgers=get(handles.burgersTable,'Data');
handles.orientation=get(handles.orientationTable,'Data');
set(handles.errorBox,'String','Successfully updated dislocation parameters.')
guidata(handles.figure1,handles)


% --- Executes on button press in covalentCheckBox.
function covalentCheckBox_Callback(hObject, eventdata, handles)
% hObject    handle to covalentCheckBox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of covalentCheckBox
if get(hObject, 'Value')
    set(handles.rigidCheckBox,'Value',0)
    set(handles.rigidCheckBox,'Enable','off')
    set(handles.orderTable,'Enable','off')
    set(handles.errorBox,'String', 'Note: Only single elements are supported currently. Bonds array must specify ALL bonds.')
else
    set(handles.rigidCheckBox,'Enable','on')
    set(handles.orderTable,'Enable','on')
    set(handles.errorBox,'String', '')
end


% --- Executes during object creation, after setting all properties.
function exportButton_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exportButton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
