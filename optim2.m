% TM-SWIFT genetic algorithm
% This code contains the GUI elements.
% Need to clean this for command line use

function varargout = optim2(varargin)
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @optim2_OpeningFcn, ...
                   'gui_OutputFcn',  @optim2_OutputFcn, ...
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


% --- Executes just before optim2 is made visible.
function optim2_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to optim2 (see VARARGIN)

% Choose default command line output for optim2
global balance_score;
global correlation_score;
global distance_score;
global waste_score
global dist_coef;
global wait_coef;
global balance_coef;
global waste_coef;
global iterations
global rjcted;
global CurrentTime;
global StartTime;
global EndTime;
global TeamStatus;

handles.output = hObject;

% Update handles structure
guidata(hObject, handles);
iterations = 200;
balance_score = 0;
distance_score = 0;
correlation_score = 0;
waste_score = 0;
dist_coef = 0.2;
wait_coef = 0;
balance_coef = 0.1;
waste_coef = 0.7;
CurrentTime = 540;
StartTime = 540;
EndTime = 1080;

rjcted = [];
set(handles.edtxtIter,'String',iterations);
set(handles.edtxtCurrTime,'String',CurrentTime);
set(handles.edtxtDcof,'String',dist_coef);
set(handles.edtxtPcof,'String',wait_coef);
set(handles.edtxtBcof,'String',balance_coef);
set(handles.edtxtWcof,'String',waste_coef);
set(handles.stxtDscore,'String',distance_score);
set(handles.stxtPScore,'String',correlation_score);
set(handles.stxtBscore,'String',balance_score);
set(handles.stxtWscore,'String',waste_score);
str = sprintf('0/%d',iterations);


% UIWAIT makes optim2 wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = optim2_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in OrderFileName.
function OrderFileName_Callback(hObject, eventdata, handles)
% hObject    handle to OrderFileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global OrderFileName;
[filename,pathname] = uigetfile;

if filename ~= 0
    OrderFileName = strcat(pathname,filename);
    set(handles.stxtOrderFile,'String',OrderFileName);
else
    set(handles.stxtOrderFile,'String','File Not Selected');
    OrderFileName = 'File Not Selected';
end



function edtxtIter_Callback(hObject, eventdata, handles)
% hObject    handle to edtxtIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtxtIter as text
%        str2double(get(hObject,'String')) returns contents of edtxtIter as a double
global iterations;
iterations = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edtxtIter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtxtIter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnLoad.
function btnLoad_Callback(hObject, eventdata, handles)
% hObject    handle to btnLoad (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

global TeamID
global JobID
global OrderFileName
global TeamFileName
global xy_o
global xy_t
global n_o_loc
global n_t_loc
global tBld
global tSkl
global oBld
global oSkl
global duration
global match_matrix
global total_wts
global oTl
global tTl
global se_no
global TeamStatus
global FxJobFileName
global FxJobs
global CurrentTime


% figure('Name','Task and Team Locations','Numbertitle','off');
scrsz = get(0,'ScreenSize');
figure('Name','Task and Team Locations','Numbertitle','off','Position',[10 scrsz(4)/8 scrsz(3)/1.5 scrsz(4)/1.5]);

if ~strcmp(TeamFileName,'File Not Selected')
    fid = fopen(TeamFileName);
    % 1.TeamID 2.Long 3.Lat 4.Building 5.Skillset 6.Toolbox
    % 7.nxt_avail_time
    B = textscan(fid, '%f %f %f %s %s %s %f', 'delimiter', '\t');
    fclose(fid);
    TeamID = B{1};
    xy_t = [B{2} B{3}];
    tBld = B{4};
    tBld = strrep(tBld,'"','');
    tBld = char(tBld);
    tSkl = B{5};
    tSkl = strrep(tSkl,'"','');
    tSkl = char(tSkl);
    tTl = B{6};
    tTl = strrep(tTl,'"','');
    tTl = char(tTl);
    

    n_t_loc = size(tSkl,1);
    TeamStatus = zeros(n_t_loc,4);
    TeamStatus(:,1) = B{1};
    TeamStatus(:,2) = xy_t(:,1);
    TeamStatus(:,3) = xy_t(:,2);
    TeamStatus(:,4) = B{7};
    plot(xy_t(:,1),xy_t(:,2),'rx');
end

hold on
if ~strcmp(OrderFileName,'File Not Selected')
    fid = fopen(OrderFileName);
    % 1.JobID 2.Long 3.Lat 4.Building 5.Skillset 6.Toolbox 7.Duration
    % 8.Weightage
    B = textscan(fid, '%s %f %f %s %s %s %f %f', 'delimiter', '\t');
    fclose(fid);
    n_o_loc = size(B{1},1);

    
    n_o_loc = size(B{1},1);
    weights = B{8};
    
    [total_wts sort_idx] = sort(weights,'descend');
    JobID = B{1}(sort_idx);
    se_no = B{1}(sort_idx);
    xy_o = [B{2}(sort_idx) B{3}(sort_idx)];
    oBld = B{4}(sort_idx);
    oBld = char(oBld);
    oSkl = B{5}(sort_idx);
    oSkl = char(oSkl);
    oTl = B{6}(sort_idx);
    oTl = char(oTl);
    duration = B{7}(sort_idx);
    plot(xy_o(:,1),xy_o(:,2),'bs');
end

if ~strcmp(FxJobFileName,'File Not Selected')
    fid = fopen(FxJobFileName);
    % 1.TeamID 2.Type 3.Long 4.Lat 5.StartTime 6.Duration 7.JobID 8.Building
    % 9.Skillset 10.Toolbox 11.Weightage
    C = textscan(fid, '%f %f %f %f %f %f %s %s %s %s %f', 'delimiter', '\t');
    fclose(fid);
    FxJobs = C;
end

 
match_matrix = cell(1,n_o_loc);

for i=1:n_o_loc
    a = StringFind(oBld(i,:),tBld);
    b = StringFind(oSkl(i,:),tSkl);
    c = StringFind(oTl(i,:),tTl);
    match_matrix{i} = findmatch(a,b,c);
end

dcm_obj = datacursormode(gcf);
set(dcm_obj,'UpdateFcn',{@dtip_function,tBld,tSkl,oBld,oSkl,xy_t})

hold off



function [idx] = StringFind(srcString,strArray)
idx = [];
l = size(strArray,1);
for i=1:l
    k = strfind(strArray(i,:),srcString);
    if size(k) ~= 0
        idx = [idx i];
    end
end




% --- Executes on button press in btnRun.
function btnRun_Callback(hObject, eventdata, handles)
% hObject    handle to btnRun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
genetic_algo(handles)

function genetic_algo(handles)
global groups
global n_t_loc
global n_o_loc
global xy_o
global xy_t
global c
global match_matrix
global duration
global oBld
global oSkl
global tBld
global tSkl
global oTl
global tTl
global balance_score;
global correlation_score;
global distance_score;
global waste_score;
global iterations
global se_no
global groups_idx
global CurrentTime
global TeamStatus
global FxJobs
global JobID

groups{n_t_loc+1} = [];
groups_idx = cell(1,n_t_loc);
temp_timeline = ones(n_t_loc,1) * inf;
greedySol = zeros(n_t_loc,1);
for i=1:n_o_loc
    tms = match_matrix{i};
    tms_size = size(tms,2);
    temp_timeline1 = inf(n_t_loc,1);
    temp_timeline1(tms) = temp_timeline(tms);
    for j=1:tms_size
        z = size(groups_idx{tms(j)},2);
        if z==0
            dst = sqrt((xy_o(i,1) - xy_t(tms(j),1)).^2+(xy_o(i,2) - xy_t(tms(j),2)).^2);
        else
            tmp_loc = groups_idx{tms(j)}(z);
            dst = sqrt((xy_o(i,1) - xy_o(tmp_loc,1)).^2+(xy_o(i,2) - xy_o(tmp_loc,2)).^2);
        end
        
        if temp_timeline1(tms(j)) ~= inf
            temp_timeline1(tms(j)) =  temp_timeline1(tms(j)) + dst;
        else
            temp_timeline1(tms(j)) = dst;
        end
    end
    [x y] = min(temp_timeline1);
    if temp_timeline(y) ~= inf
        temp_timeline(y) = temp_timeline(y) + x + duration(i);
    else
        temp_timeline(y) = x + duration(i);
    end
    groups_idx{y} = [groups_idx{y} i];
    greedySol(i) = y;
end


% GA code here
PopSize = 30;
best = [];
best_val = inf;
pop = create_population(match_matrix,PopSize);
ftns = zeros(PopSize,1);
for i=1:PopSize
   ftns(i) = get_fitness(pop{i});
end

h = waitbar(0,'Please wait...');

for j=1:iterations
    pop = select_fittest(pop);
    for i=1:PopSize
       ftns(i) = get_fitness(pop{i});
    end

    min(ftns)
    pop = crossover(pop,0.25);
    pop = select_fittest(pop);
    for i=1:PopSize
       ftns(i) = get_fitness(pop{i});
    end

    [p,q] = min(ftns);
    best = pop{q};
    best_val = [best_val p];

    get_fitness(best);
    set(handles.stxtDscore,'String',distance_score);
    set(handles.stxtPScore,'String',correlation_score);
    set(handles.stxtBscore,'String',balance_score);
    set(handles.stxtWscore,'String',waste_score);
    set(handles.stxtFtns,'String',p);

    pop = mutate(pop,0.2,match_matrix);
    pop = select_fittest(pop);
    for i=1:PopSize
       ftns(i) = get_fitness(pop{i});
    end
    [p,q] = min(ftns);
    best = pop{q};
    best_val = [best_val p];
    
    get_fitness(best);
    set(handles.stxtDscore,'String',distance_score);
    set(handles.stxtPScore,'String',correlation_score);
    set(handles.stxtBscore,'String',balance_score);
    set(handles.stxtWscore,'String',waste_score);
    set(handles.stxtFtns,'String',p);
    
    waitbar(j/iterations,h,sprintf('%d/%d Complete',j,iterations));
end 

close(h); 

grp11 = assign_jobs(best);
grp11d = get_distance(grp11);
[schedule_matrix,schedule_matrix_idx] = getScheduleMatrix(grp11,grp11d);



%schedule_matrix{n_t_loc+1} = [];
%for i=1:n_t_loc
%    s = size(groups_idx{i},1);
%    for j=1:s
%       schedule_matrix{i} = [schedule_matrix{i} groups_dist{i}(j)]; 
%       schedule_matrix{i} = [schedule_matrix{i} duration(groups_idx{i}(j))];
%    end
%end
 
scrsz = get(0,'ScreenSize');
figure('Name','Jobs Distribution','Numbertitle','off','Position',[10 scrsz(4)/8 scrsz(3)*0.8 scrsz(4)*0.8]);



for i=1:n_t_loc
    s = size(schedule_matrix{i},2);
    sum1 = CurrentTime;
    if TeamStatus(i,4) > CurrentTime
        sum1 = TeamStatus(i,4);
    end
    
    for j = 1:s
        if schedule_matrix_idx{i}(j) == 1
            line([sum1 sum1+duration(schedule_matrix{i}(j))],[i i],'Marker','none','LineStyle','-','Color','blue','LineWidth',5);
            str = sprintf('%s',char(JobID(schedule_matrix{i}(j))));
            text(sum1, i+0.2 , str, 'Color', 'r','FontSize',8,'FontWeight','light');
            sum1 = sum1 + duration(schedule_matrix{i}(j));

        elseif schedule_matrix_idx{i}(j) == 2
            line([sum1 sum1+schedule_matrix{i}(j)],[i i],'Marker','none','LineStyle','-','Color','red','LineWidth',1);
            sum1 = sum1 + schedule_matrix{i}(j);
        elseif schedule_matrix_idx{i}(j) == 3
            line([sum1 sum1+schedule_matrix{i}(j)],[i i],'Marker','none','LineStyle','-','Color','magenta','LineWidth',2);
            sum1 = sum1 + schedule_matrix{i}(j);
        elseif schedule_matrix_idx{i}(j) == 4
            line([sum1 sum1+FxJobs{6}(schedule_matrix{i}(j))],[i i],'Marker','none','LineStyle','-','Color','yellow','LineWidth',5);
            str = sprintf('%s',char(FxJobs{7}(schedule_matrix{i}(j))));
            text(sum1, i+0.2 , str, 'Color', 'r','FontSize',8,'FontWeight','light');
            sum1 = sum1 + FxJobs{6}(schedule_matrix{i}(j));
        elseif schedule_matrix_idx{i}(j) == 5
            line([sum1 sum1+FxJobs{6}(schedule_matrix{i}(j))],[i i],'Marker','none','LineStyle','-','Color','green','LineWidth',5);
            str = sprintf('%s',char(FxJobs{7}(schedule_matrix{i}(j))));
            text(sum1, i+0.2 , str, 'Color', 'r','FontSize',8,'FontWeight','light');
            sum1 = sum1 + FxJobs{6}(schedule_matrix{i}(j));
        end
    end
    hold on;
end


figure('Name','Convergence','Numbertitle','off');
plot(best_val);


figure('Name','Clustering','Numbertitle','off','Position',[10 scrsz(4)/8 scrsz(3)*0.8 scrsz(4)*0.8]);
ColorSet = varycolor(n_t_loc);
hold on
for i=1:n_t_loc
    s = size(groups_idx{i},1);
    for j=1:s
        line([xy_t(i,1) xy_o(groups_idx{i}(j),1)],[xy_t(i,2) xy_o(groups_idx{i}(j),2)],'Color',ColorSet(i,:));
    end
end

for i=1:n_t_loc
    plot(xy_t(:,1),xy_t(:,2),'rs','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
end

for i=1:n_t_loc
    plot(xy_o(groups_idx{i},1),xy_o(groups_idx{i},2),'*','Color',ColorSet(i,:));
end
hold off

figure('Name','Travelling Path','Numbertitle','off','Position',[10 scrsz(4)/8 scrsz(3)*0.8 scrsz(4)*0.8]);
hold on;

for i=1:n_t_loc
    plot([xy_t(i,1); xy_o(groups_idx{i},1)],[xy_t(i,2); xy_o(groups_idx{i},2)],'-','Color',ColorSet(i,:));
end

for i=1:n_t_loc
    plot(xy_t(:,1),xy_t(:,2),'rs','LineWidth',1,'MarkerEdgeColor','k','MarkerFaceColor','r','MarkerSize',8);
end

for i=1:n_t_loc
    plot(xy_o(groups_idx{i},1),xy_o(groups_idx{i},2),'*','Color',ColorSet(i,:));
end
hold off;

display('=========================== RESULTS =======================================');





%figure('Name','TSP_GA | Results','Numbertitle','off');
%plot([3 2],'-');


% --------------------------------------------------------------------
function uitoggletool3_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to uitoggletool3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in btnTeamFile.
function btnTeamFile_Callback(hObject, eventdata, handles)
% hObject    handle to btnTeamFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global TeamFileName;
[filename,pathname] = uigetfile;

if filename ~= 0
    TeamFileName = strcat(pathname,filename);
    set(handles.stxtTeamFile,'String',TeamFileName);
else
    set(handles.stxtTeamFile,'String','File Not Selected');
    TeamFileName = 'File Not Selected';
end


function output_txt = dtip_function(obj,event_obj,tBld,tSkl,oBld,oSkl,xy_t)
% Display the position of the data cursor
% obj          Currently not used (empty)
% event_obj    Handle to event object
% output_txt   Data cursor text string (string or cell array of strings).

pos = get(event_obj,'Position');
dindex = get(event_obj,'DataIndex');


if dindex > size(xy_t,1)
    output_txt = {['X: ',num2str(pos(1),4)],...
        ['Y: ',num2str(pos(2),4)],['I:',num2str(dindex)],...
        ['J-BLD:',oBld(dindex,:)],...
        ['J-SKL:',oSkl(dindex,:)]};
else
if ((pos(1) == xy_t(dindex,1)) && (pos(2) == xy_t(dindex,2)))
    output_txt = {['X: ',num2str(pos(1),4)],...
        ['Y: ',num2str(pos(2),4)],['I:',num2str(dindex)],...
        ['T-BLD:',tBld(dindex,:)],...
        ['T-SKL:',tSkl(dindex,:)]};
else
    output_txt = {['X: ',num2str(pos(1),4)],...
        ['Y: ',num2str(pos(2),4)],['I:',num2str(dindex)],...
        ['J-BLD:',oBld(dindex,:)],...
        ['J-SKL:',oSkl(dindex,:)]};
end
end


% If there is a Z-coordinate in the position, display it as well
if length(pos) > 2
    output_txt{end+1} = ['Z: ',num2str(pos(3),4)];
end


function Y=findmatch(A,B,C)
Y = [];
for i=1:numel(A)
    for j=1:numel(B)
        for k=1:numel(C)
            if (A(i) == B(j)) && (A(i) == C(k))
                Y = [Y A(i)];
                break;
            end
        end
    end
end

function pop = create_population(match_matrix,PopSize)
ChrLength = size(match_matrix,2);
pop = cell(PopSize,1);
for j=1:PopSize
    temp_chr = zeros(ChrLength,1);
    for i=1:ChrLength
        if size(match_matrix{i},1) ~= 0
            temp_chr(i) = match_matrix{i}(randi(numel(match_matrix{i})));
        end
    end
    pop{j} = temp_chr;
end

function grp = assign_jobs(jobs)
global n_t_loc
grp{n_t_loc} = [];
for i=1:n_t_loc
   k = find(jobs == i);
   grp{i} = k;
end

function ftns = get_fitness(individual)
global dist_coef;
global wait_coef;
global balance_coef;
global waste_coef;
global n_t_loc;
global n_o_loc;
global duration;
global balance_score;
global correlation_score;
global distance_score;
global waste_score;
global total_wts;
global match_matrix

grp = assign_jobs(individual);
dMat = get_distance(grp);
% WaitMat = get_wait_matrix(grp,dMat);


[schMat,schMat_idx] = getScheduleMatrix(grp,dMat);
WstTeam= getWasteTime(schMat,schMat_idx);
waste_score = sum(WstTeam);
total_distance = get_distance_sch(schMat,schMat_idx);
balance_score = std(WstTeam);
bal = balance_score * 10;
%[v,idx] = sort(WaitMat);
%a = [1:n_o_loc];
%cVal = corr(a',idx);
%correlation_score = cVal;
%cVal = (1/(1+cVal))*10000;
distance_score = total_distance;
% ftns = dist_coef*distance_score+wait_coef*cVal+balance_coef*bal+ltime_coef*ltime_score;
ftns = dist_coef*distance_score+balance_coef*bal+waste_coef*waste_score;

function [wTeam] = getWasteTime(schMat,schMat_idx)
global n_t_loc
wTeam(n_t_loc,1) = 0;
for i=1:n_t_loc
    wTeam(i,1) =  sum(schMat{i}(schMat_idx{i} == 3));
end



function tDist = get_distance_sch(schMat,schMat_idx)
global n_t_loc
tDist = 0;
for i=1:n_t_loc
    tDist = tDist + sum(schMat{i}(schMat_idx{i} == 2));
end



function DistMat = get_distance(groups_idx)
global xy_o;
global xy_t;
s = size(groups_idx,2);
DistMat{s} = [];
for i=1:s
    k = groups_idx{i};
    l = size(k,1);
    for j=1:l
       if j==1
           d = sqrt((xy_t(i,1)-xy_o(k(j),1)).^2+(xy_t(i,2)-xy_o(k(j),2)).^2);
       else
           d = sqrt((xy_o(k(j),1)-xy_o(k(j-1),1)).^2+(xy_o(k(j),2)-xy_o(k(j-1),2)).^2);
       end
       DistMat{i}(j) = d;
    end
end



function new_pop = select_fittest(pop)
PopSize = size(pop,1);
ftns = zeros(PopSize,1);
sel = zeros(PopSize,1);
for i=1:PopSize
    ftns(i) = get_fitness(pop{i}); 
end

for i=1:PopSize
    ftns(i) = 1/(1+ftns(i)); 
end

T = sum(ftns);
for i=1:PopSize
    ftns(i) = ftns(i)/T; 
end

cumu = zeros(PopSize,1);
for i=1:PopSize
    if i==1
        cumu(i) = ftns(i);
    else
        cumu(i) = cumu(i-1) + ftns(i);
    end
end

for i=1:PopSize
     r = rand();
     s = find(cumu>r);
     sel(i) = s(1);
end

new_pop = cell(PopSize,1);
for i=1:PopSize
     new_pop{i} = pop{sel(i)};
end

for i=1:PopSize
    ftns(i) = get_fitness(pop{i}); 
end
[p,q] = min(ftns);
for i=1:PopSize
    ftns(i) = get_fitness(new_pop{i}); 
end
[r,s] = max(ftns);
new_pop{s} = pop{q};


function NewPop = crossover(pop,c_rate)
PopSize = size(pop,1);
chr_length = size(pop{1},1);
ftns = inf(PopSize,1);
for i=1:PopSize
    ftns(i) = get_fitness(pop{i}); 
end
[p,q] = min(ftns);
sl = [];
for i=1:PopSize
    r = rand();
    if r<c_rate
        sl = [sl i];
    end
end

% cr_point = randi(PopSize);
% cr_point = (PopSize*0.2) + ceil(cr_point.*(PopSize-2*ceil(PopSize*0.2))./PopSize);
cr_point = randi(chr_length);
cr_point = (chr_length*0.2) + ceil(cr_point.*(chr_length-2*ceil(chr_length*0.2))./chr_length);

s = size(sl);
for i=1:s
    if s==1
        break;
    end
    if i==s
        pop{sl(1)}(1:cr_point) = pop{sl(s)}(1:cr_point);
    else
        p = pop{sl(i)};
        pop{sl(i)}(1:cr_point) = pop{sl(i+1)}(1:cr_point);
    end
end
NewPop = pop;

for i=1:PopSize
    ftns(i) = get_fitness(NewPop{i}); 
end
[r,s] = max(ftns);
NewPop{s} = pop{q};

function NewPop = mutate(pop,mu_rate,match_matrix)
ChrLength = size(pop{1},1);
PopSize = size(pop,1);
ftns = inf(PopSize,1);
for i=1:PopSize
    ftns(i) = get_fitness(pop{i}); 
end
[p,q] = min(ftns);
mu_points = randi(ChrLength,3,1);
for i=1:PopSize
    r = rand();
    if r<mu_rate
        s = size(mu_points);
        for j=1:s
            if match_matrix{mu_points(j)} ~= 0
                pop{i}(mu_points(j)) = match_matrix{mu_points(j)}(randi(numel(match_matrix{mu_points(j)})));
            end
        end
    end
end
NewPop = pop;
for i=1:PopSize
    ftns(i) = get_fitness(NewPop{i}); 
end
[r,s] = max(ftns);
NewPop{s} = pop{q};

function WaitMat = get_wait_matrix(grp_jobs,dist_mat)
global n_t_loc;
global n_o_loc;
global duration
WaitMat = zeros(n_o_loc,1);
for i=1:n_t_loc
    WaitVal = 0;
    s = size(grp_jobs{i},1);
    for j=1:s
        WaitMat(grp_jobs{i}(j)) = WaitVal + dist_mat{i}(j);
        WaitVal = WaitVal + dist_mat{i}(j) + duration(grp_jobs{i}(j));
    end
end


function JobsInTimeFrame = JobsForTimeFrame()


function edtxtDcof_Callback(hObject, eventdata, handles)
% hObject    handle to edtxtDcof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtxtDcof as text
%        str2double(get(hObject,'String')) returns contents of edtxtDcof as a double
global dist_coef;
dist_coef = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function edtxtDcof_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtxtDcof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtxtPcof_Callback(hObject, eventdata, handles)
% hObject    handle to edtxtPcof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtxtPcof as text
%        str2double(get(hObject,'String')) returns contents of edtxtPcof as a double
global wait_coef;
wait_coef = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edtxtPcof_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtxtPcof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edtxtBcof_Callback(hObject, eventdata, handles)
% hObject    handle to edtxtBcof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtxtBcof as text
%        str2double(get(hObject,'String')) returns contents of edtxtBcof as a double
global balance_coef;
balance_coef = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edtxtBcof_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtxtBcof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnRerun.
function btnRerun_Callback(hObject, eventdata, handles)
% hObject    handle to btnRerun (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global se_no
global rjcted
global groups_idx
global n_t_loc
for i=1:n_t_loc
    rjcted = [rjcted se_no(groups_idx{i}(1))];
end




function [schedule_matrix,schedule_matrix_idx] = getScheduleMatrix(grp,dMat)
global duration
global n_t_loc
global FxJobs
global CurrentTime
global TeamStatus
global EndTime
global JobID
global TeamID
global xy_o



schedule_matrix{n_t_loc} = [];
schedule_matrix_idx{n_t_loc} = [];
i = 0;
for t=TeamID'
    temp = find(t == FxJobs{1});
    i = (t == TeamStatus(:,1));
    if size(temp,1) == 0
        disp('getScheduleMatrix : No matching fixed job for team ID');
        continue;
    end

    sTime = CurrentTime;
    
    if TeamStatus(i,4) > CurrentTime
        sTime = TeamStatus(i,4);
    end
    
    start_point = 0;
    start_time = sTime;
    sLocation = [TeamStatus(i,2) TeamStatus(i,3)];
    w = 0;
    for j=1:numel(temp)
         if j ~= 1
             s = numel(schedule_matrix_idx{i});
             for m=(s:-1:1)
                 if ((schedule_matrix_idx{i}(m) == 1) || (schedule_matrix_idx{i}(m) == 4))
                     if schedule_matrix_idx{i}(m) == 1
                         sLocation = [xy_o(schedule_matrix{i}(m),1) xy_o(schedule_matrix{i}(m),2)];
                     else
                         sLocation = [FxJobs{3}(schedule_matrix{i}(m)) FxJobs{4}(schedule_matrix{i}(m))];
                     end
                     break;
                 end
             end
         end
         % x - current number of jobsi in team
         % sLocation - start location for current empty time blob
         % sTime1 - start time for current time blob (sTime = start time for team either next_avail_time or current time)
         % times - the time column of fixed jobs
         % grp - The jobs allocated to teams by randomness
         % dMat - the distance between jobs allocated by randomness
         % x - Location of coming fixed job
         % y - Location of coming fixed job
         [k,l,p,loc,w] = FillJobs(start_point,sLocation,start_time,FxJobs{5}(temp(j)),grp{i},dMat{i},FxJobs{3}(temp(j)),FxJobs{4}(temp(j)));
         start_point = start_point + p;
         schedule_matrix{i} = [schedule_matrix{i} k];
         schedule_matrix_idx{i} = [schedule_matrix_idx{i} l];

         schedule_matrix{i} = [schedule_matrix{i} FxJobs{5}(temp(j))-sum_sch(schedule_matrix{i},schedule_matrix_idx{i})-sTime-w w temp(j)];
         if FxJobs{2}(temp(j)) == 0
             schedule_matrix_idx{i} = [schedule_matrix_idx{i} 2 3 5];
         else
             schedule_matrix_idx{i} = [schedule_matrix_idx{i} 2 3 4];
         end
         
         start_time = sum_sch(schedule_matrix{i},schedule_matrix_idx{i}) + sTime;
         sLocation = loc;
    end
    % sending endtime as 1500 means no endtime as maximum endtime possible
    % is 24X60=1440
    [k,l,p,loc,w] = FillJobs(start_point,sLocation,start_time,1500,grp{i},dMat{i},0,0);
    schedule_matrix{i} = [schedule_matrix{i} k];
    schedule_matrix_idx{i} = [schedule_matrix_idx{i} l];
end


function SchSum = sum_sch(x,y)
global duration
global FxJobs
SchSum = 0;
SchSum = SchSum + sum(x(y == 2));
SchSum = SchSum + sum(x(y == 3));
SchSum = SchSum + sum(duration(x(y == 1)));
SchSum = SchSum + sum(FxJobs{6}(x((y == 4) | (y == 5))));


% schedule_matrix_idx 
% 1 - Normal Job
% 2 - Distance
% 3 - Unused time
% 4 - Fixed Job with location
% 5 - Fixed Job without location

function [jobset1,jobset_idx,p,loc,wTime] = FillJobs(len_sch,sLocation,sTime,eTime,grp1,dMat1,xcor,ycor)
global duration
global xy_o
global xy_t
jobset = [];
jobset1 = [];
jobset_idx = [];
loc = sLocation;
p = 0;
wTime = 0;

s = size(grp1,1);
dist = 0;


if size(grp1,1) == 0
    if (xcor ~= 0) && (ycor ~= 0)
        dist = sqrt((sLocation(1) - xcor).^2+(sLocation(2) - ycor).^2);
        wTime = wTime + (eTime - sTime - sum(jobset)-dist);
    else
        wTime = wTime + (eTime - sTime - sum(jobset));
    end
   return 
end

if eTime > 1440
    for f=len_sch+1:s
        jobset1 = [jobset1 dMat1(f) grp1(f)];
        jobset_idx = [jobset_idx 2 1];
    end
    return
end


if len_sch ~= size(grp1,1)
    dMat1(len_sch + 1) = sqrt((xy_o(grp1(len_sch+1),1)-sLocation(1)).^2 + (xy_o(grp1(len_sch+1),2)-sLocation(2)).^2);
else
    if (xcor ~= 0) && (ycor ~= 0)
        dist = sqrt((sLocation(1) - xcor).^2+(sLocation(2) - ycor).^2);
        wTime = wTime + (eTime - sTime - sum(jobset)-dist);
    else
        wTime = wTime + (eTime - sTime - sum(jobset));
    end
    return
end





% Loop through the Jobs and Distance vectors till their sum equals or less than
% available time slot
for f=len_sch+1:s
   if dist<(eTime-sTime)
      jobset = [jobset dMat1(f) duration(grp1(f))];
      jobset1 = [jobset1 dMat1(f) grp1(f)];
      jobset_idx = [jobset_idx 2 1];
      dist = dist + dMat1(f) + duration(grp1(f));
      if f==s
          if (xcor ~= 0) && (ycor ~= 0)
                if numel(jobset) == 0
                    temp_dist = sqrt((xcor - xy_t(f,1)).^2+(ycor - xy_t(f,2)).^2);
                else
                    temp_dist = sqrt((xcor - (xy_o(grp1(numel(jobset)/2),1))).^2 + (ycor - (xy_o(grp1(numel(jobset)/2),2))).^2);
                end
                if temp_dist > (eTime - sTime - sum(jobset))
                    if numel(jobset) ~= 0
                        jobset = jobset(1:numel(jobset)-2);
                        jobset1 = jobset1(1:numel(jobset1)-2);
                        jobset_idx = jobset_idx(1:numel(jobset_idx)-2);
                    end
                end
                wTime = wTime + (eTime - sTime - sum(jobset)-temp_dist);
                loc = [xcor ycor];                
          else
                if dist>(eTime-sTime)
                    jobset = jobset(1:numel(jobset)-2);
                    jobset1 = jobset1(1:numel(jobset1)-2);
                    jobset_idx = jobset_idx(1:numel(jobset_idx)-2);
                end

                wTime = wTime + (eTime - sTime - sum(jobset));
                if numel(jobset)~= 0
                    loc = [xy_o(grp1(numel(jobset)/2),1) xy_o(grp1(numel(jobset)/2),2)];
                else
                    loc = sLocation;
                end
           end
           p = numel(jobset)/2;
           return
      end
   else
       % Keep on removing 2 elements each time for Jobset till it meets the
       % criteria that, all jobs and distance fit into the time slot

       while numel(jobset) >= 2
            if (xcor ~= 0) && (ycor ~= 0)
                jobset = jobset(1:numel(jobset)-2);
                jobset1 = jobset1(1:numel(jobset1)-2);
                jobset_idx = jobset_idx(1:numel(jobset_idx)-2);
                if numel(jobset) == 0
                    temp_dist = sqrt((xcor - xy_t(f,1)).^2+(ycor - xy_t(f,2)).^2);
                else
                    temp_dist = sqrt((xcor - (xy_o(grp1(numel(jobset)/2),1))).^2 + (ycor - (xy_o(grp1(numel(jobset)/2),2))).^2);
                end
                if temp_dist > (eTime - sTime - sum(jobset))
                    if numel(jobset) ~= 0
                        continue;
                    end
                else
                    wTime = wTime + (eTime - sTime - sum(jobset)-temp_dist);
                    loc = [xcor ycor];
                end
            else
                jobset = jobset(1:numel(jobset)-2);
                jobset1 = jobset1(1:numel(jobset1)-2);
                jobset_idx = jobset_idx(1:numel(jobset_idx)-2);
                wTime = wTime + (eTime - sTime - sum(jobset));
                if numel(jobset)~= 0
                    loc = [xy_o(grp1(numel(jobset)/2),1) xy_o(grp1(numel(jobset)/2),1)];
                else
                    loc = sLocation;
                end
            end
            p = numel(jobset)/2;
            return
       end
   end
end


function edtxtWcof_Callback(hObject, eventdata, handles)
% hObject    handle to edtxtWcof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtxtWcof as text
%        str2double(get(hObject,'String')) returns contents of edtxtWcof as a double
global waste_coef;
waste_coef = str2double(get(hObject,'String'));

% --- Executes during object creation, after setting all properties.
function edtxtWcof_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtxtWcof (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in btnFxJob.
function btnFxJob_Callback(hObject, eventdata, handles)
% hObject    handle to btnFxJob (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
global FxJobFileName;
[filename,pathname] = uigetfile;

if filename ~= 0
    FxJobFileName = strcat(pathname,filename);
    set(handles.stxtFxJobFile,'String',FxJobFileName);
else
    set(handles.stxtFxJobFile,'String','File Not Selected');
    FxJobFileName = 'File Not Selected';
end



function edtxtCurrTime_Callback(hObject, eventdata, handles)
% hObject    handle to edtxtCurrTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edtxtCurrTime as text
%        str2double(get(hObject,'String')) returns contents of edtxtCurrTime as a double
global CurrentTime;
CurrentTime = str2double(get(hObject,'String'));


% --- Executes during object creation, after setting all properties.
function edtxtCurrTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edtxtCurrTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% X - Latitude, Y-Longitude
function dist = GPSDistance(x1,y1,x2,y2)
DegToRad = 0.0174532925;
x1=x1*DegToRad;
y1=y1*DegToRad;
x2=x2*DegToRad;
y1=y1*DegToRad;
% Earths mean radius 6371km
R = 6371;
DeltaX = x1-x2;
DeltaY = y1-y2;
A = (sin(DeltaX/2)).^2 + cos(x1)*cos(x2)*(sin(DeltaY)).^2;
C = 2*atan2(sqrt(A),sqrt(1-A));
dist = R*C;
