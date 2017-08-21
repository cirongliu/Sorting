%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Recognition of overlapping waveform
%select the segment which is overlapped period and then
%using the templates to find the matching period
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [New_segment] = Overlapping_recognition_optimized(segment, templates)
function [Spike_list] = Overlapping_recognition(Time_course, segment, templates, Spike_list, Sampling_Freq)

% Check of segment connectivity 
[Spike_list] = overlapping_match([Time_course;segment],Spike_list);%Checking the connectivity

% Dynamic mapping method for overlapping
A = segment;
New_segment = segment;
for Num = 1:length(templates)%For each template
    B = templates{Num}(1,:);%template
    adjustment_window_size = length(B);
    [~, d, ~, path] = dtw(A, B, adjustment_window_size);
    if isempty(path) ~= 1
        [~,New_path] = show_mapping_distance(d, path, A, B, Spike_list);
    else
        New_path = [];
    end
    if size( New_path,2 ) == 4
        temp = New_path(find(New_path(:,4)==1),1);
        Time_begin = 0.2*Sampling_Freq*2;%points
        Time_end = 0.5*Sampling_Freq*2;%points
        
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         figure,subplot(1,3,1),show_distance_matrix(d, path);
%         subplot(1,3,2),
%         plot(A,'-k');
%         hold on
%         plot((temp-Time_begin):(temp+Time_end),B,'-r');
%         hold off
%         pt_x = get(gca,'XLim');
%         pt_y = get(gca,'YLim');
%         title('Matching part');
%         box off
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        if temp<=Time_begin
            New_segment(1+Num,1:(temp+Time_end)) = Num;
            mark = 1;
        else if (temp+Time_end) > size(segment,2)
                New_segment(1+Num,(temp-Time_begin):end) = Num;
                mark = 2;
            else
                New_segment(1+Num,(temp-Time_begin):(temp+Time_end)) = Num;
                mark = 3;
            end
        end
        %Subtract the template and check the noise
        result = segment;
        if mark == 1
            result(1:(temp+Time_end)) = result(1,1:(temp+Time_end))-templates{Num}(1,(Time_begin-temp+2):end);
            temp_length = size( templates{1}, 2);
            if length(1:(temp+Time_end))<temp_length
                Column(Num,1:temp_length) = 1000;
                Column(Num,1:length(1:(temp+Time_end))) = 1:(temp+Time_end);
                Waveform_raw(Num,1:temp_length) = NaN;
                Waveform_raw(Num,1:length(1:(temp+Time_end))) = segment(1:(temp+Time_end));
                Waveform(Num,1:temp_length) = NaN;
                Waveform(Num,1:length(1:(temp+Time_end))) = result(1:(temp+Time_end));
            else
                Column(Num,1:length(1:(temp+Time_end))) = 1:(temp+Time_end);
                Waveform_raw(Num,1:length(1:(temp+Time_end))) = segment(1:(temp+Time_end));
                Waveform(Num,1:length(1:(temp+Time_end))) = result(1:(temp+Time_end));
            end
        end
        if mark == 2
            result((temp-Time_begin):end) = result(1,(temp-Time_begin):end)-templates{Num}(1,1:length(find(New_segment(1+Num,(temp-Time_begin):end) == Num)));
            temp_length = size( templates{1}, 2);
            if length((temp-Time_begin):size(New_segment,2))<temp_length
                Column(Num,1:temp_length) = 1000;
                Column(Num,1:length((temp-Time_begin):size(New_segment,2))) = (temp-Time_begin):size(New_segment,2);
                Waveform_raw(Num,1:temp_length) = NaN;
                Waveform_raw(Num,1:length((temp-Time_begin):size(New_segment,2))) = segment((temp-Time_begin):end);
                Waveform(Num,1:temp_length) = NaN;
                Waveform(Num,1:length((temp-Time_begin):size(New_segment,2))) = result((temp-Time_begin):end);
            else
                Column(Num,1:length((temp-Time_begin):size(New_segment,2))) = (temp-Time_begin):size(New_segment,2);
                Waveform_raw(Num,1:length((temp-Time_begin):size(New_segment,2))) = segment((temp-Time_begin):end);
                Waveform(Num,1:length((temp-Time_begin):size(New_segment,2))) = result((temp-Time_begin):end);
            end
        end
        if mark == 3
            result((temp-Time_begin):(temp+Time_end)) = result(1,(temp-Time_begin):(temp+Time_end))-templates{Num}(1,:);
            Column(Num,:) = (temp-Time_begin):(temp+Time_end);
            Waveform_raw(Num,:) = segment((temp-Time_begin):(temp+Time_end));
            Waveform(Num,:) = result((temp-Time_begin):(temp+Time_end));
        end
        
%         %%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         subplot(1,3,3),
%         plot(result,'-m');
%         title('After substraction');
%         ax = gca;
%         ax.YLim = pt_y;
%         ax.XLim = pt_x;
%         box off
%         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         %%%%%%%%%%%%%%%%%%%%%%%%%
    else
        Column(Num,1:size( templates{1}, 2)) = NaN;
        Waveform_raw(Num,1:size( templates{1}, 2)) = NaN;
        Waveform(Num,1:size( templates{1}, 2)) = NaN;
    end
end

%% Optimization
[C,~,ic] = unique(Column,'rows');%Search replicated value
template_match_temp = [];
if size(C,1)<size(Column,1)
    for i = 1:size(C,1)
        temp = find(ic == i);
        if length(temp)>1
            template_match_temp = [template_match_temp,i];
        end
    end
end
clear C temp i

if isempty(template_match_temp)~=1
    for i = 1:length(template_match_temp)
        temp = find(ic == template_match_temp(i));
        error = Waveform(temp,:);
        Worse_error = temp( find( nanmean(abs(error),2) ~= min(nanmean(abs(error),2)) ) );
        Column(Worse_error,:)= NaN ;
        Waveform_raw(Worse_error,:) = NaN;
        Waveform(Worse_error,:) = NaN;
    end
    clear i temp error Worse_error
else
    clear ic
end
clear template_match_temp

%Mark the highest peak
for i = 1:size(Column,1)
    if isnan( Column(i,1) ) ~= 1
        temp = find( abs( Waveform_raw(i,:) ) == max( abs(Waveform_raw(i,:)) ) );  
        Waveform_raw(i,temp) = 1000;
    end
end
clear i temp

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i = 1:size(Spike_list,1)
    for j = 1:size(Column,1)
        if isnan( Column(j,1) ) ~= 1
            temp_mark = find(Column(j,:) == Spike_list(i,8));%Try to find reasonable index
            if isempty(temp_mark)~=1
                mark(j,1) = 1;
                mark(j,2) = temp_mark;
            else
                mark(j,1) = 0;
                mark(j,2) = 0;
            end
        else
            mark(j,1) = 0;
            mark(j,2) = 0;
        end
    end
    clear temp_mark j
    temp = find(mark(:,1)==1);
    if isempty(temp) == 1
        Spike_list(i,9) = 0;%New type spike
    else
        if length(temp)>1
            for k = 1:length(temp)
                result_temp(k,:) = Waveform(temp(k),mark(temp(k),2));
            end
            temp_mark = find( abs(result_temp) == min(abs(result_temp)) );
            Spike_list(i,9) = temp(temp_mark);
            clear k temp_mark result_temp
        else
            Spike_list(i,9) = temp;%corresponding type
        end
    end
    clear mark temp
end
clear i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %Remove impossible spikes
% temp = unique(Spike_list(:,8)) ;
% if length(temp) < size(Spike_list,1)
%     for i = 1:length(temp)
%         temp1 = find( Spike_list(:,8) == temp(i));
%         if length(temp1)>1
%             temp2 = find ( Spike_list(temp1,4) ~= max( Spike_list(temp1,4) ) );
%             Spike_list(temp1(temp2),8) = NaN;
%         end
%         clear temp1 temp2
%     end
%     clear i
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% if isempty( find(isnan(Spike_list(:,8))==1) )~=1
%     Spike_list(find(isnan(Spike_list(:,8))==1),:) = [];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Spike_list(:,5:8) = [];

function [min_distance, d, g, path] = dtw(A, B, adjustment_window_size)
% Minimal time normalized dtw distance between speech patterns A and B.

r = adjustment_window_size;

% get length of speech patterns A and B
[~, I] = size(A);
[~, J] = size(B);

% local distance matrix
d = zeros(I, J);
for i = 1:I
    for j = 1:J
        feature_distances = A(:,i)-B(:,j);
        distance = sqrt(sum(feature_distances.^2));
        %         d(i,j) = distance;%Normalized Euclidean distance
        d(i,j) = distance/(I+J);%Normalized Euclidean distance
    end
end

% global distance matrix
g = zeros(I+1,J+1);
g(:,:) = inf;
g(1,1) = 2*d(1,1); % initial condition, see (19) in [SakoeChiba1978]
s = J/I; % slope from (0,0) to (I,J)
steps = zeros(I,J); % steps to take in order to reach D(i,j)

for i = 2:I+1;
    for j = 2:J+1;
        if (abs(i-(j/s)) > r)
            % we're outside the adjustment window
            continue;
        end
        
        % local distance matrix is smaller than g, translate coordinates
        i_l = i-1;
        j_l = j-1;
        
        % calculate global distances
        % (see DP-equation (20) from [SakoeChiba1978] for reference)
        [distance, step] =  min([g(i,   j-1) +   d(i_l, j_l)...
            g(i-1, j-1) + 2*d(i_l, j_l)...
            g(i-1, j)   +   d(i_l, j_l)]);
        g(i,j) = distance;
        steps(i-1,j-1) = step;
    end
end

% time normalize global distance matrix
N=I+J;
D=g/N;

% remove additional inf padded row and column from global distance matrix
D=D(2:end,2:end);

path=traceback_path(steps);

min_distance = D(end, end);

function path = traceback_path(steps)
% Traceback a path through a distance array by analyzing a steps matrix.

% (Copied almost verbatim from [Ellis2003].)

path=[];
[i, j] = size(steps);

% trace path back from bottom right to top left
while i>1 && j>1
    switch steps(i,j)
        case 1
            % we got here from g(i, j-1)
            j = j-1;
        case 2
            % we got here from g(i-1, j-1)
            i = i-1;
            j = j-1;
        case 3
            % we got here from g(i-1, j)
            i = i-1;
        otherwise
            path = [];
            break
            error('Oh noes!!1 This cannot have happened.');
    end
    path = [[i j]; path]; %#ok
end

function show_distance_matrix(D, path)

imagesc(D);
colormap(1-gray)
hold on;
plot(path(:,2), path(:,1), 'r');
hold off;
colorbar
box off

function [Mapping, path] = show_mapping_distance(D, path, Raw_data, template, Spike_list)
for i = 1:size(path,1)
    path(i,3) = D(path(i,1),path(i,2));
end
clear i
Mapping =[];
time_template = unique(path(:,2));
for i = 1:length( time_template )
    temp = find( path(:,2) == time_template(i) );
    if length(temp)~=1
        A = path(temp,3);
        temp2 = find(A == min(A));
        Mapping = [Mapping;path(temp( temp2(1) ),:)];
    else
        Mapping = [Mapping;path(temp,:)];
    end
end
%the point corresponding to the Peak of Mapping, according to the function
%of Waveform_extraction
% Time_begin = 0.2*40*2;%points
% Time_end = 0.5*40*2;%points
Mapping( find( Mapping(:,2)==17 ), 4) = 1;%Peak points

%% 
% Time_begin = 0.2*40*2;%points
% Time_end = 0.5*40*2;%points
temp = find( path(:,2) == 17 );
if isempty(temp) ~= 1
    temp2 = Raw_data(path(temp,1));
    if template(17)<0
        temp3 = find( temp2 == min(temp2) );
        temp4 = find( Spike_list(:,3)<0 );
    else
        temp3 = find( temp2 == max(temp2) );
        temp4 = find( Spike_list(:,3)>0 );
    end
    %Update peak
    if isempty(temp4)~=1
        temp5 = find( abs( temp(temp3) - Spike_list(temp4,8)) == min(abs(temp(temp3) - Spike_list(temp4,8))) );
        if length(temp5)>1
            temptemp = find( abs( template(17)- Spike_list(temp4,3) )== min( abs(template(17)- Spike_list(temp4,3)) ) );
            temp6 = find( path(:,1) == Spike_list(temp4(temptemp),8));
        else
            temp6 = find( path(:,1) == Spike_list(temp4(temp5),8));
        end
        if isempty(temp6)~=1
            path(temp6,4)=1;
        end
    end
end


function [num_list] = overlapping_match(Waveform_process,num_list)
% close all

% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% figure,
% plot(Waveform_process(1,:),Waveform_process(2,:),'-k');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%?Spike point
Spike_point = [];
for i = 1:size(num_list,1)%each spike
    Spike_point = [Spike_point,find(round(Waveform_process(1,:)*1000)/1000 == round(num_list(i,2)*1000)/1000)];
    
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     hold on
%     plot(Waveform_process(1,Spike_point(i)),Waveform_process(2,Spike_point(i)),'ob','MarkerSize',15);
%     hold off
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end

%% monotonicity check
connection = zeros(1,length(Spike_point));
temp_Spike_point = Spike_point;
num_list(:,8)=0;%Initialization
temp_Spike_point(2,:)=0;%Initialization
for i = 2:length(Spike_point)
    temp_wave = Waveform_process(2,Spike_point(i-1):Spike_point(i));%trend
    y = diff(temp_wave);
    y(1)=[];y(end)=[];
    if all(y>0)|| all(y<0)%Increase or Decrease
        connection(i-1:i) = 1;
        if length(find(connection == 1))==3 
            if length(find(connection == 1))~=size(num_list,1)
                temp = find(connection==1);
                temp_Spike_point(2,temp) = temp_Spike_point(2,temp)+1;
                num_list(temp,8)=num_list(temp,8)+1;
                temp1 = find(num_list(temp,2) == min(num_list(temp,2)));
                temp_Spike_point(2,temp(temp1)) = temp_Spike_point(2,temp(temp1))+nan;
                num_list(temp(temp1),8)=num_list(temp(temp1),8)+nan;
                clear temp temp1
                connection = zeros(1,length(Spike_point));
            else   
                temp = find(num_list(:,3) < 0); 
                if length(temp)>1
                    num_list(:,8)=num_list(:,8)+1;
                    temp_Spike_point(2,:) = temp_Spike_point(2,:)+1;
                    num_list(temp,8)=num_list(temp,8)+nan;
                    temp_Spike_point(2,temp) = temp_Spike_point(2,temp)+nan;
                else
                    num_list(temp,8)=num_list(temp,8)+nan;
                    temp_Spike_point(2,temp) = temp_Spike_point(2,temp)+nan;
                    num_list(temp+1,8)=num_list(temp+1,8)+1;
                    temp_Spike_point(2,temp+1) = temp_Spike_point(2,temp+1)+1;
                end
                connection = zeros(1,length(Spike_point));
            end
        end
    else
        if isempty(find(connection == 1))~=1
            [num_list,temp_Spike_point] = reorgnize(connection,num_list,temp_Spike_point);
            connection = zeros(1,length(Spike_point));
        end
    end
    clear y temp_wave
end
if isempty(find(connection == 1))~=1
    [num_list,temp_Spike_point] = reorgnize(connection,num_list,temp_Spike_point);
end
temp = find(isnan(num_list(:,8)) ~= 1 & num_list(:,8) ~=0 );
if isempty(temp) ~= 1
    num_list(temp,:)=[];
    temp_Spike_point(:,temp)=[];
end
num_list(:,8)=[];
temp_Spike_point(2,:)=[];
num_list = [num_list,temp_Spike_point'];
% Spike_point = temp_Spike_point;
clear i connection temp_Spike_point temp


function [num_list,temp_Spike_point] = reorgnize(connection,num_list,temp_Spike_point)
temp = find(connection == 1);
if isempty(temp)~=1
    temp1 = [0,diff(temp)];
    temp2 = find(temp1~=1);
    List_spike = temp2;
    clear temp1 temp2
    for i = 1:length(List_spike)
        if length(List_spike) == 1 || i == length(List_spike)
            Array_temp = temp(List_spike(i):end);
            %Non important points
            Array_temp1 = find( num_list(Array_temp,4) ~= max(num_list(Array_temp,4)) );
            num_list(temp(Array_temp1),8)=num_list(temp(Array_temp1),8)+1;
            temp_Spike_point(2,temp(Array_temp1))=temp_Spike_point(2,temp(Array_temp1))+1;
            %Import points
            Array_temp1 = find( num_list(Array_temp,4) == max(num_list(Array_temp,4)) );
            num_list(temp(Array_temp1),8)=num_list(temp(Array_temp1),8)+nan;
            temp_Spike_point(2,temp(Array_temp1))=temp_Spike_point(2,temp(Array_temp1))+nan;
        else
            Array_temp = temp(List_spike(i):List_spike(i+1));
            %Non important
            Array_temp1 = find( num_list(Array_temp,4) ~= max(num_list(Array_temp,4)) );
            num_list(temp(Array_temp1),8)=num_list(temp(Array_temp1),8)+1;
            temp_Spike_point(2,temp(Array_temp1))=temp_Spike_point(2,temp(Array_temp1))+1;
            %Import points
            Array_temp1 = find( num_list(Array_temp,4) == max(num_list(Array_temp,4)) );
            num_list(temp(Array_temp1),8)=num_list(temp(Array_temp1),8)+nan;
            temp_Spike_point(2,temp(Array_temp1))=temp_Spike_point(2,temp(Array_temp1))+nan;
        end
        clear Array_temp Array_temp1
    end
    clear i
end