%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Spike detection and waveform extraction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
close all

load([Datafile,'.mat'])%Load dataset

% Threshold detection
for i = 1:length(WB_band_sorting)%For each channel
    Raw_data = WB_band_sorting{i};% 1. Timecourse; 2. MUA signal; 3. Raw signal
    %Threshold, Bi-threshold
    WB = Raw_data(:,2);
    temp_WB = abs(WB);
    Thr = 4*median( temp_WB )/0.6745;%threshold
    
    %% threshold --- spike detection
    tic;
    waveform_temp1 = find(temp_WB >= Thr);
    waveform_temp1 = [waveform_temp1,Raw_data(waveform_temp1,1),Raw_data(waveform_temp1,2),temp_WB(waveform_temp1)];%Index, timepoints, voltages, absolute voltage
    clear temp_WB
    Spikes = Peak_filter(waveform_temp1);%Possible spike candidates    
    clear waveform_temp1
    timer_detect(i,:) = toc
    
    %% Preliminary classification --- single spike and overlapping spike
    tic;
    [Single_spikes, Overlap_spikes] = Spike_segregation(Spikes, Raw_data);
    timer_Preliminary(i,:) = toc
    
    %% Structure
    if size(Single_spikes,1) >= 20%Having enough single units, we can start to make template
         
        %% Template making
        tic
        % Waveform extraction
        Sampling_Freq = 40;%KHz
        [Single_spikes, ~, Waveform_single] = Waveform_extraction(Single_spikes, Raw_data, Sampling_Freq);
        % Filter some spikes
        [Single_spikes,Waveform_single] = Waveform_optimization(Single_spikes,Waveform_single,'Single');
        
        % Feature recognition
        [Spikes_S, Spikes_M, Spikes_S_wave, Spikes_M_wave] = Waveform_sorting(Waveform_single, Single_spikes);
        % Template Making
        wave_template = [];
        Num = unique(Spikes_S(:,2));
        for k = 1:length(Num)
            temp = find(Spikes_S(:,2) == Num(k));
            wave_template{k}(1,:) = mean(Spikes_S_wave(temp,:));
            wave_template{k}(2,:) = std(Spikes_S_wave(temp,:));
        end
        clear k temp Num
        Single_spikes(:,5:7) = [];%Make sure there are 5 columns
        Spikes_S_time = [Single_spikes(Spikes_S(:,1),:),Spikes_S(:,2)];%Conservative Template 
        Spikes_M_time = [Single_spikes(Spikes_M(:,1),:),Spikes_M(:,2)];
        clear Spikes_S Spikes_M
        
        timer_Template(i,:) = toc
        
        %% Overlapping
        tic
        if length(Overlap_spikes) >= 20
            
            [Possible_overlapping_Spikes] = Overlapping_Waveform_extraction(Overlap_spikes, Raw_data, wave_template, Sampling_Freq);
            % After segmentation, we organize new form
            [Possible_overlapping_Spikes, ~, Waveform_Overlapping] = Waveform_extraction(Possible_overlapping_Spikes, Raw_data, Sampling_Freq);
            Possible_overlapping_Spikes(:,6) = [];%Order,Make sure there are 5 columns
            % Filter some spikes
            [Possible_overlapping_Spikes,Waveform_Overlapping] = Waveform_optimization(Possible_overlapping_Spikes,Waveform_Overlapping,'Overlapping');
            
            for k = 1:length(wave_template)
                temp = find(Possible_overlapping_Spikes(:,end) == k);%Label mark
                if isempty(temp) ~= 1
                    Spikes_S_time = [Spikes_S_time;Possible_overlapping_Spikes(temp,:)];
                    Spikes_S_wave = [Spikes_S_wave;Waveform_Overlapping(temp,:)];
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    aa = Waveform_Overlapping(temp,:);
                    figure,
                    for j = 1:size(aa,1)
                        hold on
                        plot(aa(j,:));
                        hold off
                    end
                    clear j aa
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    Possible_overlapping_Spikes(temp,:) = [];
                    Waveform_Overlapping(temp,:) = [];
                end
            end
            clear k temp
            Spikes_M_time = [Spikes_M_time;Possible_overlapping_Spikes];
            Spikes_M_wave = [Spikes_M_wave;Waveform_Overlapping];
            clear Possible_overlapping_Spikes Waveform_Overlapping
            
            
            %% Sorting again
            [Spikes_S, Spikes_M_time_temp, Spikes_S_wave, Spikes_M_wave_temp] = Waveform_sorting_overlapping(Spikes_S_wave, Spikes_S_time);% After segment of the overlapping spikes
            % Rebuild MUA
            Spikes_M_time = [Spikes_M_time;Spikes_M_time_temp];
            Spikes_M_wave = [Spikes_M_wave;Spikes_M_wave_temp];
            clear Spikes_M_time_temp Spikes_M_wave_temp Spikes_S_time
            
            % MUA sorting
            [Spikes_M_time, Spikes_M_wave] = Waveform_optimization(Spikes_M_time, Spikes_M_wave, 'Single');
            [Spikes_S_temp, Spikes_M_final, Spikes_S_wave_temp, Spikes_M_wave_final] = Waveform_sorting_MUA(Spikes_M_wave, Spikes_M_time);% After segment of the overlapping spikes
%             Spikes_S_temp(:,5) = NaN;
            clear Spikes_M_wave Spikes_M_time
            
            if isempty(Spikes_S_temp) ~= 1
                Spikes_S = [Spikes_S;Spikes_S_temp];
                Spikes_S_wave = [Spikes_S_wave;Spikes_S_wave_temp];
            end
            clear Spikes_S_wave_temp Spikes_S_temp
        else
            Spikes_S = Single_spikes;
            Spikes_S_wave = Waveform_single;
        end
        
        timer_template_matching(i,:) = toc;
        
        %% Final isolated single unit sorting
        tic
        
        [Spikes_S_final_temp, Spikes_M_temp, Spikes_S_wave_final_temp, Spikes_M_wave_temp] = Waveform_sorting_S(Spikes_S_wave, Spikes_S);% After segment of the overlapping spikes
        
        if isempty(Spikes_M_temp) ~= 1
            if exist('Spikes_M_final')==1
                Spikes_M_req{i} = [Spikes_M_final;Spikes_M_temp];
                Spikes_M_wave_req{i} = [Spikes_M_wave_final;Spikes_M_wave_temp];
            else
                Spikes_M_req{i} = Spikes_M_temp;
                Spikes_M_wave_req{i} = Spikes_M_wave_temp;
            end
        else
            Spikes_M_req{i} = [];
            Spikes_M_wave_req{i} = [];
        end
        clear Spikes_M_wave_temp Spikes_M_temp
        
        if isempty(Spikes_S_final_temp) ~= 1
            Spikes_S_final{i} = Spikes_S_final_temp;
            Spikes_S_wave_final{i} = Spikes_S_wave_final_temp;
        else
            Spikes_S_final{i} = [];
            Spikes_S_wave_final{i} = [];
        end
        clear Spikes_S_final_temp Spikes_S_wave_final_temp
        
    else
        Spikes_S_final{i} = [];
        Spikes_S_wave_final{i} = [];
        Spikes_M_req{i} = [];
        Spikes_M_wave_req{i} = [];
    end
    timer_final(i,:) = toc
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Save each channel
    disp('Saving data to .mat file');
    save([Datafile,'_Result_',num2str(i)],'Spikes_S_final','Spikes_S_wave_final','Spikes_M_req','Spikes_M_wave_req','-v7.3');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function Spikes = Peak_filter(data)%Try to get reliable spikes
%Preprocess to the segments
temp = diff(data(:,1));
temp = [1;temp];
segment = find(temp~=1);
segment = [1;segment];
clear temp
aa = find(segment == size(data,1));
if isempty(aa)~=1
    segment(aa) = [];
end
clear aa
%For each segment to filter spikes    
label = [];
parfor i=1:length(segment) 
    if i~=length(segment)%Index, timepoints, voltages, absolute voltage
        temp = data(segment(i):(segment(i+1)-1),:);
        row = find(temp(:,4)==max(temp(:,4)));
        temp(row,5) = 1;
    else
        if segment(i)~=size(data,1)
            temp = data(segment(i):end,:);
            row = find(temp(:,4)==max(temp(:,4)));
            temp(row,5) = 1;
        end
    end
    label = [label;temp];
end
%
Spikes = label(find(label(:,5)==1),:);%Efficient spikes
% Spikes(:,1)=[];
Spikes(:,5)=[];

function [Spikes_time, Spike_waveform] = Waveform_optimization(Spikes_time, Spike_waveform, Type)%Try to remove redudent period

if strcmp (Type, 'Single') == 1
    temp = Spike_waveform;

else%Because of overlapping, try to shorten the length of waveform
    temp = Spike_waveform;
    temp(:,30:end) = [];%End other information
end
Row = [];
for k = 1:size(temp,1)
    if Spikes_time(k,3)<0 %negative value
        if temp(k,17) ~= min( temp(k,:) )
            Row = [Row;k];
        end
    else% Positive value
        if temp(k,17) ~= max( temp(k,:) )
            Row = [Row;k];
        end
    end
end
Spikes_time(Row,:) = [];
Spike_waveform(Row,:) = [];%

function [Spikes_S, Spikes_M, Spikes_S_wave, Spikes_M_wave] = Waveform_sorting(Spikes_waveform, Spikes_time)%Try to make template
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,
for j = 1:size(Spikes_waveform,1)
    hold on
    plot(Spikes_waveform(j,:));%Only select limited period
    hold off
end
clear j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fea_SPE = spe(Spikes_waveform,3);

%% SPE: pattern representation
figure,
subplot(1,3,1)
plot(Fea_SPE(:,1),Fea_SPE(:,2),'r.');
grid off
box off
subplot(1,3,2)
plot(Fea_SPE(:,1),Fea_SPE(:,3),'r.');
grid off
box off
subplot(1,3,3)
plot(Fea_SPE(:,2),Fea_SPE(:,3),'r.');
grid off
box off

% %% LPP: pattern representation
% [Fea_LPP,~] = lpp(Spikes_waveform,3);
% %%%%%%%%%%%%%%% %View from 3-dimension
% figure,
% subplot(1,3,1)
% plot(Fea_LPP(:,1),Fea_LPP(:,2),'r.');
% grid off
% box off
% subplot(1,3,2)
% plot(Fea_LPP(:,1),Fea_LPP(:,3),'r.');
% grid off
% box off
% subplot(1,3,3)
% plot(Fea_LPP(:,2),Fea_LPP(:,3),'r.');
% grid off
% box off
%  
% %% PCA: pattern representation
% [Fea_PCA,PCA_value] = pca(Spikes_waveform,3);
% %View from 3-dimension
% figure,
% subplot(1,3,1)
% plot(Fea_PCA(:,1),Fea_PCA(:,2),'r.');
% grid off
% box off
% subplot(1,3,2)
% plot(Fea_PCA(:,1),Fea_PCA(:,3),'r.');
% grid off
% box off
% subplot(1,3,3)
% plot(Fea_PCA(:,2),Fea_PCA(:,3),'r.');
% grid off
% box off

%% Spike sorting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fea_SPE = Fea_SPE*10^2;
if size(Fea_SPE,1)<1000
    opts.p = size(Fea_SPE,1);
else
    opts.p = 1000;
end
NumCluster = [1:7];
%OLPP
Value_SPE = [];
% parfor k = 1:length(NumCluster)
for k = 1:length(NumCluster)
    temp = LSC(Fea_SPE,k,opts);
    Value_SPE = [Value_SPE,temp];
end
clear temp k

%%%%%% CalinskiHarabasz Index
eva = evalclusters(Fea_SPE,Value_SPE,'Gap');%Spectral-cluster
% eva = evalclusters(Fea_SPE,Value_SPE,'CalinskiHarabasz');%Spectral-cluster
% eva = evalclusters(Fea_SPE,Value_SPE,'Silhouette');%Spectral-cluster
Max_gap_SPE = eva.OptimalK;
label_cluster_SPE = Value_SPE(:,Max_gap_SPE);

Text_function = 'try to find units and make them as the templates';
[Template_cluster,Template_label_cluster] = Modulation(Max_gap_SPE,label_cluster_SPE,Fea_SPE,Spikes_waveform,Spikes_time,Text_function);%GUI
%%%%%% Result ouput
Spikes_S=[];Spikes_S_wave=[];
Spikes_M=[];Spikes_M_wave=[];
parfor i = 1:length(Template_label_cluster)
    if strcmp( Template_label_cluster(i),'S' ) == 1%Single unit
        Spikes_S = [Spikes_S;[i,Template_cluster(i)]];%1 col, which row; 2 col, which label
        Spikes_S_wave = [Spikes_S_wave;Spikes_waveform(i,:)];
    else%Multi unit
        Spikes_M = [Spikes_M;[i,Template_cluster(i)]];
        Spikes_M_wave = [Spikes_M_wave;Spikes_waveform(i,:)];
    end
end
if isempty(Spikes_M) ~= 1
    Spikes_M(:,2) = NaN;
end
clear i
%Re-organize
temp = unique(Spikes_S(:,2));
for i = 1:length(temp)
    bb = Spikes_S;
    aa = find(bb(:,2) == temp(i));
    Spikes_S(aa,2) = i;
end
clear i aa bb
close all

function [Spikes_single, Spikes_Multi, Spikes_single_waveform, Spikes_multi_waveform] = Waveform_sorting_overlapping(Spikes_waveform, Spikes_time)%Try to isolate overlapping 
Level_template = unique(Spikes_time(:,5));%Number of template
Spikes_single = [];Spikes_single_waveform = [];
Spikes_Multi = [];Spikes_multi_waveform = [];
for kk = 1:length(Level_template)
    temp = find(Spikes_time(:,5) == Level_template(kk));
    Spikes_time_temp = Spikes_time(temp,:);
    Spikes_waveform_temp = Spikes_waveform(temp,:);%Waveform
%     Spikes_waveform_temp(:,30:end) = [];%Remove reduent waveform
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure,
    for j = 1:size(Spikes_waveform_temp,1)
        hold on
        plot(Spikes_waveform_temp(j,1:30));%Only select limited period
        hold off
    end
    clear j
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    
    % Sorting again
    Fea_SPE = spe(Spikes_waveform_temp(:,1:30),3);
    
%     %%%%SPE: pattern representation %%%%%%%%%%%%%%%%%%%%%%%%%%
%     figure,
%     subplot(1,3,1)
%     plot(Fea_SPE(:,1),Fea_SPE(:,2),'r.');
%     grid off
%     box off
%     subplot(1,3,2)
%     plot(Fea_SPE(:,1),Fea_SPE(:,3),'r.');
%     grid off
%     box off
%     subplot(1,3,3)
%     plot(Fea_SPE(:,2),Fea_SPE(:,3),'r.');
%     grid off
%     box off
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %% Spike sorting
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Fea_SPE = Fea_SPE*10^2;
    if size(Fea_SPE,1)<1000
        opts.p = size(Fea_SPE,1);
    else
        opts.p = 1000;
    end
    NumCluster = [1:7];
    %OLPP
    Value_SPE = [];
    % parfor k = 1:length(NumCluster)
    for k = 1:length(NumCluster)
        temp = LSC(Fea_SPE,k,opts);
        Value_SPE = [Value_SPE,temp];
    end
    clear temp k
    
    %%%%%% CalinskiHarabasz Index
    eva = evalclusters(Fea_SPE,Value_SPE,'CalinskiHarabasz');%Spectral-cluster
    % eva = evalclusters(Fea_SPE,Value_SPE,'Silhouette');%Spectral-cluster
    Max_gap_SPE = eva.OptimalK;
    label_cluster_SPE = Value_SPE(:,Max_gap_SPE);
    
    Text_function = 'Sorting under Template matching';
    [Template_cluster,Template_label_cluster] = Modulation(Max_gap_SPE,label_cluster_SPE,Fea_SPE,Spikes_waveform_temp(:,1:30),Spikes_time_temp,Text_function);
    %%%%%% Result ouput
    Spikes_S=[];Spikes_M=[];
    parfor i = 1:length(Template_label_cluster)
        if strcmp( Template_label_cluster(i),'S' ) == 1%Single unit
            Spikes_S = [Spikes_S;i];%1 col, which row;
        else%Multi unit
            Spikes_M = [Spikes_M;i];
        end
    end
    clear i
    if isempty(Spikes_S) ~= 1
        temp = Spikes_time_temp(Spikes_S,:);
        temp(:,5) = Template_cluster(Spikes_S,:);
        temp(:,6) = kk;
        Spikes_single = [Spikes_single;temp];
        Spikes_single_waveform = [Spikes_single_waveform;Spikes_waveform_temp(Spikes_S,:)];
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure,
        aa = Spikes_waveform_temp(Spikes_S,:);
        for j = 1:size(aa,1)
            hold on
            plot(aa(j,:));
            hold off
        end
        clear j aa
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    if isempty(Spikes_M) ~= 1
        Spikes_Multi = [Spikes_Multi;Spikes_time_temp(Spikes_M,:)];
        Spikes_multi_waveform = [Spikes_multi_waveform;Spikes_waveform_temp(Spikes_M,:)];
        Spikes_Multi(:,5) = NaN; 
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure,
        aa = Spikes_waveform_temp(Spikes_M,:);
        for j = 1:size(aa,1)
            hold on
            plot(aa(j,:));
            hold off
        end
        clear j aa
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    clear i Spikes_M Spikes_S Spikes_waveform_temp Spikes_time_temp  Template_label_cluster Template_cluster
    close all
end
% Re-organize
Num = 0;
Spikes_single_revise = [];
Spikes_single_waveform_revise = [];
for kk = 1:length(Level_template)
    temp = Spikes_single( find( Spikes_single(:,6) == kk ), : );%
    temp_wave = Spikes_single_waveform( find( Spikes_single(:,6) == kk ), : );%
    row = unique(temp(:,5));
    for i = 1:length(row)
        temptemp = temp( find(temp(:,5) == row(i)),: );
        temptemp_wave = temp_wave( find(temp(:,5) == row(i)),: );
        Num = Num+1;
        temptemp(:,7) = Num;
        Spikes_single_revise = [Spikes_single_revise;temptemp];
        Spikes_single_waveform_revise = [Spikes_single_waveform_revise;temptemp_wave];
    end
    clear i row temptemp temp temptemp_wave
end
Spikes_single = [];
Spikes_single = Spikes_single_revise;
Spikes_single(:,5:6) = [];
Spikes_single_waveform = [];
Spikes_single_waveform = Spikes_single_waveform_revise;

function [Spikes_single, Spikes_Multi, Spikes_single_waveform, Spikes_multi_waveform] = Waveform_sorting_MUA(Spikes_waveform, Spikes_time)%Try to isolate single unit from MUA
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,
for j = 1:size(Spikes_waveform,1)
    hold on
    plot(Spikes_waveform(j,:));%Only select limited period
    hold off
end
clear j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sorting again
Fea_SPE = spe(Spikes_waveform,3);
%% Spike sorting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fea_SPE = Fea_SPE*10^2;
if size(Fea_SPE,1)<1000
    opts.p = size(Fea_SPE,1);
else
    opts.p = 1000;
end
NumCluster = [1:7];
%OLPP
Value_SPE = [];
% parfor k = 1:length(NumCluster)
for k = 1:length(NumCluster)
    temp = LSC(Fea_SPE,k,opts);
    Value_SPE = [Value_SPE,temp];
end
clear temp k

%%%%%% CalinskiHarabasz Index
eva = evalclusters(Fea_SPE,Value_SPE,'CalinskiHarabasz');%Spectral-cluster
% eva = evalclusters(Fea_SPE,Value_SPE,'Silhouette');%Spectral-cluster
Max_gap_SPE = eva.OptimalK;
label_cluster_SPE = Value_SPE(:,Max_gap_SPE);

Text_function = 'Sorting for Multi units';
[MUA_cluster,MUA_label_cluster] = Modulation(Max_gap_SPE,label_cluster_SPE,Fea_SPE,Spikes_waveform,Spikes_time,Text_function);
close all
%%%%% 
Spikes_S = [];Spikes_M=[];
parfor i = 1:length(MUA_label_cluster)
    if strcmp( MUA_label_cluster(i),'S' ) == 1%Single unit
        Spikes_S = [Spikes_S;i];%1 col, which row;
    else%Multi unit
        Spikes_M = [Spikes_M;i];
    end
end
clear i
% Single unit
if isempty(Spikes_S) ~= 1
    Spikes_single = Spikes_time(Spikes_S,:);
    Spikes_single_waveform = Spikes_waveform(Spikes_S,:);
    
    %Re-orgnize
    if isempty( isnan(Spikes_single(:,5)) ) ~= 1
        Spikes_single( find(isnan(Spikes_single(:,5))==1), 5) = 1000;
    end
    row = unique(Spikes_single(:,5));
    Num = 1010;
    Spikes_single_revise = [];
    Spikes_single_waveform_revise = [];
    for i = 1:length(row)
        temptemp = Spikes_single( find(Spikes_single(:,5) == row(i)),: );
        temptemp_wave = Spikes_single_waveform( find(Spikes_single(:,5) == row(i)),: );
        Num = Num+1;
        temptemp(:,5) = Num;
        Spikes_single_revise = [Spikes_single_revise;temptemp];
        Spikes_single_waveform_revise = [Spikes_single_waveform_revise;temptemp_wave];
    end
    clear i Num temptemp_wave Spikes_single Spikes_single_waveform
    Spikes_single = Spikes_single_revise;
    Spikes_single_waveform = Spikes_single_waveform_revise;
    clear Spikes_single_revise Spikes_single_waveform_revise row
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure,
    for j = 1:size(Spikes_single_waveform,1)
        hold on
        plot(Spikes_single_waveform(j,:));
        hold off
    end
    clear j
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Spikes_single =[];
    Spikes_single_waveform =[];
end
%MUA
if isempty(Spikes_M) ~= 1
    Spikes_Multi = Spikes_time(Spikes_M,:);
    Spikes_Multi(:,5) = nan;
    if size(Spikes_Multi,2)>5
        Spikes_Multi(:,6:end)=[];
    end
    Spikes_multi_waveform = Spikes_waveform(Spikes_M,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure,
    for j = 1:size(Spikes_multi_waveform,1)
        hold on
        plot(Spikes_multi_waveform(j,:));
        hold off
    end
    clear j
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Spikes_Multi = [];
    Spikes_multi_waveform = [];
end

function [Spikes_single, Spikes_Multi, Spikes_single_waveform, Spikes_multi_waveform] = Waveform_sorting_S(Spikes_waveform, Spikes_time)%Try to validata single unit
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure,
for j = 1:size(Spikes_waveform,1)
    hold on
    plot(Spikes_waveform(j,1:40));%Only select limited period
    hold off
end
clear j
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Sorting again
Fea_SPE = spe(Spikes_waveform(:,1:40),3);
%% Spike sorting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if size(Fea_SPE,1)<1000
    opts.p = size(Fea_SPE,1);
else
    opts.p = 1000;
end
NumCluster = [1:7];
%OLPP
Value_SPE = [];
% parfor k = 1:length(NumCluster)
for k = 1:length(NumCluster)
    temp = LSC(Fea_SPE,k,opts);
%     temp = kmeans(Fea_SPE,k);
    Value_SPE = [Value_SPE,temp];
end
clear temp k

%%%%%% Gap Index
eva = evalclusters(Fea_SPE,Value_SPE,'Gap');%Spectral-cluster
% eva = evalclusters(Fea_SPE,Value_SPE,'CalinskiHarabasz');%Spectral-cluster
% eva = evalclusters(Fea_SPE,Value_SPE,'Silhouette');%Spectral-cluster
Max_gap_SPE = eva.OptimalK;
label_cluster_SPE = Value_SPE(:,Max_gap_SPE);

Text_function = 'Sorting for the final single units';
[S_cluster,S_label_cluster] = Modulation(Max_gap_SPE,label_cluster_SPE,Fea_SPE,Spikes_waveform(:,1:40),Spikes_time,Text_function);
close all
%%%%% 
Spikes_S = [];Spikes_M=[];
parfor i = 1:length(S_label_cluster)
    if strcmp( S_label_cluster(i),'S' ) == 1%Single unit
        Spikes_S = [Spikes_S;i];%1 col, which row;
    else%Multi unit
        Spikes_M = [Spikes_M;i];
    end
end
clear i
% Single unit
if isempty(Spikes_S) ~= 1
    Spikes_single = Spikes_time(Spikes_S,:);
    Spikes_single(:,5) = S_cluster(Spikes_S,:);
    Spikes_single_waveform = Spikes_waveform(Spikes_S,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure,
    for j = 1:size(Spikes_single_waveform,1)
        hold on
        plot(Spikes_single_waveform(j,:));
        hold off
    end
    clear j
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    temp = unique( Spikes_single(:,5) );
    Num = 0;
    for ii = 1:length(temp)
        aa = find(Spikes_single(:,5) == temp(ii));
        Num = Num+1;
        Spikes_single(aa,6) = Num;
    end
    Spikes_single(:,5) = [];
else
    Spikes_single = [];
    Spikes_single_waveform = [];
end
% MUA
if isempty(Spikes_M) ~= 1
    Spikes_Multi = Spikes_time(Spikes_M,:);
    Spikes_Multi(:,5) = nan;
    if size(Spikes_Multi,2)>5
        Spikes_Multi(:,6:end)=[];
    end
    Spikes_multi_waveform = Spikes_waveform(Spikes_M,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure,
    for j = 1:size(Spikes_multi_waveform,1)
        hold on
        plot(Spikes_multi_waveform(j,:));
        hold off
    end
    clear j
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Spikes_Multi = [];
    Spikes_multi_waveform = [];
end

