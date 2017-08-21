%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Waveform extraction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Spike_List] = Overlapping_Waveform_extraction(Spikes, Raw_data, templates, Sampling_Freq)

%% Waveform
Waveform_segment = [];
Waveform_time_segment = [];
%Segment: the starting time and the end time
Time_begin = 0.8*Sampling_Freq;%0.8ms---->points*40
Time_end = 1*Sampling_Freq;%1ms------>point*40

%% Segment process I
%Every 10ms bin to run
temp = find(isnan( Spikes(:,7) ) == 1);
Num = 0;
for i = 1:length(temp)%Aligned on spike time
    Num = Num+1;
    if i ~= length(temp)
        Waveform_segment(Num,1) = Spikes(temp(i),1)-Time_begin;%starting index
        Waveform_segment(Num,2) = Spikes(temp(i+1)-1,1)+Time_end;%ending index
    else
        Waveform_segment(Num,1) = Spikes(temp(i),1)-Time_begin;%starting index
        Waveform_segment(Num,2) = Spikes(end,1)+Time_end;%ending index
    end

end
clear i
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CHECK = find( Waveform_segment(:,1) < 1 | Waveform_segment(:,2) > size(Raw_data,1) );
if isempty( CHECK ) ~= 1
    for i = 1:length(CHECK)
        if CHECK(i) ~= length(temp) 
            Spikes(temp(CHECK(i),:):temp(CHECK(i)+1)-1,8)=NaN;
        else
            Spikes(temp(CHECK(i),:):end,8) = NaN;
        end
    end
    Spikes(isnan(Spikes(:,8))==1,:)=[];
    Spikes(:,8)=[];
end
clear CHECK Waveform_segment
%Update again
temp = find(isnan( Spikes(:,7) ) == 1);
Num = 0;
for i = 1:length(temp)%Aligned on spike time
    Num = Num+1;
    if i ~= length(temp)
        Waveform_segment(Num,1) = Spikes(temp(i),1)-Time_begin;%starting index
        Waveform_segment(Num,2) = Spikes(temp(i+1)-1,1)+Time_end;%ending index
    else
        Waveform_segment(Num,1) = Spikes(temp(i),1)-Time_begin;%starting index
        Waveform_segment(Num,2) = Spikes(end,1)+Time_end;%ending index
    end

end
clear i Num
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Segment process II
Spike_List = [];
for i = 1:length(temp)%Aligned on spike time
    Waveform_time_segment(i,1) = Raw_data( Waveform_segment(i,1),1 );
    Waveform_time_segment(i,2) = Raw_data( Waveform_segment(i,2),1 );
    if i ~= length(temp)
        num_list = Spikes(temp(i):temp(i+1)-1,:);%Number of spikes list
    else
        num_list = Spikes(temp(i):end,:);%Number of spikes list
    end
    if size(num_list,1) ~= 1
        %Segment with overlapping
        Segment = Raw_data(Waveform_segment(i,1):Waveform_segment(i,2),1:2);%Time points and Voltage period
        Time_course = Segment(1,1):1/(Sampling_Freq*2):Segment(end,1);%New timecourse for mapping
        Waveform_process = mean([interp1(Segment(:,1),Segment(:,2),Time_course,'spline');interp1(Segment(:,1),Segment(:,2),Time_course,'spline')]);%2 times polaration
        %Recognition
        [possible_spikes] = Overlapping_recognition(Time_course, Waveform_process, templates, num_list, Sampling_Freq);
    else
        num_list(:,5:7) = [];
        num_list(:,5) = length(templates)+1;
        possible_spikes = num_list;
    end
    possible_spikes(:,6) = i;
    Spike_List = [Spike_List;possible_spikes];
    clear possible_spikes
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%     i
%     close all
%     %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end