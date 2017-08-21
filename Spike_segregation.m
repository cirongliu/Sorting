%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Waveform extraction
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function [Single_spikes, Interval_single, Overlap_spikes] = Spike_segregation(spikes, Raw_data)
function [Single_spikes, Overlap_spikes] = Spike_segregation(spikes, Raw_data)
%% Filter by Spike interval
Spike_interval = round( diff(spikes(:,2))*1000 );%%transfer from ms to us
spikes(1,5) = NaN;
spikes(2:end,5) = Spike_interval;
clear Spike_interval
Number = find(spikes(:,5)<=1000);%Refractory period: <=1000us
for i = 1:length(Number)
    Row1 = Number(i);
    if sign(spikes(Row1-1,3))~=sign(spikes(Row1,3))%Opposite indicator
        temp = find( spikes(Row1-1:Row1,4) == max(spikes(Row1-1:Row1,4)) );
        if temp == 1
            if i == 1 
                spikes(Row1-1,6)=nan;%Peak
            else
                if spikes(Row1-1,6) == 0 
                    spikes(Row1-1,6)=nan;%Peak
                end
            end
            spikes(Row1,6)=1;%Useless
        else
            if i == 1
                spikes(Row1-1,6)=1;%Useless
            else
                if spikes(Row1-1,6) == 0
                    spikes(Row1-1,6)=1;%Useless
                end
            end
            spikes(Row1,6)=nan;%Peak
        end 
        clear temp
    else
        if Row1 == 2
            spikes(Row1-1,6)=nan;
            spikes(Row1,6)=nan;
        else
            spikes(Row1,6)=nan;
        end
    end
end
clear Number

Number = find(spikes(:,5)<=1000);%Refractory period: <=1000us
if length(Number) > 1
    temptemp = [Number,[NaN;diff(Number)]];
    Rows1 = find(temptemp(:,2)==1);
    Rows2 = find(temptemp(:,2)==1)-1;
    Rows_continue1 = temptemp(Rows1(:,1),1);
    Rows_continue1(:,2) = 0;
    Rows_continue2 = temptemp(Rows2(:,1),1);
    Rows_continue2(:,2) = 0;
    Rows_continue3 = Rows_continue2 -1;
    Rows_continue3(:,2) = NaN;
    
    overlapping_segment = unique(sortrows([Rows_continue1;Rows_continue2;Rows_continue3]),'rows');
    spikes(overlapping_segment(:,1),7) = overlapping_segment(:,2);
    
    
    %% Overlapping spikes
    Overlap_spikes = spikes(overlapping_segment(:,1),:);
    %Remove some spikes
    temp = unique(Overlap_spikes(:,1));
    if length(temp) ~= size(Overlap_spikes,1)
        for i = 1:length(temp)
            aa = find(Overlap_spikes(:,1)==temp(i));
            if length(aa)>1
                Overlap_spikes(aa(2:end),:)=[];
            end
        end
    end
    clear temp i aa
    %% Single unit
    spikes(overlapping_segment(:,1),:) = [];
    Single_spikes = spikes;
    Useless_spike = find( Single_spikes(:,6)==1 & Single_spikes(:,7)==0 );
    Single_spikes(Useless_spike,:) = [];
    %Remove some spikes
    temp = unique(Single_spikes(:,1));
    if length(temp) ~= size(Single_spikes,1)
        for i = 1:length(temp)
            aa = find(Single_spikes(:,1)==temp(i));
            if length(aa)>1
                Single_spikes(aa(2:end),:)=[];
            end
        end
    end
    clear temp i aa
else
    %% Single unit
    spikes(Number,:)=[];
    Single_spikes = spikes;  
    Single_spikes(:,7)=0;
    Overlap_spikes = [];
end

