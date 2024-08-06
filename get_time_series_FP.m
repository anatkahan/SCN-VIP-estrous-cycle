function get_time_series_FP
%% this function is used to analyse the time dependent time-series FP data, 10 minutes per hour
% Taken with SynapseTDT

%% get experimental information: 
%my_path='D:\DATA_Glab\fiberphotometry\';
my_path='Z:\Anat\DATA_Glab\fiberphotometry\';


%%%%%%%%%%%%%%%%%%
% Females. 1-6
ind=0;
ind =ind+1; mouse_info{ind}.ID='198R';mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1; mouse_info{ind}.sex='Female';
ind=ind+1; mouse_info{ind}.ID='200LL';mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1; mouse_info{ind}.sex='Female';%
ind=ind+1; mouse_info{ind}.ID='246RL';mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='Female';
ind=ind+1; mouse_info{ind}.ID='247RRL';mouse_info{ind}.side='L';mouse_info{ind}.load_new_trials=1 ;mouse_info{ind}.sex='Female';%
ind=ind+1; mouse_info{ind}.ID='259R'; mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1 ;mouse_info{ind}.sex='Female';
ind=ind+1; mouse_info{ind}.ID='260L'; mouse_info{ind}.side='L';mouse_info{ind}.load_new_trials=1; mouse_info{ind}.sex='Female';
ind=ind+1; mouse_info{ind}.ID='261RL'; mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1; mouse_info{ind}.sex='Female'; %X=8; % has 3 completed cycles, remove the last 8 sessions to look at complete 3 sets
%ind=ind+1; mouse_info{ind}.ID='313RL';
%mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1;
%mouse_info{ind}.sex='Female';% signal was very low

n_females=ind;
% % Males 7-12
% ind=ind+1;mouse_info{ind}.ID='262R'; mouse_info{ind}.side='R'; mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='Male'; %male
% ind=ind+1; mouse_info{ind}.ID='273RL'; mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='Male';%male
% ind=ind+1; mouse_info{ind}.ID='286R'; mouse_info{ind}.side='R'; mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='Male';%male
% ind=ind+1;mouse_info{ind}.ID='287L';mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='Male';%male not good
% ind=ind+1; mouse_info{ind}.ID='288RL'; mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='Male';%male
% ind=ind+1; mouse_info{ind}.ID='296R'; mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='Male';%male
% n_males=6;
% %OVX 13-16
% ind=ind+1; mouse_info{ind}.ID='247RRL_OVX';mouse_info{ind}.side='L';mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='OVX';% OVX
% ind=ind+1; mouse_info{ind}.ID='246RL_OVX'; mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='OVX';% OVX
% ind=ind+1;mouse_info{ind}.ID='261RL_OVX'; mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='OVX';% OVX
% ind=ind+1; mouse_info{ind}.ID='259R_OVX'; mouse_info{ind}.side='R';mouse_info{ind}.load_new_trials=1;mouse_info{ind}.sex='OVX'; % OVX
% n_OVX=4;



% get data for each mouse
% for idi=1:length(mouse_info)
%     [output{idi}] = get_time_series_FP_per_mouse (mouse_info{idi});
% end
% get data for each mouse
for idi=1:length(mouse_info)
    switch mouse_info{idi}.sex
        case 'Female'
            %save([my_path 'time_series_output' newf '_male_female_' mouse_info.ID], 'output')
            if ~exist ([my_path 'time_series_output_newF_general_' mouse_info{idi}.ID '.mat' ])
                disp(['get data for ' mouse_info{idi}.ID])
                output= get_time_series_FP_per_mouse_v2 (mouse_info{idi});
            else
                load([my_path 'time_series_output_newF_general_' mouse_info{idi}.ID '.mat' ])% load 'output'
            end
%         case {'Male','OVX'}
%             if ~exist ([my_path 'time_series_output_general_' mouse_info{idi}.ID '.mat' ])
%                 disp(['get data for ' mouse_info{idi}.ID])
%                 output = get_time_series_FP_per_mouse (mouse_info{idi});
%             else
%                 load([my_path 'time_series_output_general_' mouse_info{idi}.ID '.mat' ])%all_output{idi}=output;
%             end 
    end
  %  output.new_estrus_states=estrus_to_receptive(output.estrus_states);
    all_output{idi}=output;
end
clear output
output=all_output; 
clear all_output
%%%%%%%%%%%%%%%    
for i=1:length(mouse_info)
    mouse_info{i}.analysis_type='estrous_cycle';
end

estrus_states_titles={'P-2','P-1','P+0','P+1','P+2','OVX'};
[ALL_colors,color_ind]=get_estrus_colors(estrus_states_titles);


L=[];
for idi=1:length(mouse_info)
    L=[L size(output{idi}.mean_dF_over_all_samples,2)];
end


% arrange data
clear mean_dF_over_all_samples mean_dF_over_estrus_samples mean_dF_by_state mean_rate_over_estrus_samples
clear median_dF_over_all_samples median_dF_over_estrus_samples median_dF_by_state median_rate_over_estrus_samples
for idi=1:length(mouse_info)
    if size(output{idi}.mean_dF_over_estrus_samples,1)<6 % no OVX for this female- add nan arrays to match size 
        output{idi}.mean_dF_over_estrus_samples=cat(1, output{idi}.mean_dF_over_estrus_samples(:,:,1:min(L)),nan(1,24,min(L)));
        output{idi}.mean_rate_over_estrus_samples=cat(1,output{idi}.mean_rate_over_estrus_samples,nan(1,24));
        output{idi}.mean_width_over_estrus_samples=cat(1,output{idi}.mean_width_over_estrus_samples,nan(1,24));
        output{idi}.median_dF_over_estrus_samples=cat(1, output{idi}.median_dF_over_estrus_samples(:,:,1:min(L)),nan(1,24,min(L)));
        output{idi}.median_rate_over_estrus_samples=cat(1,output{idi}.median_rate_over_estrus_samples,nan(1,24));
        output{idi}.median_width_over_estrus_samples=cat(1,output{idi}.median_width_over_estrus_samples,nan(1,24));
        output{idi}.mean_dF_by_state=cat(1,output{idi}.mean_dF_by_state,nan(1,24));
        output{idi}.median_dF_by_state=cat(1,output{idi}.median_dF_by_state,nan(1,24));
        output{idi}.median_dF_dark_light=[output{idi}.median_dF_dark_light nan nan];% add two nans for OVX light and dark phasese
        output{idi}.median_rate_dark_light=[output{idi}.median_rate_dark_light nan nan];% add two nans for OVX light and dark phasese
        output{idi}.median_width_dark_light=[output{idi}.median_width_dark_light nan nan];% add two nans for OVX light and dark phasese
    end
    mean_dF_over_all_samples(idi,:,:)=output{idi}.mean_dF_over_all_samples(:,1:min(L));
    median_dF_over_all_samples(idi,:,:)=output{idi}.median_dF_over_all_samples(:,1:min(L));
    %a=output{idi}.mean_dF_over_estrus_samples;
    mean_dF_over_estrus_samples(idi,:,:,:)=output{idi}.mean_dF_over_estrus_samples(:,:,1:min(L));% ID, estrus, hour, time
    median_dF_over_estrus_samples(idi,:,:,:)=output{idi}.median_dF_over_estrus_samples(:,:,1:min(L));% ID, estrus, hour, time
    %b=output{idi}.mean_dF_by_state;
    mean_dF_by_state(idi,:,:)=output{idi}.mean_dF_by_state;% ID, estrus, hour
    mean_rate_over_estrus_samples(idi,:,:)=output{idi}.mean_rate_over_estrus_samples;% ID, estrus, hour
    mean_width_over_estrus_samples(idi,:,:)=output{idi}.mean_width_over_estrus_samples;% ID, estrus, hour

    median_dF_by_state(idi,:,:)=output{idi}.median_dF_by_state;% ID, estrus, hour
    median_rate_over_estrus_samples(idi,:,:)=output{idi}.median_rate_over_estrus_samples;% ID, estrus, hour
    median_width_over_estrus_samples(idi,:,:)=output{idi}.median_width_over_estrus_samples;% ID, estrus, hour
    
    med_dF_dark_light_over_estrus_samples(idi,:)=output{idi}.median_dF_dark_light;% ID, estrus/dark/light
    med_rate_dark_light_over_estrus_samples(idi,:)=output{idi}.median_rate_dark_light;% ID,  estrus/dark/light
    med_width_dark_light_over_estrus_samples(idi,:)=output{idi}.median_width_dark_light;% ID,  estrus/dark/light

end

% average property over specific time window
% results an array of ID, estrus
clear med_dF_over_estrus_samples_window med_rate_over_estrus_samples_window med_width_over_estrus_samples_window

CF=0.8/length(mouse_info);% color factor
full_time_estrus_states_titles={'P-2 D','P-1 L','P-1 D','P+0 L','P+0 D','P+1 L','P+1 D','P+2 L','P+2 D' 'OVX L', 'OVX D'};

% test the results
figure
subplot(1,3,1)
for idi=1:length(mouse_info)
    plot(med_dF_dark_light_over_estrus_samples(idi,:)); hold on
end
ylabel('df')
subplot(1,3,2)
for idi=1:length(mouse_info)
    plot(med_rate_dark_light_over_estrus_samples(idi,:)); hold on
end
ylabel('rate')
subplot(1,3,3)
for idi=1:length(mouse_info)
    plot(med_width_dark_light_over_estrus_samples(idi,:)); hold on
end
ylabel('width')


% statistics
figure; [p,tbl,stats] = kruskalwallis(med_dF_dark_light_over_estrus_samples,[],'off');
c_dF= multcompare(stats,'CType','hsd');% Tukey's honest significant difference criterion correction
figure; [p,tbl,stats_r] = kruskalwallis(med_rate_dark_light_over_estrus_samples,[],'off');
c_r = multcompare(stats_r,'CType','hsd');
figure; [p,tbl,stats_w] = kruskalwallis(med_width_dark_light_over_estrus_samples,[],'off');
c_w = multcompare(stats_w,'CType','hsd');
% plot the median dF and rate over light/dark
%dF

figure
subplot(2,1,1)
bar(nanmedian(med_dF_dark_light_over_estrus_samples)); hold on;
for idi=1:length(mouse_info) 
    lh1=plot(med_dF_dark_light_over_estrus_samples(idi,:),'-*'); hold on;
    set(lh1,'Color', [CF*idi CF*idi CF*idi])
end
xticklabels(full_time_estrus_states_titles)
ylabel ('median dF')
%rates
subplot(2,1,2)
bar(nanmean(med_rate_dark_light_over_estrus_samples)); hold on;
for idi=1:length(mouse_info) 
    lh1=plot(med_rate_dark_light_over_estrus_samples(idi,:),'-*'); hold on;
    set(lh1,'Color', [CF*idi CF*idi CF*idi])
end
xticklabels(full_time_estrus_states_titles)
ylabel ('median rates')
%width
% subplot(3,1,3)
% bar(median(med_width_dark_light_over_estrus_samples)); hold on;
% for idi=1:length(mouse_info) 
%     lh1=plot(med_width_dark_light_over_estrus_samples(idi,:),'-*'); hold on;
%     set(lh1,'Color', [CF*idi CF*idi CF*idi])
% end
% xticklabels(full_time_estrus_states_titles)
% ylabel ('median width')

%% plot specific window data
clear lh1
% shift data: original data taken at strating at 9am, while at 13 light
% turned off. Data shifted that 1am (Light) will be the start time
% setting the first one: 
dF_to_plot=median_dF_by_state;
rate_to_plot=median_rate_over_estrus_samples;
width_to_plot=median_width_over_estrus_samples;
dF_by_state_DL(:,1,:)=cat(3,dF_to_plot(:,1,17:24),dF_to_plot(:,1,1:16));
rate_over_estrus_samples_DL(:,1,:)=cat(3,rate_to_plot(:,1,17:24),rate_to_plot(:,1,1:16));
width_over_estrus_samples_DL(:,1,:)=cat(3,width_to_plot(:,1,17:24),width_to_plot(:,1,1:16));
% tricky point here- the first 4 hours are fake- to make them nan 
dF_by_state_DL(:,1,1:4)=nan;
rate_over_estrus_samples_DL(:,1,1:4)=nan;
width_over_estrus_samples_DL(:,1,1:4)=nan;
% set all the other estrous states: 
for esi=1:4
    dF_by_state_DL(:,esi+1,:)=cat(3,dF_to_plot(:,esi,17:24),dF_to_plot(:,esi+1,1:16));
    rate_over_estrus_samples_DL(:,esi+1,:)=cat(3,rate_to_plot(:,esi,17:24),rate_to_plot(:,esi+1,1:16));
    width_over_estrus_samples_DL(:,esi+1,:)=cat(3,width_to_plot(:,esi,17:24),width_to_plot(:,esi+1,1:16));
end
% dF_by_state_DL(:,6,1:16)=nan;
% rate_over_estrus_samples_DL(:,6,1:16)=nan;
% width_over_estrus_samples_DL(:,6,1:16)=nan;

% now shift the OVX to fit 24h DL
dF_by_state_DL(:,6,:)=cat(3,dF_to_plot(:,6,17:24),dF_to_plot(:,6,1:16));
rate_over_estrus_samples_DL(:,6,:)=cat(3,rate_to_plot(:,6,17:24),rate_to_plot(:,6,1:16));
width_over_estrus_samples_DL(:,6,:)=cat(3,width_to_plot(:,6,17:24),width_to_plot(:,6,1:16));

%window=[1:24];
for i=1:12
    window=[i:i+1];% LH onset
    med_dF_over_estrus_samples_window(:,:)=nanmedian(dF_by_state_DL(:,:,window),3);% ID, estrus, hour
    med_rate_over_estrus_samples_window(:,:)=nanmedian(rate_over_estrus_samples_DL(:,:,window),3);% ID, estrus, hour
    med_width_over_estrus_samples_window(:,:)=nanmedian(width_over_estrus_samples_DL(:,:,window),3);% ID, estrus, hour
    
    mean_dF_over_estrus_samples_window(:,:)=nanmean(dF_by_state_DL(:,:,window),3);% ID, estrus, hour
    mean_rate_over_estrus_samples_window(:,:)=nanmean(rate_over_estrus_samples_DL(:,:,window),3);% ID, estrus, hour
    mean_width_over_estrus_samples_window(:,:)=nanmean(width_over_estrus_samples_DL(:,:,window),3);% ID, estrus, hour
    
    %%% pay attaention here - median or mean?
    %med_rate_over_estrus_samples_window=med_rate_over_estrus_samples_window([1,2,4:6],:);
    %med_dF_over_estrus_samples_window=med_dF_over_estrus_samples_window([1,2,4:6],:);
    figure
    subplot(2,1,1)
    bar(nanmedian(mean_dF_over_estrus_samples_window)); hold on;
    for idi=1:size(med_dF_over_estrus_samples_window,1)
        lh1=plot(med_dF_over_estrus_samples_window(idi,:),'-*'); hold on;
        set(lh1,'Color', [CF*idi CF*idi CF*idi])
    end
    xticklabels(estrus_states_titles)
    ylim([0 25])
    ylabel ('median/mean dF')
    title(['window : ' num2str(window(1)) ' to ' num2str(window(end))])
    %rates
    subplot(2,1,2)
    bar(nanmedian(med_rate_over_estrus_samples_window)); hold on;
    for idi=1:size(med_rate_over_estrus_samples_window,1)
        lh1=plot(med_rate_over_estrus_samples_window(idi,:),'-*'); hold on;
        set(lh1,'Color', [CF*idi CF*idi CF*idi])
    end
    xticklabels(estrus_states_titles)
    ylim([0 1])
    ylabel ('median/mean rates')
    title(['window : ' num2str(window(1)) ' to ' num2str(window(end))])
    
    
    % statistics
%     figure; [p2,tbl,stats_df] = kruskalwallis(med_dF_over_estrus_samples_window,[],'off');
%     c_dF2= multcompare(stats_df,'CType','hsd');
%     figure; [p1,tbl,stats_r] = kruskalwallis(med_rate_over_estrus_samples_window,[],'off');
%     c_r2 = multcompare(stats_r,'CType','hsd');
end
% compare time windows
window1=[1:3];
window2=[6:8];
window3=[10:12];
mean_rate_over_estrus_samples_window1(:,:)=nanmean(rate_over_estrus_samples_DL(:,:,window1),3);% ID, estrus, hour
mean_rate_over_estrus_samples_window2(:,:)=nanmean(rate_over_estrus_samples_DL(:,:,window2),3);% ID, estrus, hour
mean_rate_over_estrus_samples_window3(:,:)=nanmean(rate_over_estrus_samples_DL(:,:,window3),3);% ID, estrus, hour
figure
% window 1
subplot(3,1,1)
bar(nanmedian(mean_rate_over_estrus_samples_window1)); hold on;
for idi=1:size(mean_rate_over_estrus_samples_window1,1) 
    lh1=plot(mean_rate_over_estrus_samples_window1(idi,:),'-*'); hold on;
    set(lh1,'Color', [CF*idi CF*idi CF*idi])
end
xticklabels(estrus_states_titles)
ylabel ('median/mean rates')
ylim([0 0.8])
title(['window : ' num2str(window1(1)) ' to ' num2str(window1(end))])
%window 2
subplot(3,1,2)
bar(nanmedian(mean_rate_over_estrus_samples_window2)); hold on;
for idi=1:size(mean_rate_over_estrus_samples_window2,1) 
    lh1=plot(mean_rate_over_estrus_samples_window2(idi,:),'-*'); hold on;
    set(lh1,'Color', [CF*idi CF*idi CF*idi])
end
xticklabels(estrus_states_titles)
ylabel ('median/mean rates')
ylim([0 0.8])
title(['window : ' num2str(window2(1)) ' to ' num2str(window2(end))])
%window 3
subplot(3,1,3)
bar(nanmedian(mean_rate_over_estrus_samples_window3)); hold on;
for idi=1:size(mean_rate_over_estrus_samples_window3,1) 
    lh1=plot(mean_rate_over_estrus_samples_window3(idi,:),'-*'); hold on;
    set(lh1,'Color', [CF*idi CF*idi CF*idi])
end
xticklabels(estrus_states_titles)
ylabel ('median/mean rates')
ylim([0 0.8])
title(['window : ' num2str(window3(1)) ' to ' num2str(window3(end))])
%%%

% plot rates
clear tmp
figure

for sti=1:length(estrus_states_titles)
    subplot(length(estrus_states_titles),1,sti)
    for idi=1:length(mouse_info)
        tmp(:)=rate_over_estrus_samples_DL(idi,sti,:);
        lh=plot(tmp); hold on
        set(lh,'Color', ALL_colors(color_ind(sti),:))
    end
    clear tmp
    ylim([0 0.75])
    title ('rates')
end
% plot width
figure
clear tmp

for sti=1:length(estrus_states_titles)
    subplot(length(estrus_states_titles),1,sti)
    for idi=1:length(mouse_info)
        tmp(:)=width_over_estrus_samples_DL(idi,sti,:);
        lh=plot(tmp); hold on
        set(lh,'Color', ALL_colors(color_ind(sti),:))
    end
    clear tmp
    ylim([0 30])
    title ('width')
end




figure
% plot mean rates

for sti=1:length(estrus_states_titles)
    subplot(length(estrus_states_titles),1,sti)
    tmp(:)=nanmedian(rate_over_estrus_samples_DL(:,sti,:));
    lh=bar(tmp); hold on
    set(lh,'FaceColor', ALL_colors(color_ind(sti),:))
    for idi=1:length(mouse_info)
        tmp(:)=rate_over_estrus_samples_DL(idi,sti,:);
        lh=plot(tmp,'*'); hold on
        set(lh,'Color', [0.7 0.7 0.7])
    end
    title (estrus_states_titles{sti})
    ylabel('median rates')
    ylim([0 0.75])
end
legend (estrus_states_titles)

% plot mean width
figure
subplot(2,1,2)
for sti=1:length(estrus_states_titles)
    subplot(length(estrus_states_titles),1,sti)
    tmp(:)=nanmedian(width_over_estrus_samples_DL(:,sti,:));
    lh=bar(tmp); hold on
    set(lh,'FaceColor', ALL_colors(color_ind(sti),:))
    for idi=1:length(mouse_info)
        tmp(:)=width_over_estrus_samples_DL(idi,sti,:);
        lh=plot(tmp,'*'); hold on
        set(lh,'Color', [0.7 0.7 0.7])
    end
    title (estrus_states_titles{sti})
    ylabel('median width')
    ylim([0 9])
end


% plot dF
figure
for sti=1:length(estrus_states_titles)
    subplot(length(estrus_states_titles),1,sti)
    for idi=1:length(mouse_info)
    tmp(:)=dF_by_state_DL(idi,sti,:);
    lh=plot(tmp); hold on
    set(lh,'Color', ALL_colors(color_ind(sti),:))
    end
    ylim([0 15])
    ylabel ('dF')
    title (estrus_states_titles{sti})
end




% plot mean dF
figure
for sti=1:length(estrus_states_titles)
    subplot(length(estrus_states_titles),1,sti)
    tmp(:)=nanmedian(dF_by_state_DL(:,sti,:));
    lh=bar(tmp); hold on
    set(lh,'FaceColor', ALL_colors(color_ind(sti),:))
    for idi=1:length(mouse_info)
        tmp(:)=dF_by_state_DL(idi,sti,:);
        lh=plot(tmp,'*'); hold on
        set(lh,'Color', [0.7 0.7 0.7])
    end
    ylim([-0.5 15])
   ylabel ('median dF')
    title (estrus_states_titles{sti})
end

%legend (estrus_states_titles)


figure
for idi=1:length(mouse_info)
    for hi=1:size(output{idi}.mean_dF_over_all_samples,1)
        plot(output{idi}.mean_dF_over_all_samples(hi,:)+25*hi); hold on
    end
    ylim([-20 (hi+3)*25])
end
title('all mice')


% plot mean dF per ID
CF=0.8/length(mouse_info);% color factor
figure
for sti=1:length(estrus_states_titles)
    subplot(length(estrus_states_titles),1,sti)
     for idi=1:length(mouse_info)
        tmp(:)=dF_by_state_DL(idi,sti,:);
        lh=plot(tmp,'-*'); hold on
        set(lh,'Color', [CF*idi CF*idi CF*idi])
    end
    ylim([-2 20])
   title ([ estrus_states_titles{sti} ', median dF'])
end
legend (estrus_states_titles)


for sti=1:length(estrus_states_titles)
    all_df_by_state(sti,:)=nanmedian(dF_by_state_DL(:,sti,:));
    %histcounts(med_all_df_by_state(sti,:),24);
end





