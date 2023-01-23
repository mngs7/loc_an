clear;
warning off;

%time window for analysis
tw = [10; 30; 50; 150; 300];%different time windows
%stim conditions
cnd_name = ['p'; 'w'; 'm'; 's'];
%color index setting
color_idx = [239 69 40; %red
    142 181 66; %green
    250 208 26; %yellow
    63 155 177; %blue
    245 124 33; %orange
    93 63 182]; %purple
color_idx = color_idx/255;
col2 = [255 192 203; %pink
    255 0 0; %red
    135 206 235; %sky blue
    0 0 255; %blue
    152 251 152; %pale green
    0 128 0]; %green
col2 = col2/255;
%ask whether to remove previous results (remove and make)
%update = input('remove previous results? y/n   ','s');
update = 'y';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set analysis_dir
analysis_dir = uigetdir();
cd(analysis_dir)
%scan APT exp list (will also check if the files for analysis exist)
%list_read = input('Do you have a list for APT analysis? y/n   ','s');
list_read = 'n';
if list_read == 'y'
    APT_list_dir = uigetdir();
    cd(APT_list_dir)
    list = uigetfile();
    [~, ~, APT_list] = xlsread(list);
else
    cd(analysis_dir)
    APT_list1 = struct2cell(dir('VNC_*'))';
    APT_list2 = struct2cell(dir('VNC2_*'))';
    APT_list = [APT_list1; APT_list2];
    clear APT_list1 APT_list2
end
error_dir = uigetdir();
clear list_read
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%set target files to check before running analysis on each experiment
target_files = cell(1,8);
target_files{1} = 'apt.trk';
target_files{2} = 'indicatordata.mat';
target_files{3} = 'registered_trx.mat';
target_files{4} = 'du_ctr.mat';
target_files{5} = 'dist2wall.mat';
target_files{6} = 'dcenter.mat';
target_files{7} = 'dtheta.mat';
target_files{8} = 'scores_Walk.mat';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%process each experiment
for d = 1:length(APT_list(:,1))
    tic
    disp(['exp_dir ' num2str(d) '/' num2str(length(APT_list(:,1)))])
    disp(APT_list{d,1})
    if strfind(APT_list{d,1},'VNC_') > 0
        stn = 3; %2021 protocol
    elseif strfind(APT_list{d,1},'VNC2_') > 0
        stn = 10; %2022 protocol
    end
    cd(APT_list{d,1})
    exp_dir = pwd;
    file_checker = zeros(1,8);
    for k = 1:8 %check apt.trk and other mat files for analysis
        temp_check = exist(target_files{k},'file');
        if temp_check == 2
            file_checker(1,k) = 1;
        end
        clear temp_check
    end
    if sum(file_checker(1,1:8)) == 8 %execute analysis after checking file existence
        if update == 'y'
            if isfolder('results')
                rmdir('results', 's')
            end
        end
        mkdir results
        cd results
        save_dir = pwd;
        cd(exp_dir)
        load('apt.trk','-mat','pTrk')
        fly_num = size(pTrk,2); %extract fly number
        %check data length
        dl = zeros(1,fly_num);
        for id = 1:fly_num
            dl(id) = length(rmmissing(squeeze(pTrk{id}(1,1,:))));
        end
        if min(dl) == max(dl) %if data size is equal
            expname = APT_list{d,1}; %extract only exp name
            load('indicatordata.mat','indicatorLED') %get LED log
            load('registered_trx.mat','trx') %read ctrax perframe data
            %extract ctrax perframe parameters(forward v, angular v)
            load('du_ctr.mat','data');
            du_ctr = data;
            load('dist2wall.mat','data');
            dist2wall = data;
            load('dcenter.mat','data');
            dcenter = data;
            load('dtheta.mat','data');
            dtheta = data;
            %adjust data size (add one rwo for velocity: du_ctr, dtheta)
            for id = 1:fly_num
                du_ctr{id} = cat(2,0,du_ctr{id});%insert column to adjust data length
                dtheta{id} = cat(2,0,dtheta{id});%insert column to adjust data length
            end
            %load JAABA score
            if isfile('scores_Walk.mat')
                load('scores_Walk.mat','allScores')
            end
            %generate status_filter from JAABA or ctrx parameters
            status_filter = zeros(fly_num,stn*3);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %data extraction for all flies
            data = cell(fly_num,6);%collect data for all flies
            %rows:fly#,
            %columns:1.raw_vel,2.smoothed_vel,3.walk_status,
            %4.raw_phase,5.smoothed_phase
            %6.extracted positiins
            for id = 1:fly_num
                scale = trx(id).pxpermm;
                %extract pisitional information from pTrk (in .trk file)
                trk_posi = cell(1,9);
                trk_posi{1} = squeeze(pTrk{id}(17,:,:))';%left front leg
                trk_posi{2} = squeeze(pTrk{id}(16,:,:))';%left mid leg
                trk_posi{3} = squeeze(pTrk{id}(15,:,:))';%left hind leg
                trk_posi{4} = squeeze(pTrk{id}(12,:,:))';%right front leg
                trk_posi{5} = squeeze(pTrk{id}(13,:,:))';%right mid leg
                trk_posi{6} = squeeze(pTrk{id}(14,:,:))';%right hind leg
                trk_posi{7} = (squeeze(pTrk{id}(4,:,:))' + squeeze(pTrk{id}(5,:,:))' + squeeze(pTrk{id}(6,:,:))')/3;%Center of thorax
                trk_posi{8} = squeeze(pTrk{id}(5,:,:))';%left shoulder
                trk_posi{9} = squeeze(pTrk{id}(4,:,:))';%right shoulder
                data_length = length(rmmissing(trk_posi{1,1}));
                %calculate positions thx center aligned with center of origin
                for mk = 1:9 %subtract center position to align origin
                    trk_posi_r{mk} = trk_posi{mk} - trk_posi{7};
                end
                ori_vctr = (trk_posi_r{8} + trk_posi_r{9})/2; %orientation vector calculated from mid shoulder
                ang = atan2(ori_vctr(:,2),ori_vctr(:,1));
                ang = -ang + pi/2; %align angle orient to y axis
                %convert tracking points
                for mk = 1:9
                    rot_x = [cos(ang) -sin(ang)];
                    rot_y = [sin(ang) cos(ang)];
                    temp_x = sum(rot_x.*trk_posi_r{mk},2);
                    temp_y = sum(rot_y.*trk_posi_r{mk},2);
                    trk_posi_r2{mk} = [temp_x temp_y];
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %extract velocity (raw and smoothed)
                raw_vel = cell(1,9);
                s_vel = cell(1,9);
                for p = 1:9 %6legs + thorax center + L/R shoulders
                    %calculate marker position velocity(raw)
                    temp_vel = diff(trk_posi{1,p},1,1)*150/scale;
                    temp_vel = sqrt(temp_vel(:,1).^2 + temp_vel(:,2).^2);
                    temp_vel = cat(1,0,temp_vel);
                    raw_vel{1,p} = temp_vel;%assign velocity
                    s_vel{1,p} =  smoothdata(raw_vel{1,p},'movmedian',3);
                    clear temp_vel
                end
                data{id,1} =  raw_vel;
                data{id,2} =  s_vel;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %prepare walking status (o/1 for data length) from JAABA or
                %ctrax parameters
                walk_status = zeros(1,data_length);
                if isfile('scores_Walk.mat')
                    walk_status = allScores.postprocessed{1,id};
                else
                    %check ctrax perframe (for walk_status, preparation)
                    %this process required for the following status_filter generation (if there's no JAABA score)
                    %bin_vel = du_ctr{id} > 3.5;
                    temp_dist2wall = dist2wall{1,id};
                    temp_dcenter = dcenter{1,id};
                    for k = 1:data_length
                        if data{id,2}{1,7}(k) > 3.5 && data{id,2}{1,7}(k) < 100 %speed filter
                            %walk_status(k,id) = 1;
                            if temp_dist2wall(k) > 4 && temp_dcenter(k) > 3 %wall & social distance
                                walk_status(1,k) = 1;
                            end
                        end
                    end
                    %remove small fragments by calculating pulse width
                    [W, INITCROSS] = pulsewidth([walk_status 0]);
                    for k = 1:length(W)
                        if W(k) < 15
                            idx = INITCROSS(k)+0.5;
                            walk_status(1,idx:idx+W(k)-1) = 0;
                        end
                    end
                end
                data{id,3} = walk_status;
                %generate status filter from walk_status
                for stim = 1:stn*3
                    st = indicatorLED.startframe(stim) - 29;
                    ed = indicatorLED.startframe(stim);
                    temp = data{id,3}(1,st:ed);
                    %decide status
                    if sum(temp(1,1:30)) == 30 %if continuous 30 frames are all positive
                        status_filter(id,stim) = 1;
                    end
                    clear temp
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %swing/stance phase binarization
                %detect swing phase based on threshold
                th = [20 20 20 20 20 20];%threshold for each leg
                for rs = 1:2 %rs=1(raw), rs=2(smoothed)
                    phase = cell(3,6);%phase information (1:swing onset, 2:swing offset, 3:swing duration)
                    for leg = 1:6
                        %process for raw data
                        %find swing onset
                        uth = data{id,rs}{leg}' > th(leg);
                        uth_shift = [0, uth(1:end-1)];
                        %extract first row over threshold
                        uth_trig = uth.*xor(uth,uth_shift) == 1;
                        uth_trig2 = find(uth_trig);
                        %find index which have less than 7 frames to the next next index
                        uth_gap = diff(uth_trig2);
                        uth_short = find(uth_gap <= 7);
                        %remove the index next to the detected index (too-short-next, because of miss detection)
                        uth_trig2(uth_short+1) = [];
                        %search actual swing start point
                        uth2 = zeros(1,length(uth_trig2));
                        for p = 1:length(uth_trig2)
                            temp = uth_trig2(p);
                            if uth_trig2(p)>1 && data{id,rs}{leg}(uth_trig2(p)) > data{id,rs}{leg}(uth_trig2(p)-1)
                                temp = uth_trig2(p)-1;
                                if uth_trig2(p)>2 && data{id,rs}{leg}(uth_trig2(p)-1) > data{id,rs}{leg}(uth_trig2(p)-2)
                                    temp = uth_trig2(p)-2;
                                end
                            end
                            uth2(p) = temp;
                        end
                        phase{1,leg} = uth2;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %find swing offset (stance onset)
                        %detect points below threshold at first
                        dth = data{id,rs}{leg}' < th(leg);
                        dth_shift = [0, dth(1:end-1)];
                        %extract first row below threshold
                        dth_trig = dth.*xor(dth,dth_shift) == 1;
                        dth_trig2 = find(dth_trig);
                        %search swing offset from swing onset
                        dth2 = zeros(1,length(uth_trig2));
                        for p = 1:length(uth_trig2)-1
                            temp = dth_trig2(dth_trig2 > uth_trig2(p) & dth_trig2 <= uth_trig2(p+1));
                            temp = max(temp);
                            dth2(1,p) = temp;
                        end
                        temp = dth_trig2(dth_trig2 > uth_trig2(end));
                        if ~isempty(temp)
                            temp = max(temp);
                            dth2(end) = temp;
                        else
                            dth2(end) = length(data{id,rs}{leg}');
                        end
                        %search actual swing end point
                        dth2_2 = zeros(1,length(dth2));
                        for p = 1:length(dth2)
                            temp = dth2(p);
                            if dth2(p)<data_length && data{id,rs}{leg}(dth2(p)) > data{id,rs}{leg}(dth2(p)+1)
                                temp = dth2(p)+1;
                                if dth2(p)<data_length-1 && data{id,rs}{leg}(dth2(p)+1) > data{id,rs}{leg}(dth2(p)+2)
                                    temp = dth2(p)+2;
                                end
                            end
                            if p < length(dth2) && temp > uth_trig2(p+1) %if swing offset
                                temp = uth_trig2(p+1);
                            end
                            dth2_2(p) = temp;
                        end
                        %                     plot(data{id,rs}{leg})
                        %                     hold on
                        %                     plot(data{id,rs}{leg},'k*')
                        %                     plot(uth2,data{id,rs}{leg}(uth2),'m*')
                        %                     plot(dth2_2,data{id,rs}{leg}(dth2_2),'g*')
                        %                     yline(20,'-r')
                        phase{2,leg} = dth2_2;
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        %calculate length between swing phases(onset - onset)
                        phase{3,leg} = diff(uth2);
                    end
                    data{id,rs+3} = phase;%put phase in data cell
                    clear phase
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %store posision information
                data{id,6} = trk_posi;
                %store aligned position information
                data{id,7} = trk_posi_r2;
            end
            %save_swing timing
            sw_phase = data(:,4:5);
            cd(save_dir)
            save('sw_phase.mat', 'sw_phase', '-mat')
            %         clear trx
            %         clear pTrk

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % primal data analysis (body trajectories, fv, angv)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %trajectories
            %%%%%%%%%%%%%
            for tp = 1:5
                posi_thx = cell(1,2); %column1-2:raw/aligned
                %%% process for center of thorax
                %alignment (position:start-point, angle:end-point)
                posi_raw = cell(1,3); %column1-3:w/m/s
                posi_aligned = cell(2,3); %column1-3:w/m/s
                for level = 1:3 %3 level of stimuli
                    for n = 1:stn %stimuli times (3 or 10)
                        stim = n + stn*(level-1);
                        pre = indicatorLED.startframe(stim)-tw(tp);
                        on = indicatorLED.startframe(stim);
                        post = indicatorLED.startframe(stim)+tw(tp);
                        for id = 1:fly_num
                            if status_filter(id,stim) == 1
                                temp_posi{1} = data{id,6}{7}(pre:post,:);
                                posi_raw{1,level} = [posi_raw{1,level} temp_posi];
                                %calculate aligned values
                                temp_pre{1} = data{id,6}{7}(pre:on,:) - data{id,6}{7}(pre,:); %adjust starting position(pre)
                                temp_post{1} = data{id,6}{7}(on:post,:) - data{id,6}{7}(on,:); %adjust starting position(post)
                                %calculate angle at end point (pre) and align
                                ori_angle = atan2(temp_pre{1}(end,1),temp_pre{1}(end,2));
                                %ori_angle = -ori_angle + pi/2;
                                rot_mat = [cos(ori_angle) -sin(ori_angle); sin(ori_angle) cos(ori_angle)];
                                temp_pre{1} = (rot_mat * temp_pre{1}')';
                                %calculate angle at end point (post) and align
                                ori_angle = atan2(temp_post{1}(end,1),temp_post{1}(end,2));
                                %ori_angle = -ori_angle + pi/2;
                                rot_mat = [cos(ori_angle) -sin(ori_angle); sin(ori_angle) cos(ori_angle)];
                                temp_post{1} = (rot_mat * temp_post{1}')';
                                %add aligned data (pre and post)
                                posi_aligned{1,level} = [posi_aligned{1,level} temp_pre];
                                posi_aligned{2,level} = [posi_aligned{2,level} temp_post];
                            end
                        end
                    end
                end
                posi_thx{1} = posi_raw;
                posi_thx{2} = posi_aligned;
                %visualization
                if tp == 3  %plot figure for 50, 150 frames trajectories
                    thx_trj = figure('Position', [100, 100, 300, 700],'visible','off');
                    for level = 1:3
                        subplot(3,1,level)
                        for wf = 1:length(posi_thx{2}{1,level}) %overlay all walking fly
                            plot(posi_thx{2}{1,level}{wf}(1:end,1),posi_thx{2}{1,level}{wf}(1:end,2),'Color',col2(4,:))
                            hold on
                            plot(posi_thx{2}{2,level}{wf}(1:end,1),posi_thx{2}{2,level}{wf}(1:end,2),'Color',col2(2,:))
                            axis equal
                        end
                    end
                    %get fig image
                    cd(save_dir)
                    I = getframe(gcf);
                    I1 = frame2im(I);
                end
                %save
                cd(save_dir)
                save(['posi_thx_tw' num2str(tw(tp)) '.mat'], 'posi_thx', '-mat')
                %%%%%%%%%%%%%%%%%%%%%%%
                %foward velocity
                %%%%%%%%%%%%%%%%%%%%
                fv_raw = cell(1,3); %column1-3:w/m/s
                fv_mean = cell(2,3); %row1:pre, row2:post, column1-3:w/m/s
                for level = 1:3 %3 level of stimuli
                    for n = 1:stn %stimuli times (3 or 10)
                        stim = n + stn*(level-1);
                        pre = indicatorLED.startframe(stim)-tw(tp);
                        on = indicatorLED.startframe(stim);
                        post = indicatorLED.startframe(stim)+tw(tp);
                        for id = 1:fly_num
                            if status_filter(id,stim) == 1
                                y_fv = du_ctr{id}(1, pre:post)';
                                fv_raw{1,level} = [fv_raw{1,level} y_fv];
                            end
                        end
                        %collect mean value (pre-post)
                        for id = 1:fly_num
                            if status_filter(id,stim) == 1
                                pre_v = mean(du_ctr{id}(1,pre:on));
                                post_v = mean(du_ctr{id}(1,on:post));
                                fv_mean{1,level} = [fv_mean{1, level} pre_v];
                                fv_mean{2,level} = [fv_mean{2, level} post_v];
                            end
                        end
                    end
                end
                %visualization
                if tp == 3
                    fv_raw_fig = figure('Position', [100, 100, 300, 700],'visible','off');
                    for level = 1:3 %3 level of stimuli
                        subplot(3,1,level)
                        for n = 1:stn %stimuli times (3 or 10)
                            stim = n + stn*(level-1);
                            pre = indicatorLED.startframe(stim)-tw(tp);
                            post = indicatorLED.startframe(stim)+tw(tp);
                            t = -tw(tp)/150:1/150:tw(tp)/150;
                            for id = 1:fly_num
                                if status_filter(id,stim) == 1
                                    y_fv = du_ctr{id}(1, pre:post)';
                                    hold on
                                    plot(t,y_fv)
                                    xline(0,'r')
                                end
                            end
                        end
                        xlim([-tw(tp)/150 tw(tp)/150])
                        ylim([-20 60])
                        %get fig image
                        I = getframe(gcf);
                        I2 = frame2im(I);
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%
                %angular velocity
                %%%%%%%%%%%%%%%%%%%%
                angv_raw = cell(1,3); %column1-3:w/m/s
                angv_mean = cell(2,3); %row1:mean-pre, row2:mean-post, column1-3:w/m/s
                for level = 1:3 %3 level of stimuli
                    for n = 1:stn %10 times stimuli
                        stim = n + stn*(level-1);
                        pre = indicatorLED.startframe(stim)-tw(tp);
                        on = indicatorLED.startframe(stim);
                        post = indicatorLED.startframe(stim)+tw(tp);
                        for id = 1:fly_num
                            if status_filter(id,stim) == 1
                                y_angv = abs(dtheta{id}(1, pre:post)');
                                angv_raw{1,level} = [angv_raw{1,level} y_angv];
                            end
                        end
                        %collect mean value (pre-post)
                        for id = 1:fly_num
                            if status_filter(id,stim) == 1
                                pre_v = mean(abs(dtheta{id}(1,pre:on)));
                                post_v = mean(abs(dtheta{id}(1,on:post)));
                                angv_mean{1,level} = [angv_mean{1, level} pre_v];
                                angv_mean{2,level} = [angv_mean{2, level} post_v];
                            end
                        end
                    end
                end
                %visualization
                if tp == 3
                    angv_raw_fig = figure('Position', [100, 100, 300, 700],'visible','off');
                    for level = 1:3 %3 level of stimuli
                        subplot(3,1,level)
                        for n = 1:stn %stimuli times (3 or 10)
                            stim = n + stn*(level-1);
                            pre = indicatorLED.startframe(stim)-tw(tp);
                            post = indicatorLED.startframe(stim)+tw(tp);
                            t = -tw(tp)/150:1/150:tw(tp)/150;
                            for id = 1:fly_num
                                if status_filter(id,stim) == 1
                                    y_angv = abs(dtheta{id}(1, pre:post)');
                                    hold on
                                    plot(t,y_angv)
                                    xline(0,'r')
                                end
                            end
                        end
                        xlim([-tw(tp)/150 tw(tp)/150])
                        ylim([0 25])
                        %get fig image
                        I = getframe(gcf);
                        I3 = frame2im(I);
                    end
                    %concatenate trajectries/fv/angv figs
                    cd(save_dir)
                    primary_name = [expname '_trj_fv_angv_tw' num2str(tw(tp)) '.tif'];
                    I4 = [I1 I2 I3];
                    imwrite(I4, primary_name)
                end
                %save
                cd(save_dir)
                save(['fv_angv_tw' num2str(tw(tp)) '.mat'], 'fv_raw','fv_mean','angv_raw', 'angv_mean', '-mat')
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %leg tip spatial data analysis
            %data collection and visualization
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            posi_all = cell(1,6);%for 6legs
            for tp = 1:5
                for leg =1:6
                    %%% process for center of thorax
                    %alignment (position:start-point, angle:end-point)
                    posi_pre = cell(1,3); %column:w/m/s, level
                    posi_post = cell(1,3); %column:w/m/s
                    for level = 1:3 %3 level of stimuli
                        for n = 1:stn %stimuli times (3 or 10)
                            stim = n + stn*(level-1);
                            pre = indicatorLED.startframe(stim)-tw(tp);
                            on = indicatorLED.startframe(stim);
                            post = indicatorLED.startframe(stim)+tw(tp);
                            for id = 1:fly_num
                                if status_filter(id,stim) == 1
                                    %extract aligned values for pre/post
                                    %filter by stance start-end timing
                                    stc = [data{id,5}{2,leg};[data{id,5}{1,leg}(2:end) data_length]]; %stance timing row1(st),row2(ed)
                                    %%%%%
                                    %pre
                                    temp_st_idx = find(stc(1,:)>=pre & stc(1,:)<on);
                                    temp_st = stc(1,temp_st_idx);
                                    temp_ed = stc(2,temp_st_idx);
                                    if ~isempty(temp_ed)
                                        if temp_ed(end) > on
                                            temp_ed(end) = on;
                                        end
                                    end
                                    %add aligned data (pre and post)
                                    for p = 1:length(temp_st_idx)
                                        temp_pre{1} = data{id,7}{leg}(temp_st(p):temp_ed(p),:); %adjust starting position(pre)
                                        posi_pre{level} = [posi_pre{level} temp_pre];
                                    end
                                    %%%%%
                                    %post
                                    temp_st_idx = find(stc(1,:)>=on & stc(1,:)<post);
                                    temp_st = stc(1,temp_st_idx);
                                    temp_ed = stc(2,temp_st_idx);
                                    if ~isempty(temp_ed)
                                        if temp_ed(end) > post
                                            temp_ed(end) = post;
                                        end
                                    end
                                    %add aligned data (pre and post)
                                    for p = 1:length(temp_st_idx)
                                        temp_post{1} = data{id,7}{leg}(temp_st(p):temp_ed(p),:); %adjust starting position(pre)
                                        posi_post{level} = [posi_post{level} temp_post];
                                    end
                                end
                            end
                        end
                    end
                    posi_all{tp,leg} = [posi_pre;posi_post];
                end
                %save
                cd(save_dir)
                save(['posi_tw' num2str(tw(tp)) '.mat'], 'posi_all', '-mat')
            end
            %%%
            %visualization
            tp = 2; %plot figure for 30 frames trajectories
            leg_trj = figure('Position', [100, 100, 300, 700],'visible','off');
            for level = 1:3
                subplot(3,1,level)
                for leg = [1 4]
                    for wf = 1:length(posi_all{tp,leg}{1,level}) %overlay all walking fly
                        plot(posi_all{tp,leg}{1,level}{wf}(1:end,1),posi_all{tp,leg}{1,level}{wf}(1:end,2),'Color',col2(1,:))
                        hold on
                    end
                    for wf = 1:length(posi_all{tp,leg}{2,level}) %overlay all walking fly
                        plot(posi_all{tp,leg}{2,level}{wf}(1:end,1),posi_all{tp,leg}{2,level}{wf}(1:end,2),'Color',col2(2,:))
                        hold on
                    end
                end
                axis equal
                for leg = [2 5]
                    for wf = 1:length(posi_all{tp,leg}{1,level}) %overlay all walking fly
                        plot(posi_all{tp,leg}{1,level}{wf}(1:end,1),posi_all{tp,leg}{1,level}{wf}(1:end,2),'Color',col2(3,:))
                        hold on
                    end
                    for wf = 1:length(posi_all{tp,leg}{2,level}) %overlay all walking fly
                        plot(posi_all{tp,leg}{2,level}{wf}(1:end,1),posi_all{tp,leg}{2,level}{wf}(1:end,2),'Color',col2(4,:))
                        hold on
                    end
                end
                axis equal
                for leg = [3 6]
                    for wf = 1:length(posi_all{tp,leg}{1,level}) %overlay all walking fly
                        plot(posi_all{tp,leg}{1,level}{wf}(1:end,1),posi_all{tp,leg}{1,level}{wf}(1:end,2),'Color',col2(5,:))
                        hold on
                    end
                    for wf = 1:length(posi_all{tp,leg}{2,level}) %overlay all walking fly
                        plot(posi_all{tp,leg}{2,level}{wf}(1:end,1),posi_all{tp,leg}{2,level}{wf}(1:end,2),'Color',col2(6,:))
                        hold on
                    end
                end
                axis equal
            end
            %save fig
            cd(save_dir)
            leg_trj_name = [expname '_leg_trj_tw' num2str(tw(tp)) '.tif'];
            I = getframe(gcf);
            I = frame2im(I);
            imwrite(I, leg_trj_name)
            close all hidden
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % phase transition characteristic
            % stance(leg1) to swing(leg2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %pre-stim
            %swing onset latency distribution (during each stance phase)
            %from leg1 stance-onset(data{id,5}{1,leg1}), n to n+1
            %extract swing-onset delay for other 5 legs (leg2)(abs value)
            st2sw_ltc_pre = cell(stn*3,3);%row:stim#, column:3parameters(acc_ltc, acc_ltc_p, acc_st_dur)
            for stim = 1:stn*3
                st = indicatorLED.startframe(stim)-30;
                ed = indicatorLED.startframe(stim);
                ltc_dur_st = cell(fly_num,3);
                for id = 1:fly_num
                    temp_ltc = cell(6,6);
                    temp_ltc_p = cell(6,6);
                    temp_st_dur_all = cell(6,1);
                    if status_filter(id,stim) == 1
                        for leg1 = 1:6
                            temp_on = data{id,5}{2,leg1};
                            idx = find(temp_on >= st & temp_on <= ed);
                            for k = 1:length(idx)-1 %extract stance timing(leg1)
                                st_onset = data{id,5}{2,leg1}(idx(k));
                                next_st_onset = data{id,5}{2,leg1}(idx(k+1));
                                st_dur = next_st_onset - st_onset; %stance phase duration
                                for leg2 = 1:6
                                    if leg2 ~= leg1
                                        %find swing onset in the stance phase
                                        temp_sw_idx = find(data{id,5}{1,leg2} >= st_onset & data{id,5}{1,leg2} <= next_st_onset); %find idx in each leg1 stance phase
                                        %convert idx into onset frame and calculate delay
                                        temp_ltc{leg1,leg2} = [temp_ltc{leg1,leg2} data{id,5}{1,leg2}(temp_sw_idx) - st_onset];
                                        %convert absolute delay into proportional value
                                        temp_ltc_p{leg1,leg2} = [temp_ltc_p{leg1,leg2} (data{id,5}{1,leg2}(temp_sw_idx) - st_onset)/st_dur];
                                    end
                                end
                                temp_st_dur_all{leg1,1} = [temp_st_dur_all{leg1,1} st_dur];
                            end
                        end
                    end
                    ltc_dur_st{id,1} = temp_ltc;
                    ltc_dur_st{id,2} = temp_ltc_p;
                    ltc_dur_st{id,3} = temp_st_dur_all;
                end
                %accumulate ltc for all flies
                acc_ltc = cell(6,6);
                acc_ltc_p = cell(6,6);
                acc_st_dur = cell(6,1);
                for id = 1:fly_num
                    if status_filter(id,stim) == 1
                        for leg1 = 1:6
                            for leg2 = 1:6
                                acc_ltc{leg1,leg2} = [acc_ltc{leg1,leg2} ltc_dur_st{id,1}{leg1,leg2}];
                                acc_ltc_p{leg1,leg2} = [acc_ltc_p{leg1,leg2} ltc_dur_st{id,2}{leg1,leg2}];
                            end
                        end
                        acc_st_dur{leg1,1} = [acc_st_dur{leg1,1} ltc_dur_st{id,3}{leg1,1}];
                    end
                end
                %put acc data into all stim data (st2sw_ltc)
                st2sw_ltc_pre{stim,1} = acc_ltc;
                st2sw_ltc_pre{stim,2} = acc_ltc_p;
                st2sw_ltc_pre{stim,3} = acc_st_dur;
            end
            %accumulate ltc for all stimuli
            acc_ltc_pre = cell(6,6);
            for leg1 = 1:6
                for leg2 = 1:6
                    for stim =1:stn*3
                        acc_ltc_pre{leg1,leg2} = [acc_ltc_pre{leg1,leg2} st2sw_ltc_pre{stim,1}{leg1,leg2}];
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %stim-on
            %swing onset latency distribution (during each stance phase)
            %from leg1 stance-onset(data{id,5}{1,leg1}), n to n+1
            %extract swing-onset delay for other 5 legs (leg2)(abs value)
            for tp = 1:5
                st2sw_ltc_stim = cell(stn*3,3);%row:stim#, column:3parameters(acc_ltc, acc_ltc_p, acc_st_dur)
                for stim = 1:stn*3
                    st = indicatorLED.startframe(stim);
                    ed = indicatorLED.startframe(stim)+tw(tp);
                    ltc_dur_st = cell(fly_num,3);
                    for id = 1:fly_num
                        temp_ltc = cell(6,6);
                        temp_ltc_p = cell(6,6);
                        temp_st_dur_all = cell(6,1);
                        if status_filter(id,stim) == 1
                            for leg1 = 1:6
                                temp_on = data{id,5}{2,leg1};
                                idx = find(temp_on >= st & temp_on <= ed);
                                for k = 1:length(idx)-1 %extract stance timing(leg1)
                                    st_onset = data{id,5}{2,leg1}(idx(k));
                                    next_st_onset = data{id,5}{2,leg1}(idx(k+1));
                                    st_dur = next_st_onset - st_onset; %stance phase duration
                                    for leg2 = 1:6
                                        if leg2 ~= leg1
                                            %find swing onset in the stance phase
                                            temp_sw_idx = find(data{id,5}{1,leg2} >= st_onset & data{id,5}{1,leg2} <= next_st_onset); %find idx in each leg1 stance phase
                                            %convert idx into onset frame and calculate delay
                                            temp_ltc{leg1,leg2} = [temp_ltc{leg1,leg2} data{id,5}{1,leg2}(temp_sw_idx) - st_onset];
                                            %convert absolute delay into proportional value
                                            temp_ltc_p{leg1,leg2} = [temp_ltc_p{leg1,leg2} (data{id,5}{1,leg2}(temp_sw_idx) - st_onset)/st_dur];
                                        end
                                    end
                                    temp_st_dur_all{leg1,1} = [temp_st_dur_all{leg1,1} st_dur];
                                end
                            end
                        end
                        ltc_dur_st{id,1} = temp_ltc;
                        ltc_dur_st{id,2} = temp_ltc_p;
                        ltc_dur_st{id,3} = temp_st_dur_all;
                    end
                    %accumulate ltc for all flies
                    acc_ltc = cell(6,6);
                    acc_ltc_p = cell(6,6);
                    acc_st_dur = cell(6,1);
                    for leg1 = 1:6
                        for leg2 = 1:6
                            for id = 1:fly_num
                                acc_ltc{leg1,leg2} = [acc_ltc{leg1,leg2} ltc_dur_st{id,1}{leg1,leg2}];
                                acc_ltc_p{leg1,leg2} = [acc_ltc_p{leg1,leg2} ltc_dur_st{id,2}{leg1,leg2}];
                            end
                        end
                        acc_st_dur{leg1,1} = [acc_st_dur{leg1,1} ltc_dur_st{id,3}{leg1,1}];
                    end
                    %put acc data into all stim data (st2sw_ltc)
                    st2sw_ltc_stim{stim,1} = acc_ltc;
                    st2sw_ltc_stim{stim,2} = acc_ltc_p;
                    st2sw_ltc_stim{stim,3} = acc_st_dur;
                end
                %accumulate ltc for stim levels
                acc_ltc_weak = cell(6,6);
                acc_ltc_mid = cell(6,6);
                acc_ltc_str = cell(6,6);
                for level = 1:3
                    for leg1 = 1:6
                        for leg2 = 1:6
                            for stim =1+stn*(level-1):stn*level
                                if level == 1
                                    acc_ltc_weak{leg1,leg2} = [acc_ltc_weak{leg1,leg2} st2sw_ltc_stim{stim,1}{leg1,leg2}];
                                elseif level == 2
                                    acc_ltc_mid{leg1,leg2} = [acc_ltc_mid{leg1,leg2} st2sw_ltc_stim{stim,1}{leg1,leg2}];
                                elseif level == 3
                                    acc_ltc_str{leg1,leg2} = [acc_ltc_str{leg1,leg2} st2sw_ltc_stim{stim,1}{leg1,leg2}];
                                end
                            end
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                acc_ltc_4cd = cell(1,4);
                acc_ltc_4cd{1} = acc_ltc_pre;
                acc_ltc_4cd{2} = acc_ltc_weak;
                acc_ltc_4cd{3} = acc_ltc_mid;
                acc_ltc_4cd{4} = acc_ltc_str;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %visualize (post-stim)
                %plot latency distribution for 6legs(row) pre/w/m/s(column)
                %hist_p1 = figure('Position', [100, 100, 1400, 700],'visible','off');
                figure('Position', [100, 100, 1400, 700],'visible','off');
                for cnd = 1:4
                    for leg1 = 1:6
                        row_posi = [1 3 5 2 4 6];
                        subplot(6,4,4*(row_posi(leg1)-1) + cnd)%subplot for 6legs
                        for leg2 = 1:6
                            %if leg2~=leg1
                            temp_hist = histcounts(acc_ltc_4cd{cnd}{leg1,leg2},(0:1:30));
                            hold on
                            plot(temp_hist,'Color',color_idx(leg2,:),'LineWidth',2)
                            %end
                        end
                    end
                end
                %             %save fig (generate fig onlyfor tw 150)
                %             if tp == 4
                %                 cd(save_dir)
                %                 st2sw_hist_name = [expname '_st2sw_hist_tw' num2str(tw(tp)) '.tif'];
                %                 I = getframe(gcf);
                %                 imwrite(I.cdata, st2sw_hist_name)
                %             end
                %save variables
                cd(save_dir)
                save(['st2sw_ltc_tw' num2str(tw(tp)) '.mat'],'acc_ltc_4cd','-mat')
                close all hidden
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % phase transition characteristic
            % swing(leg1) to swing(leg2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %pre-stim
            %swing onset latency in swing phase
            %from leg1 swing-onset(data{id,5}{1,leg1}), n to n+1
            %extract swing-onset delay for other 5 legs (leg2)(abs value)
            sw2sw_ltc_pre = cell(stn*3,3);%row:stim#, column:3parameters(acc_ltc, acc_ltc_p, acc_sw_dur)
            for stim = 1:stn*3
                st = indicatorLED.startframe(stim) - 30;
                ed = indicatorLED.startframe(stim);
                ltc_dur_st = cell(fly_num,3);
                for id = 1:fly_num
                    temp_ltc = cell(6,6);
                    temp_ltc_p = cell(6,6);
                    temp_sw_dur_all = cell(6,1);
                    if status_filter(id,stim) == 1
                        for leg1 = 1:6
                            temp_on = data{id,5}{1,leg1};
                            idx = find(temp_on >= st & temp_on <= ed);
                            for k = 1:length(idx)-1 %extract stance timing(leg1)
                                sw_onset = data{id,5}{1,leg1}(idx(k));
                                next_sw_onset = data{id,5}{1,leg1}(idx(k+1));
                                sw_dur = next_sw_onset - sw_onset; %stance phase duration
                                for leg2 = 1:6
                                    if leg2 ~= leg1
                                        %find swing onset in the stance phase
                                        temp_sw_idx = find(data{id,5}{1,leg2} >= sw_onset & data{id,5}{1,leg2} <= next_sw_onset); %find idx in each leg1 stance phase
                                        %convert idx into onset frame and calculate delay
                                        temp_ltc{leg1,leg2} = [temp_ltc{leg1,leg2} data{id,5}{1,leg2}(temp_sw_idx) - sw_onset];
                                        %convert absolute delay into proportional angle
                                        temp_ltc_p{leg1,leg2} = [temp_ltc_p{leg1,leg2} pi*2*(data{id,5}{1,leg2}(temp_sw_idx) - sw_onset)/sw_dur];
                                    end
                                end
                                temp_sw_dur_all{leg1,1} = [temp_sw_dur_all{leg1,1} sw_dur];
                            end
                        end
                    end
                    ltc_dur_st{id,1} = temp_ltc;
                    ltc_dur_st{id,2} = temp_ltc_p;
                    ltc_dur_st{id,3} = temp_sw_dur_all;
                end
                %accumulate ltc for all flies
                acc_ltc = cell(6,6);
                acc_ltc_p = cell(6,6);
                acc_sw_dur = cell(6,1);
                for id = 1:fly_num
                    if status_filter(id,stim) == 1
                        for leg1 = 1:6
                            for leg2 = 1:6
                                acc_ltc{leg1,leg2} = [acc_ltc{leg1,leg2} ltc_dur_st{id,1}{leg1,leg2}];
                                acc_ltc_p{leg1,leg2} = [acc_ltc_p{leg1,leg2} ltc_dur_st{id,2}{leg1,leg2}];
                            end
                        end
                        acc_sw_dur{leg1,1} = [acc_sw_dur{leg1,1} ltc_dur_st{id,3}{leg1,1}];
                    end
                end
                %put acc data into all stim data (sw2sw_ltc)
                sw2sw_ltc_pre{stim,1} = acc_ltc;
                sw2sw_ltc_pre{stim,2} = acc_ltc_p;
                sw2sw_ltc_pre{stim,3} = acc_sw_dur;
            end
            %accumulate ltc for all stimuli
            acc_ltc_pre = cell(6,6);
            acc_ltc_p_pre = cell(6,6);
            for leg1 = 1:6
                for leg2 = 1:6
                    for stim =1:stn*3
                        acc_ltc_pre{leg1,leg2} = [acc_ltc_pre{leg1,leg2} sw2sw_ltc_pre{stim,1}{leg1,leg2}];
                        acc_ltc_p_pre{leg1,leg2} = [acc_ltc_p_pre{leg1,leg2} sw2sw_ltc_pre{stim,2}{leg1,leg2}];
                    end
                end
            end
            %%%
            %visualize
            %polarhistgram for neighboring pairs in 4 conditions (p/w/m/s)
            sw2sw_polar1 = figure('Position', [100, 100, 500, 700],'visible','off');
            %Intra pairs (F/M/H)
            for T = 1:3
                subplot(6,6,[3+12*(T-1),4+12*(T-1),9+12*(T-1),10+12*(T-1)])
                leg1 = T;
                leg2 = T+3;
                polarhistogram(acc_ltc_p_pre{leg1,leg2},30)
            end
            %Inter pairs
            %Left
            for T = 1:2
                subplot(6,6,[7+12*(T-1),8+12*(T-1),13+12*(T-1),14+12*(T-1)])
                leg1 = T;
                leg2 = T+1;
                polarhistogram(acc_ltc_p_pre{leg1,leg2},30)
            end
            %Right
            for T = 1:2
                subplot(6,6,[11+12*(T-1),12+12*(T-1),17+12*(T-1),18+12*(T-1)])
                leg1 = T+3;
                leg2 = T+4;
                polarhistogram(acc_ltc_p_pre{leg1,leg2},30)
            end
            %extract im data from fig
            sw2sw_polar_name = [expname '_polar_sw2sw_nb_p_tw30.tif'];
            Ip = getframe(gcf);
            Ip = frame2im(Ip);

            %polar histogram for tripod pairs (pre)
            sw2sw_polar2 = figure('Position', [100, 100, 500, 500],'visible','off');
            leg_pair = [1 5; 4 2; 2 6; 5 3];
            for k = 1:4
                subplot(2,2,k)
                polarhistogram(acc_ltc_p_pre{leg_pair(k,1),leg_pair(k,2)},30)
            end
            %save figs
            cd(save_dir)
            sw2sw_polar_name = [expname '_polar_sw2sw_tri_p_tw30.tif'];
            I = getframe(gcf);
            imwrite(I.cdata, sw2sw_polar_name)
            close all hidden
            for tp = 1:5
                I_nb = Ip;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %stim-on
                %swing onset latency distribution (during each stance phase)
                %from leg1 stance-onset(data{id,5}{1,leg1}), n to n+1
                %extract swing-onset delay for other 5 legs (leg2)(abs value)
                sw2sw_ltc_stim = cell(stn*3,3);%row:stim#, column:3parameters(acc_ltc, acc_ltc_p, acc_sw_dur)
                for stim = 1:stn*3
                    st = indicatorLED.startframe(stim);
                    ed = indicatorLED.startframe(stim)+tw(tp);
                    ltc_dur_st = cell(fly_num,3);
                    for id = 1:fly_num
                        temp_ltc = cell(6,6);
                        temp_ltc_p = cell(6,6);
                        temp_sw_dur_all = cell(6,1);
                        if status_filter(id,stim) == 1
                            for leg1 = 1:6
                                temp_on = data{id,5}{1,leg1};
                                idx = find(temp_on >= st & temp_on <= ed);
                                for k = 1:length(idx)-1 %extract stance timing(leg1)
                                    sw_onset = data{id,5}{1,leg1}(idx(k));
                                    next_sw_onset = data{id,5}{1,leg1}(idx(k+1));
                                    sw_dur = next_sw_onset - sw_onset; %stance phase duration
                                    for leg2 = 1:6
                                        if leg2 ~= leg1
                                            %find swing onset in the stance phase
                                            temp_sw_idx = find(data{id,5}{1,leg2} >= sw_onset & data{id,5}{1,leg2} <= next_sw_onset); %find idx in each leg1 stance phase
                                            %convert idx into onset frame and calculate delay
                                            temp_ltc{leg1,leg2} = [temp_ltc{leg1,leg2} data{id,5}{1,leg2}(temp_sw_idx) - sw_onset];
                                            %convert absolute delay into proportional angle
                                            temp_ltc_p{leg1,leg2} = [temp_ltc_p{leg1,leg2} pi*2*(data{id,5}{1,leg2}(temp_sw_idx) - sw_onset)/sw_dur];
                                        end
                                    end
                                    temp_sw_dur_all{leg1,1} = [temp_sw_dur_all{leg1,1} sw_dur];
                                end
                            end
                        end
                        ltc_dur_st{id,1} = temp_ltc;
                        ltc_dur_st{id,2} = temp_ltc_p;
                        ltc_dur_st{id,3} = temp_sw_dur_all;
                    end
                    %accumulate ltc for all flies
                    acc_ltc = cell(6,6);
                    acc_ltc_p = cell(6,6);
                    acc_sw_dur = cell(6,1);
                    for leg1 = 1:6
                        for leg2 = 1:6
                            for id = 1:fly_num
                                acc_ltc{leg1,leg2} = [acc_ltc{leg1,leg2} ltc_dur_st{id,1}{leg1,leg2}];
                                acc_ltc_p{leg1,leg2} = [acc_ltc_p{leg1,leg2} ltc_dur_st{id,2}{leg1,leg2}];
                            end
                        end
                        acc_sw_dur{leg1,1} = [acc_sw_dur{leg1,1} ltc_dur_st{id,3}{leg1,1}];
                    end
                    %put acc data into all stim data (sw2sw_ltc)
                    sw2sw_ltc_stim{stim,1} = acc_ltc;
                    sw2sw_ltc_stim{stim,2} = acc_ltc_p;
                    sw2sw_ltc_stim{stim,3} = acc_sw_dur;
                end
                %accumulate ltc for stim levels
                acc_ltc_weak = cell(6,6);
                acc_ltc_mid = cell(6,6);
                acc_ltc_str = cell(6,6);
                acc_ltc_p_weak = cell(6,6);
                acc_ltc_p_mid = cell(6,6);
                acc_ltc_p_str = cell(6,6);
                for level = 1:3
                    for leg1 = 1:6
                        for leg2 = 1:6
                            for stim =1+stn*(level-1):stn*level
                                if level == 1
                                    acc_ltc_weak{leg1,leg2} = [acc_ltc_weak{leg1,leg2} sw2sw_ltc_stim{stim,1}{leg1,leg2}];
                                    acc_ltc_p_weak{leg1,leg2} = [acc_ltc_p_weak{leg1,leg2} sw2sw_ltc_stim{stim,2}{leg1,leg2}];
                                elseif level == 2
                                    acc_ltc_mid{leg1,leg2} = [acc_ltc_mid{leg1,leg2} sw2sw_ltc_stim{stim,1}{leg1,leg2}];
                                    acc_ltc_p_mid{leg1,leg2} = [acc_ltc_p_mid{leg1,leg2} sw2sw_ltc_stim{stim,2}{leg1,leg2}];
                                elseif level == 3
                                    acc_ltc_str{leg1,leg2} = [acc_ltc_str{leg1,leg2} sw2sw_ltc_stim{stim,1}{leg1,leg2}];
                                    acc_ltc_p_str{leg1,leg2} = [acc_ltc_p_str{leg1,leg2} sw2sw_ltc_stim{stim,2}{leg1,leg2}];
                                end
                            end
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                acc_ltc_4cd = cell(1,4);
                acc_ltc_4cd{1} = acc_ltc_pre;
                acc_ltc_4cd{2} = acc_ltc_weak;
                acc_ltc_4cd{3} = acc_ltc_mid;
                acc_ltc_4cd{4} = acc_ltc_str;
                acc_ltc_p_4cd = cell(1,4);
                acc_ltc_p_4cd{1} = acc_ltc_p_pre;
                acc_ltc_p_4cd{2} = acc_ltc_p_weak;
                acc_ltc_p_4cd{3} = acc_ltc_p_mid;
                acc_ltc_p_4cd{4} = acc_ltc_p_str;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %visualize
                %polarhistgram for neighboring pairs in 4 conditions (p/w/m/s)
                if tp == 2 || tp == 4
                    for cnd = 2:4
                        sw2sw_polar1 = figure('Position', [100, 100, 500, 700],'visible','off');
                        %Intra pairs (F/M/H)
                        for T = 1:3
                            subplot(6,6,[3+12*(T-1),4+12*(T-1),9+12*(T-1),10+12*(T-1)])
                            leg1 = T;
                            leg2 = T+3;
                            polarhistogram(acc_ltc_p_4cd{cnd}{leg1,leg2},30)
                        end
                        %Inter pairs
                        %Left
                        for T = 1:2
                            subplot(6,6,[7+12*(T-1),8+12*(T-1),13+12*(T-1),14+12*(T-1)])
                            leg1 = T;
                            leg2 = T+1;
                            polarhistogram(acc_ltc_p_4cd{cnd}{leg1,leg2},30)
                        end
                        %Right
                        for T = 1:2
                            subplot(6,6,[11+12*(T-1),12+12*(T-1),17+12*(T-1),18+12*(T-1)])
                            leg1 = T+3;
                            leg2 = T+4;
                            polarhistogram(acc_ltc_p_4cd{cnd}{leg1,leg2},30)
                        end
                        %extract imdata from fig
                        %sw2sw_polar_name = [expname '_polar_sw2sw_nb_' cnd_name(cnd) '_tw' num2str(tw(tp)) '.tif'];
                        I = getframe(gcf);
                        I = frame2im(I);
                        I_nb = [I_nb I];
                        %imwrite(I.cdata, sw2sw_polar_name)
                    end
                    %save fig
                    cd(save_dir)
                    sw2sw_polar_name = [expname '_polar_sw2sw_nb_tw' num2str(tw(tp)) '.tif'];
                    imwrite(I_nb, sw2sw_polar_name)
                end
                %save variables
                save(['sw2sw_p_ltc_tw' num2str(tw(tp)) '.mat'],'acc_ltc_p_4cd','-mat')
                close all hidden
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % phase transition characteristic
            % stance(leg1) to stance(leg2)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %pre-stim
            %stance onset latency in swing phase
            %from leg1 stance-onset(data{id,5}{2,leg1}), n to n+1
            %extract stance-onset delay for other 5 legs (leg2)(abs value)
            st2st_ltc_pre = cell(stn*3,3);%row:stim#, column:3parameters(acc_ltc, acc_ltc_p, acc_st_dur)
            for stim = 1:stn*3
                st = indicatorLED.startframe(stim)-30;
                ed = indicatorLED.startframe(stim);
                ltc_dur_st = cell(fly_num,3);
                for id = 1:fly_num
                    temp_ltc = cell(6,6);
                    temp_ltc_p = cell(6,6);
                    temp_st_dur_all = cell(6,1);
                    if status_filter(id,stim) == 1
                        for leg1 = 1:6
                            temp_on = data{id,5}{2,leg1};
                            idx = find(temp_on >= st & temp_on <= ed);
                            for k = 1:length(idx)-1 %extract stance timing(leg1)
                                st_onset = data{id,5}{2,leg1}(idx(k));
                                next_st_onset = data{id,5}{2,leg1}(idx(k+1));
                                st_dur = next_st_onset - st_onset; %stance phase duration
                                for leg2 = 1:6
                                    if leg2 ~= leg1
                                        %find swing onset in the stance phase
                                        temp_st_idx = find(data{id,5}{2,leg2} >= st_onset & data{id,5}{2,leg2} <= next_st_onset); %find idx in each leg1 stance phase
                                        %convert idx into onset frame and calculate delay
                                        temp_ltc{leg1,leg2} = [temp_ltc{leg1,leg2} data{id,5}{2,leg2}(temp_st_idx) - st_onset];
                                        %convert absolute delay into proportional angle
                                        temp_ltc_p{leg1,leg2} = [temp_ltc_p{leg1,leg2} pi*2*(data{id,5}{2,leg2}(temp_st_idx) - st_onset)/st_dur];
                                    end
                                end
                                temp_st_dur_all{leg1,1} = [temp_st_dur_all{leg1,1} st_dur];
                            end
                        end
                    end
                    ltc_dur_st{id,1} = temp_ltc;
                    ltc_dur_st{id,2} = temp_ltc_p;
                    ltc_dur_st{id,3} = temp_st_dur_all;
                end
                %accumulate ltc for all flies
                acc_ltc = cell(6,6);
                acc_ltc_p = cell(6,6);
                acc_st_dur = cell(6,1);
                for id = 1:fly_num
                    if status_filter(id,stim) == 1
                        for leg1 = 1:6
                            for leg2 = 1:6
                                acc_ltc{leg1,leg2} = [acc_ltc{leg1,leg2} ltc_dur_st{id,1}{leg1,leg2}];
                                acc_ltc_p{leg1,leg2} = [acc_ltc_p{leg1,leg2} ltc_dur_st{id,2}{leg1,leg2}];
                            end
                        end
                        acc_st_dur{leg1,1} = [acc_st_dur{leg1,1} ltc_dur_st{id,3}{leg1,1}];
                    end
                end
                %put acc data into all stim data (st2st_ltc)
                st2st_ltc_pre{stim,1} = acc_ltc;
                st2st_ltc_pre{stim,2} = acc_ltc_p;
                st2st_ltc_pre{stim,3} = acc_st_dur;
            end
            %accumulate ltc for all stimuli
            acc_ltc_pre = cell(6,6);
            acc_ltc_p_pre = cell(6,6);
            for leg1 = 1:6
                for leg2 = 1:6
                    for stim =1:stn*3
                        acc_ltc_pre{leg1,leg2} = [acc_ltc_pre{leg1,leg2} st2st_ltc_pre{stim,1}{leg1,leg2}];
                        acc_ltc_p_pre{leg1,leg2} = [acc_ltc_p_pre{leg1,leg2} st2st_ltc_pre{stim,2}{leg1,leg2}];
                    end
                end
            end
            %         %%%%%%
            %         %visualize
            %         %polarhistgram for neighboring pairs
            %         st2st_polar1 = figure('Position', [100, 100, 500, 700],'visible','off');
            %         %Intra pairs (F/M/H)
            %         for T = 1:3
            %             subplot(6,6,[3+12*(T-1),4+12*(T-1),9+12*(T-1),10+12*(T-1)])
            %             leg1 = T;
            %             leg2 = T+3;
            %             polarhistogram(acc_ltc_p_4cd{cnd}{leg1,leg2},30)
            %         end
            %         %Inter pairs
            %         %Left
            %         for T = 1:2
            %             subplot(6,6,[7+12*(T-1),8+12*(T-1),13+12*(T-1),14+12*(T-1)])
            %             leg1 = T;
            %             leg2 = T+1;
            %             polarhistogram(acc_ltc_p_4cd{cnd}{leg1,leg2},30)
            %         end
            %         %Right
            %         for T = 1:2
            %             subplot(6,6,[11+12*(T-1),12+12*(T-1),17+12*(T-1),18+12*(T-1)])
            %             leg1 = T+3;
            %             leg2 = T+4;
            %             polarhistogram(acc_ltc_p_4cd{cnd}{leg1,leg2},30)
            %         end
            %         %extract im data from fig
            %         st2st_polar_name = [expname '_polar_st2st_nb_p_tw30.tif'];
            %         Ip = getframe(gcf);
            %         Ip = frame2im(Ip);
            %
            %         %polar histogram for tripod pairs (pre)
            %         st2st_polar2 = figure('Position', [100, 100, 500, 500],'visible','off');
            %         leg_pair = [1 5; 4 2; 2 6; 5 3];
            %         for k = 1:4
            %             subplot(2,2,k)
            %             polarhistogram(acc_ltc_p_4cd{cnd}{leg_pair(k,1),leg_pair(k,2)},30)
            %         end
            %         %save figs
            %         cd(save_dir)
            %         st2st_polar_name = [expname '_polar_st2st_tri_p_tw30.tif'];
            %         I = getframe(gcf);
            %         imwrite(I.cdata,st2st_polar_name)
            %         close all hidden
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %stim-on
            %swing onset latency distribution (during each stance phase)
            %from leg1 stance-onset(data{id,5}{2,leg1}), n to n+1
            %extract swing-onset delay for other 5 legs (leg2)(abs value)
            for tp = 1:5
                %I_nb = Ip;
                st2st_ltc_stim = cell(stn*3,3);%row:stim#, column:3parameters(acc_ltc, acc_ltc_p, acc_st_dur)
                for stim = 1:stn*3
                    st = indicatorLED.startframe(stim);
                    ed = indicatorLED.startframe(stim)+tw(tp);
                    ltc_dur_st = cell(fly_num,3);
                    for id = 1:fly_num
                        temp_ltc = cell(6,6);
                        temp_ltc_p = cell(6,6);
                        temp_st_dur_all = cell(6,1);
                        if status_filter(id,stim) == 1
                            for leg1 = 1:6
                                temp_on = data{id,5}{2,leg1};
                                idx = find(temp_on >= st & temp_on <= ed);
                                for k = 1:length(idx)-1 %extract stance timing(leg1)
                                    st_onset = data{id,5}{2,leg1}(idx(k));
                                    next_st_onset = data{id,5}{2,leg1}(idx(k+1));
                                    st_dur = next_st_onset - st_onset; %stance phase duration
                                    for leg2 = 1:6
                                        if leg2 ~= leg1
                                            %find swing onset in the stance phase
                                            temp_st_idx = find(data{id,5}{2,leg2} >= st_onset & data{id,5}{2,leg2} <= next_st_onset); %find idx in each leg1 stance phase
                                            %convert idx into onset frame and calculate delay
                                            temp_ltc{leg1,leg2} = [temp_ltc{leg1,leg2} data{id,5}{2,leg2}(temp_st_idx) - st_onset];
                                            %convert absolute delay into proportional angle
                                            temp_ltc_p{leg1,leg2} = [temp_ltc_p{leg1,leg2} pi*2*(data{id,5}{2,leg2}(temp_st_idx) - st_onset)/st_dur];
                                        end
                                    end
                                    temp_st_dur_all{leg1,1} = [temp_st_dur_all{leg1,1} st_dur];
                                end
                            end
                        end
                        ltc_dur_st{id,1} = temp_ltc;
                        ltc_dur_st{id,2} = temp_ltc_p;
                        ltc_dur_st{id,3} = temp_st_dur_all;
                    end
                    %accumulate ltc for all flies
                    acc_ltc = cell(6,6);
                    acc_ltc_p = cell(6,6);
                    acc_st_dur = cell(6,1);
                    for leg1 = 1:6
                        for leg2 = 1:6
                            for id = 1:fly_num
                                acc_ltc{leg1,leg2} = [acc_ltc{leg1,leg2} ltc_dur_st{id,1}{leg1,leg2}];
                                acc_ltc_p{leg1,leg2} = [acc_ltc_p{leg1,leg2} ltc_dur_st{id,2}{leg1,leg2}];
                            end
                        end
                        acc_st_dur{leg1,1} = [acc_st_dur{leg1,1} ltc_dur_st{id,3}{leg1,1}];
                    end
                    %put acc data into all stim data (st2st_ltc)
                    st2st_ltc_stim{stim,1} = acc_ltc;
                    st2st_ltc_stim{stim,2} = acc_ltc_p;
                    st2st_ltc_stim{stim,3} = acc_st_dur;
                end
                %accumulate ltc for stim levels
                acc_ltc_weak = cell(6,6);
                acc_ltc_mid = cell(6,6);
                acc_ltc_str = cell(6,6);
                acc_ltc_p_weak = cell(6,6);
                acc_ltc_p_mid = cell(6,6);
                acc_ltc_p_str = cell(6,6);
                for level = 1:3
                    for leg1 = 1:6
                        for leg2 = 1:6
                            for stim =1 + stn*(level-1) : stn*level
                                if level == 1
                                    acc_ltc_weak{leg1,leg2} = [acc_ltc_weak{leg1,leg2} st2st_ltc_stim{stim,1}{leg1,leg2}];
                                    acc_ltc_p_weak{leg1,leg2} = [acc_ltc_p_weak{leg1,leg2} st2st_ltc_stim{stim,2}{leg1,leg2}];
                                elseif level == 2
                                    acc_ltc_mid{leg1,leg2} = [acc_ltc_mid{leg1,leg2} st2st_ltc_stim{stim,1}{leg1,leg2}];
                                    acc_ltc_p_mid{leg1,leg2} = [acc_ltc_p_mid{leg1,leg2} st2st_ltc_stim{stim,2}{leg1,leg2}];
                                elseif level == 3
                                    acc_ltc_str{leg1,leg2} = [acc_ltc_str{leg1,leg2} st2st_ltc_stim{stim,1}{leg1,leg2}];
                                    acc_ltc_p_str{leg1,leg2} = [acc_ltc_p_str{leg1,leg2} st2st_ltc_stim{stim,2}{leg1,leg2}];
                                end
                            end
                        end
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                acc_ltc_4cd = cell(1,4);
                acc_ltc_4cd{1} = acc_ltc_pre;
                acc_ltc_4cd{2} = acc_ltc_weak;
                acc_ltc_4cd{3} = acc_ltc_mid;
                acc_ltc_4cd{4} = acc_ltc_str;
                acc_ltc_p_4cd = cell(1,4);
                acc_ltc_p_4cd{1} = acc_ltc_p_pre;
                acc_ltc_p_4cd{2} = acc_ltc_p_weak;
                acc_ltc_p_4cd{3} = acc_ltc_p_mid;
                acc_ltc_p_4cd{4} = acc_ltc_p_str;
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %             %visualize
                %             %polarhistgram for neighboring pairs in 4 conditions (p/w/m/s)
                %             if tp == 2 || tp == 4
                %                 for cnd = 2:4
                %                     st2st_polar1 = figure('Position', [100, 100, 500, 700],'visible','off');
                %                     %Intra pairs (F/M/H)
                %                     for T = 1:3
                %                         subplot(6,6,[3+12*(T-1),4+12*(T-1),9+12*(T-1),10+12*(T-1)])
                %                         leg1 = T;
                %                         leg2 = T+3;
                %                         polarhistogram(acc_ltc_p_4cd{cnd}{leg1,leg2},30)
                %                     end
                %                     %Inter pairs
                %                     %Left
                %                     for T = 1:2
                %                         subplot(6,6,[7+12*(T-1),8+12*(T-1),13+12*(T-1),14+12*(T-1)])
                %                         leg1 = T;
                %                         leg2 = T+1;
                %                         polarhistogram(acc_ltc_p_4cd{cnd}{leg1,leg2},30)
                %                     end
                %                     %Right
                %                     for T = 1:2
                %                         subplot(6,6,[11+12*(T-1),12+12*(T-1),17+12*(T-1),18+12*(T-1)])
                %                         leg1 = T+3;
                %                         leg2 = T+4;
                %                         polarhistogram(acc_ltc_p_4cd{cnd}{leg1,leg2},30)
                %                     end
                %                     %extract imdata from fig
                %                     st2st_polar_name = [expname '_polar_st2st_nb_' cnd_name(cnd) '_tw' num2str(tw(tp)) '.tif'];
                %                     I = getframe(gcf);
                %                     I = frame2im(I);
                %                     I_nb = [I_nb I];
                %                 end
                %                 %save fig
                %                 cd(save_dir)
                %                 st2st_polar_name = [expname '_polar_st2st_nb_tw' num2str(tw(tp)) '.tif'];
                %                 imwrite(I_nb, st2st_polar_name)
                %             end
                %             %save variables
                %             save(['st2st_p_ltc_tw' num2str(tw(tp)) '.mat'],'acc_ltc_p_4cd','-mat')
                %             close all hidden
            end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %scatter plot for combination of leg tip velocities
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %data preparation
            lt_pre = cell(fly_num,stn*3);
            lt_post = cell(fly_num,stn*3);
            for stim = 1 : stn*3
                for id = 1:fly_num
                    lt_pre{id,stim} = cell(1,6);
                    lt_post{id,stim} = cell(1,6);
                end
            end
            %assign values
            for leg = 1:6 %L1-3:1-3, R1-3:4-6
                for stim = 1 : stn*3
                    st = indicatorLED.startframe(stim);
                    for id = 1:fly_num
                        lt_pre{id,stim}{leg} = data{id,2}{leg}(st-30:st);
                        lt_post{id,stim}{leg} = data{id,2}{leg}(st:st+30);
                    end
                end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %setting of graph
            sz = 6;%size of scatter plot
            ax_max = 80;
            out_th = [8 8 8]; %out of L region threshold (T1/2/3)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %intra-seg scatter plot
            %classify and count values in 3 category
            intra_classi = cell(3,4); %prepare cell to put category counts, T1-3 x pre/w/m/s
            %pre_stim
            %figure('Position', [100, 100, 1000, 1000],'visible','off');
            for T = 1:3
                cat1 = 0; %around origin of coordinates
                cat2_1 = 0; %along the axis
                cat2_2 = 0; %along the axis
                cat3 = 0; %out of L region
                subplot(3,4,1+4*(T-1))
                axis([0 ax_max 0 ax_max])
                hold on
                for stim = stn*1+1 : stn*2
                    for id = 1:fly_num
                        if status_filter(id,stim) == 1
                            hold on
                            %scatter(lt_pre{id,stim}{T},lt_pre{id,stim}{T+3},sz,'filled','k');
                            for k = 1:31
                                %classify and count values
                                if lt_pre{id,stim}{T}(k) <= out_th(T) && lt_pre{id,stim}{T+3}(k) <= out_th(T)
                                    cat1 = cat1 +1;
                                elseif lt_pre{id,stim}{T}(k) > out_th(T) && lt_pre{id,stim}{T+3}(k) > out_th(T)
                                    cat3 = cat3 +1;
                                elseif lt_pre{id,stim}{T}(k) <= out_th(T) && lt_pre{id,stim}{T+3}(k) > out_th(T)
                                    cat2_1 = cat2_1 + 1;
                                elseif lt_pre{id,stim}{T}(k) > out_th(T) && lt_pre{id,stim}{T+3}(k) <= out_th(T)
                                    cat2_2 = cat2_2 + 1;
                                end
                            end
                        end
                    end
                end
                classi = [cat2_1 cat3; cat1 cat2_2];
                classi_rel = classi/(cat1+cat2_1+cat2_2+cat3);
                classi = [classi classi_rel];
                intra_classi{T,1} = classi;
                %add referential line in scatter plot
                % xline(out_th(T),'--r')
                % yline(out_th(T),'--r')
            end
            %weak/middle/strong
            for level = 1:3
                for T = 1:3
                    cat1 = 0; %around origin of coordinates
                    cat2_1 = 0; %along the axis
                    cat2_2 = 0; %along the axis
                    cat3 = 0; %out of L region\
                    subplot(3,4,level+1+(T-1)*4)
                    axis([0 ax_max 0 ax_max])
                    hold on
                    for stim = 1+(level-1)*stn:level*stn
                        for id = 1:fly_num
                            if status_filter(id,stim) == 1
                                hold on
                                %scatter(lt_post{id,stim}{T},lt_post{id,stim}{T+3},sz,'filled','k');
                                for k = 1:31
                                    %classify and count values
                                    if lt_post{id,stim}{T}(k) <= out_th(T) && lt_post{id,stim}{T+3}(k) <= out_th(T)
                                        cat1 = cat1 +1;
                                    elseif lt_post{id,stim}{T}(k) > out_th(T) && lt_post{id,stim}{T+3}(k) > out_th(T)
                                        cat3 = cat3 +1;
                                    elseif lt_post{id,stim}{T}(k) <= out_th(T) && lt_post{id,stim}{T+3}(k) > out_th(T)
                                        cat2_1 = cat2_1 + 1;
                                    elseif lt_post{id,stim}{T}(k) > out_th(T) && lt_post{id,stim}{T+3}(k) <= out_th(T)
                                        cat2_2 = cat2_2 + 1;
                                    end
                                end
                            end
                        end
                    end
                    classi = [cat2_1 cat3; cat1 cat2_2];
                    classi_rel = classi/(cat1+cat2_1+cat2_2+cat3);
                    classi = [classi classi_rel];
                    intra_classi{T,level+1} = classi;
                    %xline(out_th(T),'--r')
                    %yline(out_th(T),'--r')
                end
            end
            %         %save
            %         cd(save_dir)
            %         fig_sc_name = [expname '_fig_sc_intraseg.tif'];
            %         I = getframe(gcf);
            %         imwrite(I.cdata, fig_sc_name)
            %         close all hidden

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %intersegmental combination
            %accumulate data for all flies
            inter_classi = cell(4,4); %prepare cell to put category counts
            %pre_stim
            %         fig_sc = figure('Position', [100, 100, 1000, 1000],'visible','off');
            for LR = 1:2
                for p = 1:2 %pro-meso/meso-meta
                    cat1 = 0; %around origin of coordinates
                    cat2_1 = 0; %along the axis
                    cat2_2 = 0; %along the axis
                    cat3 = 0; %out of L region
                    subplot(4,4,1+(LR-1)*4+(p-1)*8)
                    axis([0 ax_max 0 ax_max])
                    hold on
                    for stim = stn*1+1 : stn*2
                        for id = 1:fly_num
                            if status_filter(id,stim) == 1
                                hold on
                                %scatter(lt_pre{id,stim}{(LR-1)*3+p},lt_pre{id,stim}{(LR-1)*3+p+1},sz,'filled','k');
                                for k = 1:31
                                    %classify and count values
                                    if lt_pre{id,stim}{(LR-1)*3+p}(k) <= out_th(p) && lt_pre{id,stim}{(LR-1)*3+p+1}(k) <= out_th(p)
                                        cat1 = cat1 +1;
                                    elseif lt_pre{id,stim}{(LR-1)*3+p}(k) > out_th(p) && lt_pre{id,stim}{(LR-1)*3+p+1}(k) > out_th(p)
                                        cat3 = cat3 +1;
                                    elseif lt_pre{id,stim}{(LR-1)*3+p}(k) <= out_th(p) && lt_pre{id,stim}{(LR-1)*3+p+1}(k) > out_th(p)
                                        cat2_1 = cat2_1+1;
                                    elseif lt_pre{id,stim}{(LR-1)*3+p}(k) > out_th(p) && lt_pre{id,stim}{(LR-1)*3+p+1}(k) <= out_th(p)
                                        cat2_2 = cat2_2+1;
                                    end
                                end
                            end
                        end
                    end
                    classi = [cat2_1 cat3; cat1 cat2_2];
                    classi_rel = classi/(cat1+cat2_1+cat2_2+cat3);
                    classi = [classi classi_rel];
                    inter_classi{LR+(p-1)*2,1} = classi;
                    %xline(out_th(p),'--r')
                    %yline(out_th(p),'--r')
                end
            end
            %weak/middle/strong
            for level = 1:3
                for LR = 1:2
                    for p = 1:2
                        cat1 = 0; %around origin of coordinates
                        cat2_1 = 0; %along the axis
                        cat2_2 = 0; %along the axis
                        cat3 = 0; %out of L region
                        subplot(4,4,level+1+(LR-1)*4+(p-1)*8)
                        axis([0 ax_max 0 ax_max])
                        hold on
                        for stim = 1+(level-1)*stn:level*stn
                            for id = 1:fly_num
                                if status_filter(id,stim) == 1
                                    hold on
                                    %scatter(lt_post{id,stim}{(LR-1)*3+p},lt_post{id,stim}{(LR-1)*3+p+1},sz,'filled','k');
                                    for k = 1:31
                                        %classify and count values
                                        if lt_post{id,stim}{(LR-1)*3+p}(k) <= out_th(p) && lt_post{id,stim}{(LR-1)*3+p+1}(k) <= out_th(p)
                                            cat1 = cat1 +1;
                                        elseif lt_post{id,stim}{(LR-1)*3+p}(k) > out_th(p) && lt_post{id,stim}{(LR-1)*3+p+1}(k) > out_th(p)
                                            cat3 = cat3 +1;
                                        elseif lt_post{id,stim}{(LR-1)*3+p}(k) <= out_th(p) && lt_post{id,stim}{(LR-1)*3+p+1}(k) > out_th(p)
                                            cat2_1 = cat2_1+1;
                                        elseif lt_post{id,stim}{(LR-1)*3+p}(k) > out_th(p) && lt_post{id,stim}{(LR-1)*3+p+1}(k) <= out_th(p)
                                            cat2_2 = cat2_2+1;
                                        end
                                    end
                                end
                            end
                        end
                        classi = [cat2_1 cat3; cat1 cat2_2];
                        classi_rel = classi/(cat1+cat2_1+cat2_2+cat3);
                        classi = [classi classi_rel];
                        inter_classi{LR+(p-1)*2,level+1} = classi;
                        %xline(out_th(p),'--r')
                        %yline(out_th(p),'--r')
                    end
                end
            end
            %         %save
            %         cd(save_dir)
            %         fig_sc_name = [expname '_fig_sc_interseg.tif'];
            %         I = getframe(gcf);
            %         imwrite(I.cdata, fig_sc_name)
            %         close all hidden
            %save intra/inter-classi
            save('classi.mat','intra_classi','inter_classi')
            %adjust leg_tip_v by using status filter
            for stim = 1:stn*3
                for id = 1:fly_num
                    if status_filter(id,stim) == 0
                        lt_pre{id,stim} = [];
                        lt_post{id,stim} = [];
                    end
                end
            end
            save('leg_tip_v.mat','lt_pre','lt_post','-mat')
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            close all hidden
        else
            cd(error_dir)
            mkdir(APT_list{d,1})
            cd(APT_list{d,1})
            error_dir2 = pwd;
            cd(analysis_dir)
            copyfile(APT_list{d,1},error_dir2)
            rmdir(APT_list{d,1},'s')
        end
    end
    cd(analysis_dir)
    toc
end


