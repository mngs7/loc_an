clear;
warning off;

%time window for analysis
tw = [10; 30; 50; 150; 300];%different time windows
%set analysis_dir
analysis_dir = uigetdir();
cd(analysis_dir)
if isfolder('acc_results')
    rmdir('acc_results', 's')
end
mkdir acc_results
cd acc_results
acc_dir = pwd;

%scan APT exp list (will also check if the files for analysis exist)
cd(analysis_dir)
APT_list = struct2cell(dir('VNC*'))';
%extract line name
line_names = [];
for k = 1:length(APT_list(:,1))
    temp = APT_list{k,1};
    if strfind(temp,'VNC_') > 0
        line = extractBetween(temp,5,15);
        line_names = [line_names; line];
    elseif strfind(temp,'VNC2_') > 0
        line = extractBetween(temp,6,16);
        line_names = [line_names; line];
    end
end
uni_lines = unique(line_names,'rows');

%check protocol type (3x3stim or 10x3stim)
%not for mixed protocol (work only either type)
if strfind(APT_list{1,1},'VNC_') > 0
    stn = 3; %2021 protocol
elseif strfind(APT_list{1,1},'VNC2_') > 0
    stn = 10; %2022 protocol
end

%process each ss line
cd(analysis_dir)
tp = 2; %tp = 1:5
for ss = 1:length(uni_lines(:,1))
    ss_name = uni_lines{ss,1};
    %create ss name folder
    cd(acc_dir)
    mkdir(ss_name)
    cd(ss_name)
    save_dir = pwd;
    cd(analysis_dir)
    %cumulative fv angv
    fv_mean_all = cell(2,3);
    angv_mean_all = cell(2,3);
    %load thx position
    %acc_posi_thx

    %cumulative leg-tip position

    %cumulative leg tip velocity
    lt_pre_all = [];
    lt_post_all = [];
    %cumulative phase latency
    acc_ltc_4cd_all_st2sw = cell(1,4);%4 conditions (p/w/m/s)
    acc_ltc_4cd_all_st2sw{1} = cell(6,6);%leg 6-by-6
    acc_ltc_4cd_all_st2sw{2} = cell(6,6);
    acc_ltc_4cd_all_st2sw{3} = cell(6,6);
    acc_ltc_4cd_all_st2sw{4} = cell(6,6);
    acc_ltc_p_4cd_all_sw2sw = acc_ltc_4cd_all_st2sw;%sw2sw copied from st2sw
    %cumulative classification
    intra_classi_all = [];
    inter_classi_all = [];
    %find experiment folders for the ss line
    for d = 1:length(APT_list(:,1))
        if strfind(APT_list{d,1},ss_name) > 0
            cd(analysis_dir)
            cd(APT_list{d,1})
            cd results
            %accumulate fv angv 
            load(['fv_angv_tw' num2str(tw(tp)) '.mat'])
            for row = 1:2
                for col = 1:3
                    fv_mean_all{row,col} = [fv_mean_all{row,col} fv_mean{row,col}];
                    angv_mean_all{row,col} = [angv_mean_all{row,col} angv_mean{row,col}]; 
                end
            end
            %accumulate leg tip velocity
            load('leg_tip_v.mat')
            lt_pre_all = [lt_pre_all; lt_pre];
            lt_post_all = [lt_post_all; lt_post];

            %accumulate st2sw latency data
            load(['st2sw_ltc_tw' num2str(tw(tp)) '.mat'])
            for cnd = 1:4 %condition pre/w/m/s
                for leg1 = 1:6
                    for leg2 = 1:6
                        acc_ltc_4cd_all_st2sw{cnd}{leg1,leg2} = [acc_ltc_4cd_all_st2sw{cnd}{leg1,leg2} acc_ltc_4cd{cnd}{leg1,leg2}];
                    end
                end
            end
            %accumulate sw2sw proportional latency data
            load(['sw2sw_p_ltc_tw' num2str(tw(tp)) '.mat'])
            for cnd = 1:4 %condition pre/w/m/s
                for leg1 = 1:6
                    for leg2 = 1:6
                        acc_ltc_p_4cd_all_sw2sw{cnd}{leg1,leg2} = [acc_ltc_p_4cd_all_sw2sw{cnd}{leg1,leg2} acc_ltc_p_4cd{cnd}{leg1,leg2}];
                    end
                end
            end
        end
        cd(analysis_dir)
    end
    %calculate stats
    fv_diff_all = cell(1,3);
    angv_diff_all = cell(1,3);
    for level = 1:3
        fv_diff_all{level} = fv_mean_all{1,level} - fv_mean_all{2,level};
        angv_diff_all{level} = angv_mean_all{1,level} - angv_mean_all{2,level};
    end
    acc_fv_mean = zeros(1,3); %mean pre/post
    acc_fv_std = zeros(1,3); %std pre/post
    acc_fv_median = zeros(1,3); %median pre/post
    acc_fv_iqr = zeros(1,3); %iqr pre/post
    acc_angv_mean = zeros(1,3); %mean pre/post
    acc_angv_std = zeros(1,3); %std pre/post
    acc_angv_median = zeros(1,3); %median pre/post
    acc_angv_iqr = zeros(1,3); %iqr pre/post
    for level = 1:3
        acc_fv_mean(level) = mean(fv_diff_all{level});
        acc_fv_std(level) = std(fv_diff_all{level});
        acc_fv_median(level) = median(fv_diff_all{level});
        acc_fv_iqr(level) = iqr(fv_diff_all{level});
        acc_angv_mean(level) = mean(angv_diff_all{level});
        acc_angv_std(level) = std(angv_diff_all{level});
        acc_angv_median(level) = median(angv_diff_all{level});
        acc_angv_iqr(level) = iqr(angv_diff_all{level});
    end
    %mean-pre/post, std-pre/post, median-pre/post, iqr-pre/post
    fv_stat{1} = fv_diff_all;
    fv_stat{2} = [acc_fv_mean; acc_fv_std; acc_fv_median; acc_fv_iqr];
    angv_stat{1} = angv_diff_all;
    angv_stat{2} = [acc_angv_mean; acc_angv_std; acc_angv_median; acc_angv_iqr];
    clear acc_fv_mean acc_fv_std acc_fv_median acc_fv_iqr
    clear acc_angv_mean acc_angv_std acc_angv_median acc_angv_iqr

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %save cumulative data
    cd(save_dir)
    save(['fv_stat_tw' num2str(tw(tp)) '.mat'], 'fv_stat', '-mat')
    save(['angv_stat_tw' num2str(tw(tp)) '.mat'], 'angv_stat', '-mat')
    save(['lt_all_tw' num2str(tw(tp)) '.mat'], 'lt_pre_all', 'lt_post_all', '-mat')
    save(['ltc_all_st2sw_' num2str(tw(tp)) '.mat'],'acc_ltc_4cd_all_st2sw','-mat')
    save(['ltc_all_sw2sw_' num2str(tw(tp)) '.mat'],'acc_ltc_p_4cd_all_sw2sw','-mat')

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %visualize cumulative data
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %plot 2-D density map (2-var-histogram)
    %basic settings
    ax_max = 80;
    ax_min = -5;
    bin = [3 3];
    n_method = 'count';
    %reconstruct leg tip velocity data
    acc_lt_all = cell(4,6);%4condition(p/w/m/s)x6legs
    %pre-stim
    for leg = 1:6
        for stim = 1:stn*3
            for id  = 1:length(lt_pre_all(:,1))
                if ~isempty(lt_pre_all{id,stim})
                    acc_lt_all{1,leg} = [acc_lt_all{1,leg}; lt_pre_all{id,stim}{leg}];
                end
            end
        end
    end
    %post-stim
    for leg = 1:6
        for level = 1:3
            for stim = 1+stn*(level-1):stn*level
                for id  = 1:length(lt_post_all(:,1))
                    if ~isempty(lt_post_all{id,stim})
                        acc_lt_all{level+1,leg} = [acc_lt_all{level+1,leg}; lt_post_all{id,stim}{leg}];
                    end
                end
            end
        end
    end
    %visualization (intra-leg coordination)
    h2d = figure('Position', [100, 100, 1000, 1000],'visible','off');
    for cnd = 1:4 %p/w/m/s
        for T = 1:3 %F/M/H
            subplot(3,4,cnd+4*(T-1))
            axis([ax_min ax_max ax_min ax_max])
            h2d = histogram2(acc_lt_all{cnd,T},acc_lt_all{cnd,T+3},...
                'FaceColor','flat',...
                'Edgecolor','none',...
                'Normalization',n_method,...
                'ShowEmptyBins','on',...
                'DisplayStyle','tile');
            h2d.BinWidth = bin;
            h2d.XBinLimits = [ax_min ax_max];
            h2d.YBinLimits = [ax_min ax_max];
            if cnd == 1
                caxis([0 45])
            else
                caxis([0 15])
            end
        end
    end
    %save fig
    cd(save_dir)
    intra_lt_hm = [ss_name '_intra_lt_hm_tw' num2str(tw(tp)) '.tif'];
    I = getframe(gcf);
    imwrite(I.cdata, intra_lt_hm)
    close all hidden
    %visualization (inter-leg coordination)
    h2d = figure('Position', [100, 100, 1000, 1000],'visible','off');
    ilc = [1 2; 4 5; 2 3; 5 6];
    for cnd = 1:4 %p/w/m/s
        for idx = 1:4 %inter-leg combination
            subplot(4,4,cnd+4*(idx-1))
            axis([ax_min ax_max ax_min ax_max])
            h2d = histogram2(acc_lt_all{cnd,ilc(idx,1)},acc_lt_all{cnd,ilc(idx,2)},...
                'FaceColor','flat',...
                'Edgecolor','none',...
                'Normalization',n_method,...
                'ShowEmptyBins','on',...
                'DisplayStyle','tile');
            h2d.BinWidth = bin;
            h2d.XBinLimits = [ax_min ax_max];
            h2d.YBinLimits = [ax_min ax_max];
            if cnd == 1
                caxis([0 45])
            else
                caxis([0 15])
            end
        end
    end
    %save fig
    cd(save_dir)
    inter_lt_hm = [ss_name '_inter_lt_hm_tw' num2str(tw(tp)) '.tif'];
    I = getframe(gcf);
    imwrite(I.cdata, inter_lt_hm)
    close all hidden
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %phase latency plot for accumulated data
    %plot latency distribution for 6legs(row) pre/w/m/s(column)
    hist_p1 = figure('Position', [100, 100, 1400, 700],'visible','off');
    color_idx = [239 69 40; %red
        142 181 66; %green
        250 208 26; %yellow
        63 155 177; %blue
        245 124 33; %orange
        93 63 182]; %purple
    color_idx = color_idx/255;
    for cnd = 1:4
        for leg1 = 1:6
            row_posi = [1 3 5 2 4 6];
            subplot(6,4,4*(row_posi(leg1)-1) + cnd)%subplot for 6legs
            for leg2 = 1:6
                %if leg2~=leg1
                temp_hist = histcounts(acc_ltc_4cd_all_st2sw{cnd}{leg1,leg2},(0:1:30));
                hold on
                plot(temp_hist,'Color',color_idx(leg2,:),'LineWidth',2)
                %end
            end
        end
    end
    %save fig
    cd(save_dir)
    st2sw_hist_name = [ss_name '_st2sw_hist_tw' num2str(tw(tp)) '.tif'];
    I = getframe(gcf);
    imwrite(I.cdata, st2sw_hist_name)
    close all hidden

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % phase transition characteristic
    % swing(leg1) to swing(leg2)
    %polarhistgram for neighboring pairs in 4 conditions (p/w/m/s)
    cnd_name = ['p'; 'w'; 'm'; 's'];
    for cnd = 1:4
        sw2sw_polar1 = figure('Position', [100, 100, 500, 700],'visible','off');
        %Intra pairs (F/M/H)
        for T = 1:3
            subplot(6,6,[3+12*(T-1),4+12*(T-1),9+12*(T-1),10+12*(T-1)])
            leg1 = T;
            leg2 = T+3;
            polarhistogram(acc_ltc_p_4cd_all_sw2sw{cnd}{leg1,leg2},30)
        end
        %Inter pairs
        %Left
        for T = 1:2
            subplot(6,6,[7+12*(T-1),8+12*(T-1),13+12*(T-1),14+12*(T-1)])
            leg1 = T;
            leg2 = T+1;
            polarhistogram(acc_ltc_p_4cd_all_sw2sw{cnd}{leg1,leg2},30)
        end
        %Right
        for T = 1:2
            subplot(6,6,[11+12*(T-1),12+12*(T-1),17+12*(T-1),18+12*(T-1)])
            leg1 = T+3;
            leg2 = T+4;
            polarhistogram(acc_ltc_p_4cd_all_sw2sw{cnd}{leg1,leg2},30)
        end
        %save figs
        cd(save_dir)
        sw2sw_polar_name = [ss_name '_polar_sw2sw_nb_' cnd_name(cnd) '_tw' num2str(tw(tp)) '.tif'];
        I = getframe(gcf);
        imwrite(I.cdata, sw2sw_polar_name)
    end
    %overlay post on pre
    for cnd = 2:4
        sw2sw_polar1 = figure('Position', [100, 100, 500, 700],'visible','off');
        %Intra pairs (F/M/H)
        for T = 1:3
            subplot(6,6,[3+12*(T-1),4+12*(T-1),9+12*(T-1),10+12*(T-1)])
            leg1 = T;
            leg2 = T+3;
            polarhistogram(acc_ltc_p_4cd_all_sw2sw{1}{leg1,leg2},30,'Normalization','probability')
            hold on
            polarhistogram(acc_ltc_p_4cd_all_sw2sw{cnd}{leg1,leg2},30,'Normalization','probability')
        end
        %Inter pairs
        %Left
        for T = 1:2
            subplot(6,6,[7+12*(T-1),8+12*(T-1),13+12*(T-1),14+12*(T-1)])
            leg1 = T;
            leg2 = T+1;
            polarhistogram(acc_ltc_p_4cd_all_sw2sw{1}{leg1,leg2},30,'Normalization','probability')
            hold on
            polarhistogram(acc_ltc_p_4cd_all_sw2sw{cnd}{leg1,leg2},30,'Normalization','probability')
        end
        %Right
        for T = 1:2
            subplot(6,6,[11+12*(T-1),12+12*(T-1),17+12*(T-1),18+12*(T-1)])
            leg1 = T+3;
            leg2 = T+4;
            polarhistogram(acc_ltc_p_4cd_all_sw2sw{1}{leg1,leg2},30,'Normalization','probability')
            hold on
            polarhistogram(acc_ltc_p_4cd_all_sw2sw{cnd}{leg1,leg2},30,'Normalization','probability')
        end
        %save figs
        cd(save_dir)
        sw2sw_polar_name = [ss_name '_polar_sw2sw_nb_' cnd_name(cnd) '_ovl_tw' num2str(tw(tp)) '.tif'];
        I = getframe(gcf);
        imwrite(I.cdata, sw2sw_polar_name)
    end

    %polar histogram for tripod pairs (pre)
    cnd = 1;
    sw2sw_polar2 = figure('Position', [100, 100, 500, 500],'visible','off');
    leg_pair = [1 5; 4 2; 2 6; 5 3];
    for k = 1:4
        subplot(2,2,k)
        polarhistogram(acc_ltc_p_4cd_all_sw2sw{cnd}{leg_pair(k,1),leg_pair(k,2)},30)
    end
    %save figs
    cd(save_dir)
    sw2sw_polar_name = [ss_name '_polar_sw2sw_tri_' cnd_name(cnd) '_tw' num2str(tw(tp)) '.tif'];
    I = getframe(gcf);
    imwrite(I.cdata, sw2sw_polar_name)
    close all hidden    
end








