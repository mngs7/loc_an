clear;
warning off;

%time window for analysis
tw = [10; 30; 50; 150; 300];%different time windows
tp = 2;
%set analysis_dir
analysis_dir = uigetdir();
cd(analysis_dir)

SS_lines = struct2cell(dir('JRC_*'))';
posi_ctrl = struct2cell(dir('*VGLUT*'))';
nega_ctrl = struct2cell(dir('*YNA*'))';
lines = [SS_lines; posi_ctrl; nega_ctrl];
clear SS_lines posi_ctrl nega_ctrl
lines(:, 2:end) = [];
stat = lines;
for k = 1:length(lines(:,1))
    cd(analysis_dir)
    cd(lines{k,1})
    load(['fv_stat_tw' num2str(tw(tp)) '.mat'])
    stat{k,2} = fv_stat;
    load(['angv_stat_tw' num2str(tw(tp)) '.mat'])
    stat{k,3} = angv_stat;
end

%visualization
cd(analysis_dir)
level = {'w','m','s'};
for n = 1:3 %w/m/s
    %fv
    temp_fv = lines;
    for k = 1:length(lines(:,1))
        temp_fv{k,2} = stat{k,2}{1}{n};
        temp_fv{k,3} = stat{k,2}{2}(3,n);
        temp_fv{k,4} = stat{k,2}{2}(4,n);
    end    
    temp_fv = sortrows(temp_fv,3);
    x = [];
    g = [];
    for k = 1:length(temp_fv(:,1))
        x = [x; temp_fv{k,2}'];
        g = [g; repmat(temp_fv{k,1},length(temp_fv{k,2}),1)];
    end
    figure
    boxplot(x,g)
    fv_bp = ['fv_bp_tw' num2str(tw(tp)) '_' level{n} '.tif'];
    I = getframe(gcf);
    imwrite(I.cdata, fv_bp)

    %angv
    temp_angv = lines;
    for k = 1:length(lines(:,1))
        temp_angv{k,2} = stat{k,3}{1}{n};
        temp_angv{k,3} = stat{k,3}{2}(3,n);
        temp_angv{k,4} = stat{k,3}{2}(4,n);
    end
    temp_angv = sortrows(temp_angv,3);
    x = [];
    g = [];
    for k = 1:length(temp_angv(:,1))
        x = [x; temp_angv{k,2}'];
        g = [g; repmat(temp_angv{k,1},length(temp_angv{k,2}),1)];
    end
    figure
    boxplot(x,g)
    angv_bp = ['angv_bp_tw' num2str(tw(tp)) '_' level{n} '.tif'];
    I = getframe(gcf);
    imwrite(I.cdata, angv_bp)
end








