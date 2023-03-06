close all
clear variables
% clc

% Allison Cutri & Lauren Eckermann
% Creation date: January 28, 2022
% Analyze potential modulation in bacteria and determine if fluorescence
% signals are statistically significant

format long
font = 12;
linewidth = 1.5;

tic % Start timer for code runtime

% Load text file of compiled processed time traces

addpath('E:\Wide Field\2023\February\02132023\100pc succinate 24hr\MPS6');

filename1 = '100pc succinate laser 20mW MPS6.tiff_data.txt';

% Load the time trace structure into a numerical matrix

time_trace = uiimport(filename1); % Open the import wizard
time_trace = cell2mat(struct2cell(time_trace)); % Transform structure array to numerical matrix

time = time_trace(:,1);
time_trace (:,1) = []; % Delete the first column because it's only the time data
  
[m,n] = size(time_trace); % number of rows and columns in the time trace matrix

%%
% Extract each time trace into a histogram

for i = 1:n
    g = figure (); % opens a blank figure
    
    h = histogram(time_trace(:,i),'FaceColor','k'); % plots a histogram in g of the time trace data
    bins(:,i)=h.NumBins;
    xlabel('Intensity (a.u.)')
    ylabel('Occurances')
    
    saveas(g,sprintf('histogram pore-%d.png',i)); % saves plots to current folder
    close(g);
    
end

%%

% Determine mean value, separate data based on where it falls compared to
% the mean

mean_vals=mean(time_trace);

for column=1:n
    time_trace_on_i=1;
    time_trace_off_i=1;
    for row=1:m
        if time_trace(row,column) > mean_vals(:,column) % test if number is greater than the mean value and "on"
            time_trace_on(time_trace_on_i,column)=time_trace(row,column); % store values to new matrix
            time_trace_on_i=time_trace_on_i+1; % add one to avoid indexing errors
            
        else
            time_trace_off(time_trace_off_i,column)=time_trace(row,column); % all other data is considered to be less than the mean and "off"
            time_trace_off_i=time_trace_off_i+1; % add one to avoid indexing errors
        end
    end
end

time_trace_off(time_trace_off==0)=NaN;
time_trace_on(time_trace_on==0)=NaN;

%%
% Create histogram plot overlaying data sets

int_quart_range=iqr(time_trace); % calculate the interquartile range
bin_width=(2*int_quart_range)/(600^(1/3)); % calculate the estimated bin width based on Freedman-Diaconis rule

for i=1:n
    
    g=figure();   
    
    h1=histogram(time_trace_off(:,i),'BinWidth',bin_width(:,i),'FaceColor','k','EdgeColor','k');
    hold on
    h2=histogram(time_trace_on(:,i),'BinWidth',bin_width(:,i),'FaceColor','r','EdgeColor','r');
    xlabel('Intensity (a.u.)')
    ylabel('Occurances')
    
    h1.Normalization='probability';
    h2.Normalization='probability';
    
    saveas(g,sprintf('separated histogram pore-%d.png',i)); % saves plots to current folder
    close(g);
    
end


%%
% 
% for i=1:n
%    
%     pd(:,i)=fitdist(time_trace(:,i),'kernel');
%     [h(:,i),p(:,i),st(:,i)]=chi2gof(time_trace(:,i),'CDF',pd(:,i));
%     chi_results(:,i)=st(:,i).chi2stat;
%     
% end

%%

% Calculate the moving average based on the datasets and plot over the
% original timetrace

mu_mov=movmean(time_trace,5);
stdev_mov=movstd(time_trace,5);

up_bound=mu_mov+stdev_mov;
lo_bound=mu_mov-stdev_mov;

for i=1:n
    
    g = figure (); % opens a blank figure
  
    plot(time, time_trace(:,i), 'color','k','linewidth',1.5)
    hold on
    plot(time, mu_mov(:,i), 'color','r','linewidth',1.5)
    hold off
    
    xlabel('Time (s)')
    ylabel('Intensity (a.u.)')
    legend('Intensity timetrace','Moving average','location','best')
    
    saveas(g,sprintf('moving average timetrace pore-%d.png',i)); % saves plots to current folder
    close(g)
    
end

%%

time_trace_3_sig=time_trace;

logical_array=time_trace_3_sig >= up_bound | time_trace_3_sig <= lo_bound;
time_trace_3_sig(~logical_array)=NaN;
    

for i=1:n
    
    g = figure ();
    
    plot(time, mu_mov(:,i), 'color','r','linewidth',1.5)
    hold on
    
    plot(time, time_trace_3_sig(:,i),'o','color','b')
    hold off
    
    xlabel('Time (s)')
    ylabel('Intensity (a.u.)')
    legend('Moving average','Outliers at 3\sigma','location','best')
    
    saveas(g,sprintf('only moving average with 3sigma points pore-%d.png',i)); % saves plots to current folder
    close(g)
    
end

%%
% Calculate means at potential extremes, store as one var

stdev_mov_3=3*stdev_mov;
sterr_mov=stdev_mov_3/sqrt(m);

[h,p,ci,stats]=ttest(mu_mov);

M10=mean(mu_mov(97:103,:));
M20=mean(mu_mov(197:203,:));
M30=mean(mu_mov(297:303,:));
M40=mean(mu_mov(397:403,:));
M50=mean(mu_mov(497:503,:));
M60=mean(mu_mov(597:600,:));

M = [M10; M20; M30; M40; M50; M60];

%%
% Compare means at potential extremes to confidence intervals calculated earlier

Y = M < ci(1,:) | M > ci(2,:);

% for column=1:n
%     outside_ci_i=1;
%     for row=1:m
%         if time_trace(row,column) > ci(2,column) | time_trace(row,column) < ci(1,column) % test if number is greater than the mean value and "on"
%             outside_ci(outside_ci_i,column)=time_trace(row,column); % store values to new matrix
%             outside_ci_i=outside_ci_i+1; % add one to avoid indexing errors
%             
%         end
%     end
%     
%     outside-ci_time=time;
% end


% creating logical matrix with above conditions

%%

% Perform a baseline correction using 

mean_vals=mean(time_trace);
stdev_vals=std(time_trace);
stdev_3=3*stdev_vals;

for i=1:n

    [peak,loc]=findpeaks(time_trace(:,i),'MinPeakHeight',stdev_3(:,i));
    peaks(1:size(peak),i)=peak;    
    time_peaks(1:size(loc),i)=time(loc);
    
end

peaks(peaks==0)=NaN;
time_peaks(time_peaks==0)=NaN;

%%

for i=1:n
    
    g = figure();
    
    plot(time,time_trace(:,i),'k','linewidth',1.5)
    hold on
    plot(time_peaks(:,i),peaks(:,i),'o','MarkerFaceColor','r','MarkerEdgeColor','r')
    
    xlabel('Time (s)')
    ylabel('Intensity (a.u.)')
    legend('Intensity timetrace','Transient peaks','location','best')
    
    saveas(g,sprintf('peak analysis pore-%d.png',i)); % saves plots to current folder
    close(g)

end


%%

toc