%This Script calculates relative power for Channesl Fz, F3, F4, for six frequency bands
%it also gathers demographic information: site, pid, gender, risk, outcome
%it writes all of this information to a csv file
%NOTE: you have to have run Scotts preprocessing scripts for this to work.
%otherwise the file wont have the variables needed for this script

%use system() to establish the path to files you want to run
[~,filenames] = system('find /Volumes/seh33@uw.ed/2019_RelativePower/1BiologicalPsychiatry/Data/qcr/washington/derivatives/lossless/sub-s166/ses-m06/eeg/*_Seg*.set');
filenames = strsplit(filenames,'/Volumes');
filenames = strcat('/Volumes',filenames);
filenames = filenames(2:end);

%ROI are the channels that you want to extract power from.
%I am selecting F3, Fz, F4, which corresopnds to chan numbers 4, 5, 6
ROI = [4, 5, 6]; %4, 9, 14, 5, 10, 15, 6, 11, 16

    
for  k=1:length(filenames);
     EEG = pop_loadset(filenames{k});
     
     %the section below pulls participant information from the EEG data
     %structure, to be written to a csv file
     
     EEG.subject = EEG.filename(6:9); %loads pid
     subj = EEG.subject; %participant id
     group = EEG.group; %group is risk + outcome, high risk asd, high risk no asd, etc.
     
    Trials = num2str(EEG.trials); %converts number of trials to string so that it can be written to the csv
     
    risk = 'na';
    if strncmp(group,'HR',2);       %if EEG.group contains HR
        risk = 'HR';                %then participiant is HR: high risk
    elseif strncmp(group,'LR',2);   %if EEG.group contains LR
        risk = 'LR';                %then participant is LR: low risk
    elseif strncmp(group,'lang',2); %if EEG.group contains 'lang'
        risk = 'lang';              %then participant is lang: language risk (occurs in boston dataset only)
    else risk = 'UKN';              %if none of the above then risk is unknown
    end;

    outcome = 'na';
    if strfind(group,'UKN');    %if EEG.group contains UKN
        outcome = 'unknown';    %then outcome is unknown
    elseif strfind(group,'+');  %if EEG.group contains + then
        outcome = 'asd';        %outcome is ASD
    else;
        outcome = 'no-asd';     %otherwise outcome is no ASD
    end;
    
    if any(strcmpi(strsplit(EEG.filepath,'/'),'boston')); %if EEG.filepath contains boston
        site = 'boston'; %then site is boston
        gender = EEG.condition; %EEG.condition contains gender info
        session = num2str(EEG.session); %session is the visit age, 6m, 12m, etc
        if length(session)==1; session = ['0' session]; end;
        video = 'social'; %for boston videotype is always 'social'
    elseif any(strcmpi(strsplit(EEG.filepath,'/'),'washington')); %if EEG.filepath contains 'washington'
        site = 'washington'; %then site is washington
        gender = EEG.gender; %EEG.condition contains gender info
        session = num2str(EEG.session); %session is the visit age, 6m, 12m, etc
        if length(session)==1; session = ['0' session]; end; 
        if any(strcmpi(strsplit(EEG.filename,'_'),'socl')); %if EEG.filename contains socl
            video = 'social'; %then videotype is social
        elseif any(strcmpi(strsplit(EEG.filename,'_'),'toys')); %if EEG.filename contains 'toys'
            video = 'non-social'; %then videotype is non-social (vid of toys)
        else
            video = 'combined'; %otherwise the videotype is combined. this is if you used a seg script that combined all rest videotypes into one resting marker
        end;
    elseif any(strcmpi(strsplit(EEG.filepath,'/'),'london')); %if EEG.filepath contains 'london'
        site = 'london'; %then site is london
        if strcmp(EEG.condition,'0'); %if EEG.condition is 0 (london participant log used a 0, 1 gender coding)
            gender = 'M'; %then gender is male
        else;
            gender = 'F'; %otherwise gender is female (1)
        end;
        if any(strcmpi(strsplit(EEG.filename,'_'),'ses-m06')); %if EEG.filename contains ses-m06
            session = '06'; %visit age is 6m
            video = 'social'; %videotype is always ;social' for 6m visit
        else;
            session = '12'; %otherwise visit age is 12m
            if any(strcmpi(strsplit(EEG.filename,'_'),'rest1')); %if EEG.filename contains 'rest1'
                video = 'social'; %then videotype is social (video of person)
            elseif any(strcmpi(strsplit(EEG.filename,'_'),'rest2')); %if EEG.filename contains 'rest2'
                video = 'non-social'; %then videotype is 'non-social' (video of toy)
            elseif any(strcmpi(strsplit(EEG.filename,'_'),'rest3')); %if EEG.filename contains 'rest3'
                video = 'semi-social'; %then videotype was semi-social (hand activating a toy)
            else 
                video = 'combined'; %if none of the above, then its combined, you ran a seg script that combined 3 videotypes into a single rest marker
            end;
        end;       
    end;
     
    sRate = EEG.srate; %sampling rate
    window = sRate*4;
    NFFT = sRate*4; %making window and NFFT variables the same size
    h=hamming(size(EEG.data,2));
 
    TD_dat_power_ROI = []; %clearing this variable before each case is run 
    TD_dat_power_ROI_Rel = []; %^^
    TD_tmp_dat_power = []; %^^
    
    %TD_dat_power_ROI = zeros(size(EEG.data,2),length(ROI));
    %TD_dat_power_ROI_Rel = zeros(size(EEG.data,2),length(ROI));
    for i=1:length(ROI); %For every channel
        %TD_tmp_dat_power = zeros(size(EEG.data,2),size(EEG.data,3));
        for j=1:size(EEG.data,3); %EEG.data,3 is the number of trials. "for every trial"
            
            %pwelch call, welch method of calculating PSD
            [Pxx,F] = pwelch(EEG.data(ROI(i),:,j),window,[],NFFT,sRate); %calculating absolute power on each channel and trial. doing Fft on 256 points
            TD_tmp_dat_power(:,j) = Pxx; %frequencies rows x trials column? accumulating the spectrum for each trial
            
            %fft call, fft method of calculating PSD (for establishing replicability only, not to be written to csv)
            %dh=h.*squeeze(EEG.data(ROI(i),:,j))';
            %fq=fft(dh);
            %TD_tmp_dat_power_f(:,j)=abs(fq(1:length(TD_tmp_dat_power(:,j)))).^2;% squared (^2) to match "power" values like Pxx of pwelch.
            
        end
        
         %Get frequency bins of interest   
         if i == 1
             fndx = find(F >= 2 & F <= 50); %F is freqs    
         end
         
         TD_dat_power_ROI(:,i) = mean(TD_tmp_dat_power(:,:),2); %calculating average spectrum across trials for channel i. frequencies rows? x channels column
         
         %RelativePower= freqofinterest(:,:)./sum(freqofinterest,1) % ^^
         TD_dat_power_ROI_Rel(:,i)=TD_dat_power_ROI(:,i)./trapz(F(fndx),(TD_dat_power_ROI(fndx,i))); % F is freqs, ichan is the channel index, fndx is the vector of frequency bins
         
        %all of the rows are mean spectrum across trials for that channel
        
        %TD_dat_power_ROI_f(:,i) = mean(TD_tmp_dat_power_f(:,:),2); %doing the same but for fft (replicability only)
    end

    
    %TD_dat_power_ROI_Rel=TD_dat_power_ROI./sum(TD_dat_power_ROI(9:201,:),1); %converting to relative power. summing each channel (column).
    %freqofinterest=TD_dat_power_ROI(2:50,:); %2019/05/011 testing out code
   
    %FFT to test for replicability of pwelch
    %TD_dat_power_ROI_Rel_f(:,i)=TD_dat_power_ROI_f(:,i_)./trapz(F(fndx),(TD_dat_power_ROI(fndx,i)));
    
    %plot psd for both pwelch (blue) and fft (red)on same graph. only do this during testing, otherwise comment out 128 (figure) & 129 (hold). when plotting pwelch and fft PSDs, they should overlap almost perfectly
    %fqs=0:length(TD_dat_power_ROI_Rel)-1;
    %fqs=fqs*.25;                
    %figure;plot(fqs,TD_dat_power_ROI_Rel);
    %hold on;plot(fqs,TD_dat_power_ROI_Rel_f,'r');
    
    %plot psd for absolute power using both the welch method and fft (red)
    %figure;plot(fqs,TD_dat_power_ROI);
    %hold on;plot(fqs,TD_dat_power_ROI_f,'r');
    
    %plot LOG TRANSFORMED psd for absolute power using both the welch method and fft (red)
    %figure;plot(fqs,log10(TD_dat_power_ROI),'b');
    %hold on;plot(fqs,log10(TD_dat_power_ROI_f),'r');
    
    
    
    %log transform the relative power, and average across frontal
    %channels.It is then creating a 4th channel which is the values of the
    %log average signal
    TD_pow_Rel_avg = mean(TD_dat_power_ROI_Rel,2); % average across channels
    logTDpowRelAvg = log10(TD_pow_Rel_avg); % log10 transform 
    TD_dat_power_ROI_Rel(:,4) =  logTDpowRelAvg; %this is adding logTDpowRelAvg to the 4th column of TD_dat_power_ROI_REL so that it can be written to the csv without substantially chaning the code.
    
    %define frequency bands for ABSOLUTE power
    
    delta_ABS=(TD_dat_power_ROI(9:17,:)); %2-4hz, technically 2-3.75hz, corresponds to bins 9:17
    
    delta_ABS=sum(delta_ABS,1); %sums the values in bins above
    
    delta_ABS=num2str(sum(delta_ABS,1),'%f,'); %converts integers to string so that it can be written to the csv
    delta_ABS=delta_ABS(1:end-1); %values for chans 4,5,6 are in one cell, so this tells it how to index the one cell to get the 3 channels delta values and write it to different cells in the csv
    
    theta_ABS=(TD_dat_power_ROI(18:25,:));
    
    theta_ABS=sum(theta_ABS,1);
    
    theta_ABS=num2str(sum(theta_ABS,1),'%f,');
    theta_ABS=theta_ABS(1:end-1);
    
    alpha1_ABS=(TD_dat_power_ROI(26:37,:)); %6-9hz, technically 6.25
    
    alpha1_ABS=sum(alpha1_ABS,1);
    
    alpha1_ABS=num2str(sum(alpha1_ABS,1),'%f,');
    alpha1_ABS=alpha1_ABS(1:end-1);
    
    
    alpha2_ABS=(TD_dat_power_ROI_Rel(38:53,:)); %9-13hz, technically 9.25
    
    alpha2_ABS=sum(alpha2_ABS,1);
    
    alpha2_ABS=num2str(sum(alpha2_ABS,1),'%f,');
    alpha2_ABS=alpha2_ABS(1:end-1);
    
    beta_ABS=(TD_dat_power_ROI(54:121,:)); %13-30hz, technically 30.25
    
    beta_ABS=sum(beta_ABS,1);
    
    beta_ABS=num2str(sum(beta_ABS,1),'%f,');
    beta_ABS=beta_ABS(1:end-1);
    
    gamma_ABS=(TD_dat_power_ROI(122:201,:)); %30hz - 50hz, technically 30.25
    
    gamma_ABS=sum(gamma_ABS,1);
    
    gamma_ABS=num2str(sum(gamma_ABS,1),'%f,');
    gamma_ABS=gamma_ABS(1:end-1);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %define frequency bins for RELATIVE power
    delta=(TD_dat_power_ROI_Rel(9:17,:)); %2-4hz, technically 2-3.75hz, corresponds to bins 9:17
    
    delta=sum(delta,1); %sums the values in bins above
    
    delta=num2str(sum(delta,1),'%f,'); %converts integers to string so that it can be written to the csv
    delta=delta(1:end-1); %values for chans 4,5,6 are in one cell, so this tells it how to index the one cell to get the 3 channels delta values and write it to different cells in the csv
    
    theta=(TD_dat_power_ROI_Rel(18:25,:));
    
    theta=sum(theta,1);
    
    theta=num2str(sum(theta,1),'%f,');
    theta=theta(1:end-1);
    
    alpha1=(TD_dat_power_ROI_Rel(26:37,:)); %6-9hz, technically 6.25
    
    alpha1=sum(alpha1,1);
    
    alpha1=num2str(sum(alpha1,1),'%f,');
    alpha1=alpha1(1:end-1);
    
    
    alpha2=(TD_dat_power_ROI_Rel(38:53,:)); %9-13hz, technically 9.25
    
    alpha2=sum(alpha2,1);
    
    alpha2=num2str(sum(alpha2,1),'%f,');
    alpha2=alpha2(1:end-1);
    
    beta=(TD_dat_power_ROI_Rel(54:121,:)); %13-30hz, technically 30.25
    
    beta=sum(beta,1);
    
    beta=num2str(sum(beta,1),'%f,');
    beta=beta(1:end-1);
    
    gamma=(TD_dat_power_ROI_Rel(122:201,:)); %30hz - 50hz, technically 30.25
    
    gamma=sum(gamma,1);
    
    gamma=num2str(sum(gamma,1),'%f,');
    gamma=gamma(1:end-1);
                    
    % RELATIVE POWER CSV FILE: create headers for the csv file if the csv file doesn't exist yet
    if ~exist('PowerChecks166.csv');    
         fid = fopen('Powerchecks166.csv','w');

        header = {'site','pid','age','gender','risk', ...
                  'outcome','group','video','nTrials', ...
                  'Delta_LF','Delta_MF','Delta_RF','Delta_LA', ...
                  'Theta_LF','Theta_MF','Theta_RF','Theta_LA', ...
                  'Alpha1_LF','Alpha1_MF','Alpha1_RF','Alpha1_LA', ...
                  'Alpha2_LF','Alpha2_MF','Alpha2_RF','Alpha2_LA', ...
                  'Beta_LF','Beta_MF','Beta_RF','Beta_LA', ...
                  'Gamma_LF','Gamma_MF','Gamma_RF','Gamma_LA'};
        fprintf(fid,'%s,',header{1,1:end-1});
        fprintf(fid,'%s\n',header{1,end});
    end;
    
    
    
    %this section assigns the variables to be written to the corresponding
    %headers. make sure order of headers and order of data is the same
    data={site,subj,session,gender,risk,outcome,group,video,Trials,delta,theta,alpha1,alpha2,beta,gamma};
    fprintf(fid,'%s,',data{1,1:end-1});
    fprintf(fid,'%s\n',data{1,end});
    
      % ABSOLUTE POWER CSV FILE: create headers for the csv file if the csv file doesn't exist yet
    %if ~exist('absolutePower.csv');    
     %    fid = fopen('absolutepower.csv','w');

      %  header = {'site','pid','age','gender','risk', ...
      %            'outcome','group','video','nTrials' ...
      %            'Delta_LF','Delta_MF','Delta_RF', ...
      %            'Theta_LF','Theta_MF','Theta_RF', ...
      %            'Alpha1_LF','Alpha1_MF','Alpha1_RF', ...
      %            'Alpha2_LF','Alpha2_MF','Alpha2_RF', ...
      %            'Beta_LF','Beta_MF','Beta_RF', ...
      %            'Gamma_LF','Gamma_MF','Gamma_RF'};
      %  fprintf(fid,'%s,',header{1,1:end-1});
      %  fprintf(fid,'%s\n',header{1,end});
    %end;
    
    
    
    %this section assigns the variables to be written to the corresponding
    %headers. make sure order of headers and order of data is the same
   % data={site,subj,session,gender,risk,outcome,group,video,Trials,delta_ABS,theta_ABS,alpha1_ABS,alpha2_ABS,beta_ABS,gamma_ABS};
   % fprintf(fid,'%s,',data{1,1:end-1});
   % fprintf(fid,'%s\n',data{1,end});
    
    % TD_dat_power_ROI_subs(:,:,k) = TD_dat_power_ROI;
end
    fclose(fid);
