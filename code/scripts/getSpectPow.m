%This Script calculates relative power for Channesl Fz, F3, F4, for six frequency bands
%it also gathers demographic information: site, pid, gender, risk, outcome
%it writes all of this information to a csv file
%NOTE: you have to have run Scotts preprocessing scripts for this to work.
%otherwise the file wont have the variables needed for this script

%This function loads EEGLAB *.set files and extracts features that are saved 
%to CSV file. The fnamefile input is a text file containing the filenames
%with paths to each of the *.set files to be loaded (one file per line).
%From a bash terminal these file name text files can be created as follows:
%for example...
%find . -type f -name "*-SEGrest*.set" > derivatives/seg_rest_spect/code/misc/fnames.txt
%The outfname input is the name of the file containing the output features.

function getSpectPow(fnamefile,outfname)


fidIn=fopen(fnamefile,'r');
fnd=fread(fidIn);
fclose(fidIn);
cr_ind=find(fnd==10);
for i=1:length(cr_ind);
   if i==1;
      c_fname{i}=deblank(char(fnd(3:cr_ind(i))'));
   else
      c_fname{i}=deblank(char(fnd(cr_ind(i-1)+3:cr_ind(i))'));
   end
end

%use system() to establish the path to files you want to run
%[~,filenames] = system('find /Volumes/seh33@uw.ed/2019_RelativePower/1BiologicalPsychiatry/Data/qcr/boston/derivatives/lossless/sub-s3069/ses-m06/eeg/*_Seg*.set');
%filenames = strsplit(filenames,'/Volumes');
%filenames = strcat('/Volumes',filenames);
%filenames = filenames(2:end);

%ROI are the channels that you want to extract power from.
%I am selecting F3, Fz, F4, which corresopnds to chan numbers 4, 5, 6
ROI = [4, 5, 6]; %4, 9, 14, 5, 10, 15, 6, 11, 16



fidOut = fopen(outfname,'w');       
header = ['site,pid,age,gender,risk,', ...
          'outcome,group,', ...
          'video,nTrials,sRate,',...
          'absDeltaLog,absThetaLog,absAlpha1Log,absAlpha2Log,absBetaLog,absGammaLog,',...
          'relDeltaLog,relThetaLog,relAlpha1Log,relAlpha2Log,relBetaLog,relGammaLog,',...
          'trapzAbsDelta,trapzAbsTheta,trapzAbsAlpha1,trapzAbsAlpha2,trapzAbsBeta,trapzAbsGamma,',...
          'trapzRelDelta,trapzRelTheta,trapzRelAlpha1,trapzRelAlpha2,trapzRelBeta,trapzRelGamma,',...
          'F3AbsDelta,FzAbsDelta,F4AbsDelta,',...
          'F3AbsTheta,FzAbsTheta,F4AbsTheta,',...
          'F3AbsAlpha1,FzAbsAlpha1,F4AbsAlpha1,',...
          'F3AbsAlpha2,FzAbsAlpha2,F4AbsAlpha2,',...
          'F3AbsBeta,FzAbsBeta,F4AbsBeta,',...
          'F3AbsGamma,FzAbsGamma,F4AbsGamma,',...
          'F3RelDelta,FzRelDelta,F4RelDelta,',...
          'F3RelTheta,FzRelTheta,F4RelTheta,',...
          'F3RelAlpha1,FzRelAlpha1,F4RelAlpha1,',...
          'F3RelAlpha2,FzRelAlpha2,F4RelAlpha2,',...
          'F3RelBeta,FzRelBeta,F4RelBeta,',...
          'F3RelGamma,FzRelGamma,F4RelGamma,'];     
fprintf(fidOut,'%s\n',header);

    
for  k=1:length(c_fname);
     EEG = pop_loadset(c_fname{k});
     
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
    
    p=what(EEG.filepath);
    inAbsPath=p.path;
    if strfind(inAbsPath,'boston'); %if EEG.filepath contains boston
        site = 'boston'; %then site is boston
        gender = EEG.condition; %EEG.condition contains gender info
        session = num2str(EEG.session); %session is the visit age, 6m, 12m, etc
        if length(session)==1; session = ['0' session]; end;
        video = 'social'; %for boston videotype is always 'social'
    elseif strfind(EEG.filepath,'washington'); %if EEG.filepath contains 'washington'
        site = 'washingtinAbsPathon'; %then site is washington
        gender = EEG.gender; %EEG.condition contains gender info
        session = num2str(EEG.session); %session is the visit age, 6m, 12m, etc
        if length(session)==1; session = ['0' session]; end; 
        if strfind(EEG.filename,'socl'); %if EEG.filename contains socl
            video = 'social'; %then videotype is social
        elseif strfind(EEG.filename,'toys'); %if EEG.filename contains 'toys'
            video = 'non-social'; %then videotype is non-social (vid of toys)
        else
            video = 'combined'; %otherwise the videotype is combined. this is if you used a seg script that combined all rest videotypes into one resting marker
        end;
    elseif strfind(inAbsPath,'london'); %if EEG.filepath contains 'london'
        site = 'london'; %then site is london
        if strcmp(EEG.condition,'0'); %if EEG.condition is 0 (london participant log used a 0, 1 gender coding)
            gender = 'M'; %then gender is male
        else;
            gender = 'F'; %otherwise gender is female (1)
        end;
        if strfind(EEG.filename,'ses-m06'); %if EEG.filename contains ses-m06
            session = '06'; %visit age is 6m
            video = 'social'; %videotype is always ;social' for 6m visit
        else;
            session = '12'; %otherwise visit age is 12m
            if strfind(EEG.filename,'rest1'); %if EEG.filename contains 'rest1'
                video = 'social'; %then videotype is social (video of person)
            elseif strfind(EEG.filename,'rest2'); %if EEG.filename contains 'rest2'
                video = 'non-social'; %then videotype is 'non-social' (video of toy)
            elseif strfind(EEG.filename,'rest3'); %if EEG.filename contains 'rest3'
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
    
    %%CLEAR VARIABLES BETWEEN PARTICIPANTS
    TD_dat_power_ROI = []; %clearing this variable before each case is run 
    TD_dat_power_ROI_Rel = []; %^^
    TD_tmp_dat_power = []; %^^

    TD_dat_power_ROI_f = []; %clearing this variable before each case is run 
    TD_dat_power_ROI_Rel_f = []; %^^
    TD_tmp_dat_power_f = []; %^^

    relDenom = [];
    
    delta = [];
    deltaRel = [];
    avgDelta = [];
    avgDeltaRel = [];
    deltaLog10 = [];
    deltaRelLog10 = [];

    theta = [];
    thetaRel = [];
    avgTheta = [];
    avgThetaRel = [];
    thetaLog10 = [];
    thetaRelLog10 = [];
    
    alpha1 = [];
    alpha1Rel = [];
    avgAlpha1 = [];
    avgAlpha1Rel = [];
    alpha1Log10 = [];
    alpha1RelLog10 = [];
  
    alpha2 = [];
    alpha2Rel = [];
    avgAlpha2 = [];
    avgAlpha2Rel = [];
    alpha2Log10 = [];
    alpha2RelLog10 = [];

    beta = [];
    betaRel = [];
    avgBeta =  [];
    avgBetaRel = [];
    betaLog10 = [];
    betaRelLog10 = [];

    gamma = [];
    gammaRel = [];
    avgGamma = [];
    avgGammaRel = [];
    gammaLog10 = [];
    gammaRelLog10 = [];

    absLog2Str = [];
    relLog2Str = [];
    avgAbsTrapz2Str = [];
    avgRelTrapz2Str = [];
    absPowerAllChans2Str = [];
    relPowerAllChans2Str = [];

    
    for i=1:length(ROI); %For every channel
        
        for j=1:size(EEG.data,3); %EEG.data,3 is the number of trials. "for every trial"
            
            %pwelch call, welch method of calculating PSD
            [Pxx,F] = pwelch(EEG.data(ROI(i),:,j),window,[],NFFT,sRate); %calculating absolute power on each channel and trial. doing Fft on 256 points
            TD_tmp_dat_power(:,j) = Pxx; %frequencies rows x trials column? accumulating the spectrum for each trial
            
            %fft call, fft method of calculating PSD (for establishing replicability only, not to be written to csv)
            dh=h.*squeeze(EEG.data(ROI(i),:,j))';
            fq=fft(dh);
            TD_tmp_dat_power_f(:,j)=abs(fq(1:length(TD_tmp_dat_power(:,j)))).^2;% squared (^2) to match "power" values like Pxx of pwelch.
            
        end
        
         %Get frequency bins of interest   
         if i == 1
             fndx = find(F >= 2 & F <= 50); %F is freqs    
         end
         
         TD_dat_power_ROI(:,i) = mean(TD_tmp_dat_power(:,:),2); %calculating average spectrum across trials for channel i. frequencies rows? x channels column
         TD_dat_power_ROI_f(:,i) = mean(TD_tmp_dat_power_f(:,:),2); %doing the same but for fft (replicability only)
         
         %RelativePower that will be used for plotting figures... = freqofinterest(:,:)./sum(freqofinterest,1) % ^^
         TD_dat_power_ROI_Rel(:,i)=TD_dat_power_ROI(:,i)./trapz(F(fndx),(TD_dat_power_ROI(fndx,i))); % F is freqs, ichan is the channel index, fndx is the vector of frequency bins
         TD_dat_power_ROI_Rel_f(:,i)= TD_dat_power_ROI_f(:,i)./trapz(F(fndx),(TD_dat_power_ROI_f(fndx,i))); %^^ for FFT replicability test
         
        %all of the rows are mean spectrum across trials for that channel
        
        
    end

    %plot psd for both pwelch (blue) and fft (red)on same graph. 
    %only do this during testing, otherwise COMMENT OUT figure & hold so it doesnt generate a plot for every single file
    %when plotting pwelch and fft PSDs, they should overlap almost perfectly
    fqs=0:length(TD_dat_power_ROI_Rel)-1;
    fqs=fqs*.25;                
    
    %Relative Power using both the welch method (blue) and fft (red)
    figure;plot(fqs,mean(TD_dat_power_ROI_Rel,2));
    hold on;plot(fqs,mean(TD_dat_power_ROI_Rel_f,2),'r');
    title(['Relative Power ', EEG.filename]);
    
    
    %absolute power using both the welch method (blue) and fft (red)
    %James can you help here? I can't figure out how to plot the
    %trapezoidal absolute power!
       
    figure;plot(fqs,mean(TD_dat_power_ROI,2));
    hold on;plot(fqs,mean(TD_dat_power_ROI,2),'r');
    title(['absolute power ', EEG.filename]); 
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %define frequency bins for power
    %delta=(TD_dat_power_ROI_Rel(9:17,:)); %2-4hz, technically 2-3.75hz, corresponds to bins 9:17
    %delta=sum(delta,1); %sums the values in bins above
    
    relDenom = trapz(F(fndx),TD_dat_power_ROI(fndx,:));  %the power of 2-50hz (our freq range of interest) using the trapezoidal method. It will be the denominator of every bands relative power calculation    
  
    %DELTA
    delta = trapz(F(9:17),TD_dat_power_ROI(9:17,:)); %the absolute power of delta, using the trapezoidal method. %2-4hz corresponds to bins 9:17
    deltaRel = delta./relDenom; %the relative power of delta. dividing delta power by the power of our freq range of interest
    avgDelta = mean(delta,2); %average across channels for absolute delta
    avgDeltaRel = mean(deltaRel,2); %average across channels for relative delta
    deltaLog10 = log10(avgDelta); %log 10 transformation of averaged absolute power of delta
    deltaRelLog10 = log10(avgDeltaRel); %log 10 transformation of averaged relative power of delta
          
    %THETA
    theta = trapz(F(18:25),TD_dat_power_ROI(18:25,:)); %the absolute power of theta, using the trapezoidal method. %4-6hz corresponds to bins 18:25
    thetaRel = theta./relDenom;  %the relative power of theta. dividing theta power by the power of our freq range of interest
    avgTheta = mean(theta,2); %average across channels for absolute theta
    avgThetaRel = mean(thetaRel,2); %average across channels for relative theta
    thetaLog10 = log10(avgTheta); %log 10 transformation of averaged absolute power of theta
    thetaRelLog10 = log10(avgThetaRel); %log 10 transformation of averaged relative power of theta
    
    %LOW ALPHA
    alpha1 = trapz(F(26:37),TD_dat_power_ROI(26:37,:)); %the absolute power of Low Alpha, using the trapezoidal method. %6-8hz corresponds to bins 26:37
    alpha1Rel = alpha1./relDenom; %the relative power of Low Alpha, using the trapezoidal method. dividing alpha1 power by the power of our freq range of interest
    avgAlpha1 = mean(alpha1,2); %average across channels for absolute Low Alpha
    avgAlpha1Rel = mean(alpha1Rel,2); %average across channels for relative theta
    alpha1Log10 = log10(avgAlpha1); %log 10 transformation of averaged absolute power of Low Alpha
    alpha1RelLog10 = log10(avgAlpha1Rel); %log 10 transformation of averaged relative power of Low Alpha
    
    %HIGH AlPHA
    alpha2 = trapz(F(38:53),TD_dat_power_ROI(38:53,:)); %the absolute power of High Alpha, using the trapezoidal method. %9-13hz corresponds to bins 38:53
    alpha2Rel = alpha2./relDenom; %the relative power of High Alpha, using the trapezoidal method. dividing high alpha power by the power of our freq range of interest
    avgAlpha2 = mean(alpha2,2); %average across channels for absolute High Alpha
    avgAlpha2Rel = mean(alpha2Rel,2); %average across channels for relative High Alpha
    alpha2Log10 = log10(avgAlpha2); %log 10 transformation of averaged absolute power of High Alpha
    alpha2RelLog10 = log10(avgAlpha2Rel); %log 10 transformation of averaged relative power of High Alpha
    
    %BETA
    beta = trapz(F(54:121),TD_dat_power_ROI(54:121,:));%the absolute power of beta, using the trapezoidal method. 13-30hz corresponds to bins 54:121
    betaRel = beta./relDenom; %the relative power of beta, using the trapezoidal method. dividing beta power by the power of our freq range of interest
    avgBeta =  mean(beta,2); %average across channels for absolute beta
    avgBetaRel = mean(betaRel,2); %average across channels for relative Beta
    betaLog10 = log10(avgBeta); %log 10 transformation of averaged absolute power of Beta
    betaRelLog10 = log10(avgBetaRel);
    
    %GAMMA
    gamma = trapz(F(122:201),TD_dat_power_ROI(122:201,:));%the absolute power of gamma, using the trapezoidal method. 30-50hz corresponds to bins 54:121
    gammaRel = gamma./relDenom; %the relative power of gamma, using the trapezoidal method. dividing beta power by the power of our freq range of interest
    avgGamma = mean(gamma,2); %average across channels for absolute gamma
    avgGammaRel = mean(gammaRel,2); %average across channels for relative Gamma
    gammaLog10 = log10(avgGamma); %log 10 transformation of averaged absolute power of Gamma
    gammaRelLog10 = log10(avgGammaRel);
         
    
    %Convert Values into string so that they can be written to the CSVfile. THis also organizes the variables in a intuitive way for when
    %they are written to the CSV
  
    absLog2Str = num2str([deltaLog10,thetaLog10,alpha1Log10,alpha2Log10,betaLog10,gammaLog10], '%f,');
    relLog2Str = num2str([deltaRelLog10,thetaRelLog10,alpha1RelLog10,alpha2RelLog10,betaRelLog10,gammaRelLog10], '%f,');
    avgAbsTrapz2Str = num2str([avgDelta,avgTheta,avgAlpha1,avgAlpha2,avgBeta,avgGamma], '%f,');
    avgRelTrapz2Str = num2str([avgDeltaRel,avgThetaRel,avgAlpha1Rel,avgAlpha2Rel,avgBetaRel,avgGammaRel], '%f,');
    absPowerAllChans2Str = num2str([delta,theta,alpha1,alpha2,beta,gamma], '%f,');
    relPowerAllChans2Str = num2str([deltaRel,thetaRel,alpha1Rel,alpha2Rel,betaRel,gammaRel], '%f,');
    sRate = num2str(EEG.srate); %convert to string so that it can be written to csv
    
    %POWER CSV FILE: create headers for the csv file if the csv file doesn't exist yet
    
    %Assigns the matlab variables that will be written to the csv file
    fprintf(fidOut,['%s,%s,%s,%s,%s,%s,%s,%s,%s,%s,',...
                '%s%s%s%s%s%s\n'],...
                site,subj,session,gender,risk,outcome,group,video,Trials,sRate,...
                absLog2Str,...
                relLog2Str,...
                avgAbsTrapz2Str,...
                avgRelTrapz2Str,...
                absPowerAllChans2Str,...
                relPowerAllChans2Str);
    
end
fclose(fidOut);
