%This function loads EEGLAB *.set files and extracts features that are saved 
%to CSV file. The fnamefile input is a text file containing the filenames
%with paths to each of the *.set files to be loaded (one file per line).
%From a bash terminal these file name text files can be created as follows:
%for example...
%find . -type f -name "*-SEG*.set" > derivatives/seg_rest_spect/code/misc/fnames.txt
%The outfname input is the name of the file containing the output features.

function cont=getSpectPVAF(fnamefile,outfname)


fid=fopen(fnamefile,'r');
fnd=fread(fid);
fclose(fid);
cr_ind=find(fnd==10);
for i=1:length(cr_ind);
   if i==1;
      c_fname{i}=deblank(char(fnd(3:cr_ind(i))'));
      %EEG=eeg_checkset(EEG);
      %ALLEEG(i)=EEG;
   else
      c_fname{i}=deblank(char(fnd(cr_ind(i-1)+3:cr_ind(i))'));
      %EEG=eeg_checkset(EEG);
      %ALLEEG(i)=EEG;
   end
end

cont=zeros(length(c_fname),100);
for fi=1:length(c_fname);
    %fprintf(fid,'\n%s,',c_fname{fi});
    EEG=pop_loadset('filename',c_fname{fi});
    EEG=eeg_checkset(EEG);
    figure;
    [spect,freq,spectc,cont(fi,[1:min(size(EEG.icaweights))]),std]=pop_spectopo(EEG,0,[0 3996],'EEG', ...
                'freq',[4], ...
                'plotchan',0, ...
                'icacomps',[1:min(size(EEG.icaweights))], ...
                'nicamaps',4, ...
                'freqrange',[1 30], ...
                'electrodes','off');

end

figure;plot(sort(cont'));

return


fid=fopen(fnamefile,'r');
fnd=fread(fid);
fclose(fid);
cr_ind=find(fnd==10);
for i=1:length(cr_ind);
   if i==1;
      c_fname{i}=deblank(char(fnd(3:cr_ind(i))'));
      %EEG=eeg_checkset(EEG);
      %ALLEEG(i)=EEG;
   else
      c_fname{i}=deblank(char(fnd(cr_ind(i-1)+3:cr_ind(i))'));
      %EEG=eeg_checkset(EEG);
      %ALLEEG(i)=EEG;
   end
end


frontal_labs={'AF8','Fpz','AF7','F3','Fz','F4'};
central_labs={'C3','Cz','C4'};
posterior_labs={'P3','Pz','P4','PO7','PO8','Oz'};
temporal_labs={'FT8','TP8','FT7','TP7'};

fid=fopen(outfname,'w');
fprintf(fid,'%s',['fname,',...
                'frontal_m_area,central_m_area,posterior_m_area,temporal_m_area,', ...
                'frontal_rms_area,central_rms_area,posterior_rms_area,temporal_rms_area,', ...
                'frontal_t_area,central_t_area,posterior_t_area,temporal_t_area,', ...
                'frontal_peak_p,central_peak_p,posterior_peak_p,temporal_peak_p,', ...
                'frontal_peak_p_lat,central_peak_p_lat,posterior_peak_p_lat,temporal_peak_p_lat,', ...
                'frontal_peak_n,central_peak_n,posterior_peak_n,temporal_peak_n,', ...
                'frontal_peak_n_lat,central_peak_n_lat,posterior_peak_n_lat,temporal_peak_n_lat,',...
                'frontal_t1_peak_p,central_t1_peak_p,posteriort_t1_peak_p,temporal_t1_peak_p,', ...
                'frontal_t1_peak_p_lat,central_t1_peak_p_lat,posterior_t1_peak_p_lat,temporal_t1_peak_p_lat,', ...
                'frontal_t2_peak_p,central_t2_peak_p,posteriort_t2_peak_p,temporal_t2_peak_p,', ...
                'frontal_t2_peak_p_lat,central_t2_peak_p_lat,posterior_t2_peak_p_lat,temporal_t2_peak_p_lat,', ...
                'frontal_t3_peak_p,central_t3_peak_p,posteriort_t3_peak_p,temporal_t3_peak_p,', ...
                'frontal_t3_peak_p_lat,central_t3_peak_p_lat,posterior_t3_peak_p_lat,temporal_t3_peak_p_lat,', ...
                'frontal_t4_peak_p,central_t4_peak_p,posteriort_t4_peak_p,temporal_t4_peak_p,', ...
                'frontal_t4_peak_p_lat,central_t4_peak_p_lat,posterior_t4_peak_p_lat,temporal_t4_peak_p_lat']);

for fi=1:length(c_fname);
    fprintf(fid,'\n%s,',c_fname{fi});
    EEG=pop_loadset('filename',c_fname{fi});
    EEG=eeg_checkset(EEG);
    
    frontal_inds=[];
    central_inds=[];
    posterior_inds=[];
    temporal_inds=[];
    
    for li=1:length(frontal_labs);
        frontal_inds(li)=find(strcmp({EEG.chanlocs.labels},frontal_labs{li}));
    end
    for li=1:length(central_labs);
        central_inds(li)=find(strcmp({EEG.chanlocs.labels},central_labs{li}));
    end
    for li=1:length(posterior_labs);
        posterior_inds(li)=find(strcmp({EEG.chanlocs.labels},posterior_labs{li}));
    end
    for li=1:length(temporal_labs);
        temporal_inds(li)=find(strcmp({EEG.chanlocs.labels},temporal_labs{li}));
    end
    
    ERP=mean(EEG.data,3);
    ERP_r=[];
    ERP_r(1,:)=mean(ERP(frontal_inds,:),1);
    ERP_r(2,:)=mean(ERP(central_inds,:),1);
    ERP_r(3,:)=mean(ERP(posterior_inds,:),1);
    ERP_r(4,:)=mean(ERP(temporal_inds,:),1);
    
    %0 to 789
    t_inds=find(EEG.times>0 & EEG.times<=789);
    
    % measures
    % area [mean amplitude]
    m_area=mean(ERP_r(:,t_inds),2);
    rms_area=sqrt(mean(ERP_r(:,t_inds).^2,2));
    t_area=trapz(ERP_r(:,t_inds),2);
    
    % peakes
    %positive
    [peak_p_amp,peak_p_lat_ind]=max(ERP_r(:,t_inds),[],2);    
    peak_p_lat=[];
    for ai=1:size(ERP_r,1);
        peak_p_lat(ai)=EEG.times(t_inds(peak_p_lat_ind(ai)));
    end
    
    %negative
    [peak_n_amp,peak_n_lat_ind]=min(ERP_r(:,t_inds),[],2);    
    peak_n_lat=[];
    for ai=1:size(ERP_r,1);
        peak_n_lat(ai)=EEG.times(t_inds(peak_n_lat_ind(ai)));
    end
    

    %100 to 219
    t1_inds=find(EEG.times>=100 & EEG.times<=219);
    [t1_peak_p_amp,t1_peak_p_lat_ind]=max(ERP_r(:,t1_inds),[],2);    
    t1_peak_p_lat=[];
    for ai=1:size(ERP_r,1);
        t1_peak_p_lat(ai)=EEG.times(t1_inds(t1_peak_p_lat_ind(ai)));
    end
    
    %220 to 319
    t2_inds=find(EEG.times>=220 & EEG.times<=319);
    [t2_peak_n_amp,t2_peak_n_lat_ind]=min(ERP_r(:,t2_inds),[],2);    
    t2_peak_n_lat=[];
    for ai=1:size(ERP_r,1);
        t2_peak_n_lat(ai)=EEG.times(t2_inds(t2_peak_n_lat_ind(ai)));
    end

    %320 to 540
    t3_inds=find(EEG.times>=320 & EEG.times<=540);
    [t3_peak_p_amp,t3_peak_p_lat_ind]=max(ERP_r(:,t3_inds),[],2);    
    t3_peak_p_lat=[];
    for ai=1:size(ERP_r,1);
        t3_peak_p_lat(ai)=EEG.times(t3_inds(t3_peak_p_lat_ind(ai)));
    end
    
    %541 to 789
    t4_inds=find(EEG.times>=541 & EEG.times<=789);
    [t4_peak_p_amp,t4_peak_p_lat_ind]=max(ERP_r(:,t4_inds),[],2);    
    t4_peak_p_lat=[];
    for ai=1:size(ERP_r,1);
        t4_peak_p_lat(ai)=EEG.times(t4_inds(t4_peak_p_lat_ind(ai)));
    end

    
    % print current input file measures
    fprintf(fid,'%6.3f,',[m_area',rms_area',t_area',peak_p_amp',peak_p_lat,peak_n_amp',peak_n_lat, ...
                            t1_peak_p_amp',t1_peak_p_lat, ...
                            t2_peak_n_amp',t2_peak_n_lat, ...
                            t3_peak_p_amp',t3_peak_p_lat, ...
                            t4_peak_p_amp',t4_peak_p_lat,]);
    
end