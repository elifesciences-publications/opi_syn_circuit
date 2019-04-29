% Script to read Axograph X generated data (.axgx) and affiliated
% metadata(.csv), extracted data is appended to dataset (user-defined).
% Axograph X: https://axograph.com/support/contact
%
%
% Author:Bart C Jongbloets
% Last version: 2019-04-24
% Tested on Matlab: r2014a - r2018b
% Contact: b.c.jongbloets@gmail.com
%
% Dependencies:
% 1) bcj_XPSCParam: extraction of waveform parameters.
% 2) bcj_slidingPeak: determine peak amplitude using a sliding window
% averaging.
% 3) bcj_struct2csv: saveing dataset as a csv.
% 
% Requirements:
% 1) Folder containing .AXGX files and related metadata .csv files.
% .csv files should have a defined organization see:
% axgx_metadata_template.xlsx (First sheet need to be saved as csv, second
% sheet contains explanation of columns and rows in the template.)
%
% 2) filename of AXGX and CSV files need to be the same!
%
% 3) an dataset should be present. If not, make a empty .csv containing the
% following columns in the first row of the .csv 
% (tab-separated below are unqiue column names):
% "cellID	animalID	genoType	cellType	recordLayer	recordRegion
% stimSource	stimChannel	conditionName	signal	samplerate	onset
% risetime	slope	amplitude	onset2peakTime	halfwidth	decay	
% slope2ms	slope3ms	area	test1DrugName	waveformSequence"
%
% Script will save the old dataset with the appended new data in 'rawdataset.csv'
% with the same path as the old dataset was stored in. Beaware that it will
% overwrite the old dataset, both .mat and .csv! The latter only if .csv
% has the same name.


%locate dataset to append the AXGX data to
[fileNameDf, pathNameDf]=uigetfile('.mat','Please select existing dataset to append data to :');
load(strcat(pathNameDf,fileNameDf))
rowsDf = size(dfEla, 2)

%locate main/parent folder with all unprocessed axgx files and related .csv
[pathNameAxgx] = uigetdir();
fileList = dir(pathNameAxgx);
fileList(1:3) = [];
indexFolder = find([fileList.isdir]==1);
%iterate through all the folders and read each axgx and csv files
for folderN = 1: length(indexFolder)
	pathAXGX = strcat(pathNameAxgx,'/',fileList(folderN).name,'/');
	fileListAx = dir(pathAXGX);
	fileListAx([1:3])= [];
	cd(pathAXGX)
	nameAXGX = nan(1,size(fileListAx,1));
	for fileListR = 1:size(fileListAx,1)
		nameAXGX(1,fileListR) = ~isempty(regexp(fileListAx(fileListR).name, '.axgx'));
		if nameAXGX(1,fileListR) == 1
			locAXGX = regexp(fileListAx(fileListR).name,'.axgx');
			fileName = fileListAx(fileListR).name(1:locAXGX-1)
			csvFile = csv2cell(strcat(fileName,'.csv'));
			[data, time, meta] = importaxo(strcat(fileName,'.axgx'));
			sampleRate = floor(1/meta.SampInt);
			idx = cellfun(@(x) strcmp(x,'NaN'), csvFile);
			csvFile(idx) = {NaN};
			header = {csvFile{:,1}};
			metaStruc = cell2struct(csvFile',header,2);
			metaStruc(1) = []; % clean up the first row, which contains the header.

% In cases where the user has been copy/pasting the metadata defaults it
% may happen that a typo is carried-over through a whole list of csv files.
% In that case, uncomment the beneath code and modify it such that it can
% catch and replace typos in a consice way. However the correct csv file is
% not saved using the code below!!
% 			% In case of a typo in .csv's. Turn 'ESPC' into 'EPSC'
% 			for rowN= 1:size(metaStruc,1)
% 				if ischar(metaStruc(rowN).signalName)
% 					if strcmp(metaStruc(rowN).signalName,'ESPC')
% 						metaStruc(rowN).signalName = 'EPSC';
% 					end
% 				end
% 			end

            % extrect the waveform IDs
			waveformIDs = [metaStruc(:).dataColumn];
			waveformIDs = waveformIDs(isfinite(waveformIDs));
			
			for i = 1:length(waveformIDs)
				% per waveform check if there is a duplication required
                if ischar(metaStruc(i+1).conditionNameDupl)
                    dupl = 2;
                else
                    dupl = 1;
                end
                for dupli = 1:dupl
                    % test how many stims there were
                    stimNumber = ischar(metaStruc(i+1).stim1Type)+ischar(metaStruc(i+1).stim2Type)+ischar(metaStruc(i+1).stim3Type);
                    for stimIt = 1:stimNumber
                        cL = size(dfEla,2)+1;
                        dfEla(cL).cellID = strcat(metaStruc(1).animalID,metaStruc(1).cellID);
                        dfEla(cL).animalID = metaStruc(1).animalID;
                        dfEla(cL).genoType = metaStruc(1).geno;
                        dfEla(cL).cellType = metaStruc(1).cellType;
                        dfEla(cL).recordLayer = metaStruc(1).recordLayer;
                        dfEla(cL).recordRegion = metaStruc(1).recordRegion;
                        if dupli == 1
                            dfEla(cL).conditionName = metaStruc(i+1).conditionName;
                        elseif dupli == 2
                            dfEla(cL).conditionName = metaStruc(i+1).conditionNameDupl;
                        end
                        dfEla(cL).signal = metaStruc(i+1).signalName;
                        dfEla(cL).samplerate = sampleRate;
                        dfEla(cL).averageTrace = data{waveformIDs(i)-1,1}';
                        dfEla(cL).test1DrugName = metaStruc(i+1).test1DrugName;
                        dfEla(cL).waveformSequence = i;
                        dfEla(cL).allTraces=  NaN;
                        dfEla(cL).meta = metaStruc;
                        
                        if stimIt == 1
                            if ~isempty(regexp(metaStruc(i+1).stim1Type, 'Opto'))
                                dfEla(cL).stimSource = metaStruc(1).virus1Location;
                                dfEla(cL).stimChannel =metaStruc(1).virus1Construct;
                            elseif ~isempty(regexp(metaStruc(i+1).stim1Type, 'Elec'))
                                dfEla(cL).stimSource = metaStruc(i+1).stim1Type;
                                dfEla(cL).stimChannel =metaStruc(i+1).stim1Type;
                            elseif ~isempty(regexp(metaStruc(i+1).stim1Type, 'Poly'))
                                dfEla(cL).stimSource = metaStruc(i+1).stim1Type;
                                dfEla(cL).stimChannel =metaStruc(i+1).stim1Type;
                            end
                            dfEla(cL).stimOnset = metaStruc(i+1).stim1Onset;
                        elseif stimIt == 2
                            if ~isempty(regexp(metaStruc(i+1).stim2Type, 'Opto'))
                                dfEla(cL).stimSource = metaStruc(1).virus1Location;
                                dfEla(cL).stimChannel =metaStruc(1).virus1Construct;
                            elseif ~isempty(regexp(metaStruc(i+1).stim2Type, 'Elec'))
                                dfEla(cL).stimSource = metaStruc(i+1).stim2Type;
                                dfEla(cL).stimChannel =metaStruc(i+1).stim2Type;
                            elseif ~isempty(regexp(metaStruc(i+1).stim2Type, 'Poly'))
                                dfEla(cL).stimSource = metaStruc(i+1).stim2Type;
                                dfEla(cL).stimChannel =metaStruc(i+1).stim2Type;
                            end
                            dfEla(cL).stimOnset = metaStruc(i+1).stim2Onset;
                        end
                        % extract all parameters from the traces. Before doing so the waveform needs to be normalized to a 50 ms baseline prior to stimonset
                        windowWidth = 0.05; %seconds
                        stimOnset = dfEla(cL).stimOnset;
                        stopWindow = ceil(double((stimOnset/1000)*sampleRate)); % datapoints
                        startWindow = ceil(double((stopWindow - (windowWidth*sampleRate)))); % datapoints
                        avTrace = dfEla(cL).averageTrace;
                        normTrace = avTrace - mean(avTrace(startWindow:stopWindow));
                        
                        % incase of electrical stimulation the stimulus artifact needs to be cut off (set to baseline level = 0 pA)
                        if ~isempty(regexp(dfEla(cL).stimSource, 'Elec'))
                            normTrace(stopWindow:stopWindow+(0.002*sampleRate)) = 0; 
                            % turn the first 2 ms of the stimulation to 0 pA
                        end
                        % now we can run the parameter extraction function
                        saveFigName = strcat(dfEla(cL).cellID,dfEla(cL).stimSource,dfEla(cL).cellType,dfEla(cL).signal);
                        [Eout] = bcj_XPSCParam(normTrace, [0 3 50], [stimOnset 1], sampleRate, 0, saveFigName);
                        artiLength = 0.002;
                        if ~isempty(regexp(dfEla(cL).signal,'EPSC'))
                            checkArti = 1;
                            while checkArti == 1
                                if isnan(Eout.amplitude(1))
                                    checkArti = 0;
                                elseif Eout.amplitude(1)*10^12 <0
                                    checkArti = 0;
                                elseif artiLength == 0.005
                                    Eout.onset = NaN;
                                    Eout.risetime = NaN;
                                    Eout.slope = NaN;
                                    Eout.amplitude = [NaN NaN];
                                    Eout.halfwidth = NaN;
                                    Eout.decay = NaN;
                                    Eout.slope2ms = NaN;
                                    Eout.slope3ms = NaN;
                                    Eout.area = NaN;
                                    checkArti = 0;
                                elseif Eout.amplitude(1)*10^12 > 0
                                    % in this case the artifact is obscuring the amplitude detection
                                    normTrace(stopWindow:stopWindow+((artiLength+0.001)*sampleRate))=0; % turn the first 3 ms of the stimulation to 0 pA
                                    [Eout] = bcj_XPSCParam(normTrace, [0 3 50], [stimOnset 1], sampleRate, 0, saveFigName);
                                end
                            end
                        end
                        dfEla(cL).onset = Eout.onset(1);
                        dfEla(cL).risetime = Eout.risetime;
                        dfEla(cL).slope = Eout.slope;
                        dfEla(cL).amplitude = Eout.amplitude(1)*10^12; % convert to pA
                        try
                            dfEla(cL).onset2peakTime = ((Eout.amplitude(2)/dfEla(cL).samplerate)*1000)-dfEla(cL).onset; 
                            % input Eout.amplitude(2) is relative to start of pulse. So correct for samplerate and calculate relative to onset time: onset - peak time
                        catch
                            dfEla(cL).onset2peakTime = NaN;
                        end
                        dfEla(cL).halfwidth = Eout.halfwidth;
                        dfEla(cL).decay = Eout.decay;
                        dfEla(cL).slope2ms = Eout.slope2ms;
                        dfEla(cL).slope3ms = Eout.slope3ms;
                        dfEla(cL).area = Eout.area*10^12; % convert to pA/s
                    end
                end
			end
		end
	end
end
%Fill all empty spots with NaN
for i = 1: size(dfEla,2)
    if isempty(dfEla(i).test1DrugName)
        dfEla(i).test1DrugName = NaN;
    end
    if isempty(dfEla(i).waveformSequence)
        dfEla(i).waveformSequence = NaN;
    end
    if isempty(dfEla(i).meta)
        dfEla(i).meta = NaN;
    end
    if isempty(dfEla(i).stimOnset)
        dfEla(i).stimOnset = NaN;
    end
end
dfExport = rmfield(dfEla, {'averageTrace','allTraces','meta','stimOnset'});

bcj_struct2csv(dfExport,strcat(pathNameDf,'rawdataset.csv'))
save(strcat(pathNameDf,fileNameDf),'dfEla')


