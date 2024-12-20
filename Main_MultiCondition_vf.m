%% Main_Algorithms
% Programme main visant à détecter et écrire les événements de la marche
% sur les fichiers C3D d'une session
%
% En lien avec article :
%
% Nécessaire pour faire tourner :
% 1. Matlab + btk installé + path [http://biomechanical-toolkit.github.io/]
% 2. Fichiers C3D with marker trajectories cleaned and labelled.
%
% Sorties :
% 1. Fichiers C3D avec événements détectés
% 2. Informations sur marker décisif en Foot Strike et Foot Off
% NB: Foot markers must be labelled as in MarkerList, based on pyCGM2.5 (Leboeuf et al. 2019)
% tic
clear all;
close all;
clc;

% Variables modifiables du programme :
MarkerList = {'ANK','HEE','TOE','FMH','SMH','VMH','HAX'};

%% Paramètres de calage de la routine
% frequence de coupure du filtre passe bas
fc=9;
% seuil de vitesse Vx
Seuil_Vx = 18 / 100 ;
% seuil de vitesse Vy
Seuil_Vy = 31 / 100 ;
% seuil de vitesse Vz
Seuil_Vz = 8 / 100 ;

%% Paramètres du Fenetrage
% Le fenetrage découpe le c3d par les phases oscillantes
Coeff_MPH = 0.5; MPD = 50; % MinPeakHeight = coeff * max(Vy) / MinPeakDistance

%% Initialization

% Path to C3D files
path_c3d = strcat(pwd,'\c3d');
addpath(path_c3d);

% Selection of files (.c3d)  ( /!\ order in which files are selected)
C3D_files = dir(fullfile(path_c3d, '**', '*.c3d'));


%% Events detection

ComptANKs = 0; ComptHEEs = 0; ComptTOEs = 0; ComptFMHs = 0; ComptSMHs = 0; ComptVMHs = 0; ComptHAXs = 0;
ComptANKo = 0; ComptHEEo = 0; ComptTOEo = 0; ComptFMHo = 0; ComptSMHo = 0; ComptVMHo = 0; ComptHAXo = 0;

for ii = 1 : length(C3D_files)
    
    %% Get data from file with btk
    
    Filename = C3D_files(ii).name;
    clear Markers Crop f n B A
    
    btkData = btkReadAcquisition(strcat(C3D_files(ii).folder,'\',C3D_files(ii).name));
    btkClearEvents(btkData);
    Markers = btkGetMarkers(btkData);
    Crop    = btkGetFirstFrame(btkData) - 1     ; % utile si le c3d ne commence pas à frame #1 = 1
    f = btkGetPointFrequency(btkData);
    n = btkGetPointFrameNumber(btkData);
    
    
    % définition du filtre passe bas (butterworth, ordre 4, freq coupure
    % définit en début de routine
    [B,A] = butter(4,fc/(f/2),'low');
    
    %% Coordinate system correction based on walking direction
    [Markers.SACR] = (Markers.LPSI + Markers.RPSI)./2;
    
    %% Events for both legs
    for jj = 1 : 2                      % Separation by sides
        if jj == 1
            Leg = 'L';
            WriteLeg = 'Left';
        else
            Leg = 'R';
            WriteLeg = 'Right';
        end
        
        for m = 1:length(MarkerList)
            SideMarkerList(m) = strcat(Leg, MarkerList(m));
        end
        
        %% Event Detection
        
        FramesMarkFS = [];
        FramesMarkFO = [];
        FenetrageEssai = [];
        
        for mark = 1:length(SideMarkerList)
            MarkerUsed = char(SideMarkerList(mark));
            
            
            %% Filtering
            FiltMarker = [];
            FiltMarker(:,:,1) = filtfilt(B, A, Markers.(MarkerUsed));
            
            X_Mark = FiltMarker(:,1);
            Y_Mark = FiltMarker(:,2);
            Z_Mark = FiltMarker(:,3);
            
            %% Kinematiks
            
            % Velocities
            Vx= zeros(1,n-1); Vy= zeros(1,n-1);  Vz= zeros(1,n-1);
            
            Vx = diff(X_Mark)' / (1/f);
            Vy = diff(Y_Mark)' / (1/f);
            Vz = diff(Z_Mark)' / (1/f);
            %
            
            %% Windowing
            
            [~, FramesWin] = findpeaks(abs(Vy),'MinPeakHeight',Coeff_MPH*max(abs(Vy)),'MinPeakDistance',MPD);   % Fenetrage global via Vitesse Antero-Post
            if FramesWin(1) > 4
                FramesWin = [1, FramesWin];
            end
            if FramesWin(end) < length(Vy)-4
                FramesWin = [FramesWin, length(Vy)];
            end
            
            if strcmp(MarkerUsed(2:4),MarkerList{1})
                FenetrageEssai = FramesWin;
            end
            
            %% Events calculation
            
            CrossX = zeros(1,n); CrossY = zeros(1,n); CrossZ = zeros(1,n);
            Sigma = zeros(1,n-1); Bin = zeros(1,n-1);
            FramesFS = []; FramesFO = [];
            
            for i = 1 : length(FenetrageEssai)-1
                for j = (FenetrageEssai(i)+1) : (FenetrageEssai(i+1)-2)
                    if abs(Vx(j)) <  abs((Markers.SACR(end,2) - Markers.SACR(1,2)))/length(Markers.SACR)*f * Seuil_Vx
                        CrossX(j-1)   = 1;
                        CrossX(j)     = 1;
                        CrossX(j+1)   = 1;
                        CrossX(j+2)   = 1;
                    else
                        if CrossX(j) ~= 1
                            CrossX(j) = 0;
                        end
                    end
                    
                    if abs(Vy(j)) < abs((Markers.SACR(end,2) - Markers.SACR(1,2)))/length(Markers.SACR)*f * Seuil_Vy
                        CrossY(j-1)   = 1;
                        CrossY(j)     = 1;
                        CrossY(j+1)   = 1;
                        CrossY(j+2)   = 1;
                    else
                        if CrossY(j) ~= 1
                            CrossY(j) = 0;
                        end
                    end
                    
                    if abs(Vz(j)) < abs((Markers.SACR(end,2) - Markers.SACR(1,2)))/length(Markers.SACR)*f * Seuil_Vz
                        CrossZ(j-1)   = 1;
                        CrossZ(j)     = 1;
                        CrossZ(j+1)   = 1;
                        CrossZ(j+2)   = 1;
                    else
                        if CrossZ(j) ~= 1
                            CrossZ(j) = 0;
                        end
                    end
                end
                
                Recouvrement = [];
                
                for k = FenetrageEssai(i) : FenetrageEssai(i+1)-1
                    Sigma(k) = CrossX(k) + CrossY(k) + CrossZ(k);
                    if (Sigma(k) == 3)
                        Recouvrement(end+1) = k;
                        Bin(k) = 1;
                    else
                        Bin(k) = 0;
                    end
                    
                end
                
                if ~isempty(Recouvrement)
                    FramesFS(i) = Recouvrement(1);
                    FramesFO(i) = Recouvrement(end);
                else
                    FramesFS(i) = NaN;
                    FramesFO(i) = NaN;
                end
                
            end
            
            %
            %% cleaning c3d: à cause de fenetrage, les premières et dernières frames remplissent les conditions de FS et FO
            % et sont des faux positifs donc on les supprime
            
            % FS
            while size(FramesMarkFS,1) ~= length(FramesFS)
                zz = 0;
                if (size(FramesMarkFS,1) == 0)
                    break
                elseif size(FramesMarkFS,1) > length(FramesFS)
                    for zz = 1 : length(FramesFS)
                        if ~isnan(FramesFS(zz))
                            break
                        end
                    end
                    if (FenetrageEssai(zz) < FramesFS(zz)) && (FramesFS(zz) < FenetrageEssai(zz+1))
                        FramesFS(end+1) = NaN;
                    else
                        FramesFS = [NaN FramesFS];
                    end
                elseif size(FramesMarkFS,1) < length(FramesFS)
                    for zz = 1 : length(FramesFS)
                        if ~isnan(FramesFS(zz))
                            break
                        end
                    end
                    if (FenetrageEssai(zz) < FramesFS(zz)) && (FramesFS(zz) < FenetrageEssai(zz+1))
                        FramesFS = FramesFS(1:end-1);
                    else
                        FramesFS = FramesFS(2:end);
                    end
                end
            end
            
            FramesMarkFS(:,mark) = FramesFS;
            
            % FO
            while size(FramesMarkFO,1) ~= length(FramesFO)
                zz = 0;
                if (size(FramesMarkFO,1) == 0)
                    break
                elseif size(FramesMarkFO,1) > length(FramesFO)
                    for zz = 1 : length(FramesFO)
                        if ~isnan(FramesFO(zz))
                            break
                        end
                    end
                    if (FenetrageEssai(zz) < FramesFO(zz)) && (FramesFO(zz) < FenetrageEssai(zz+1))
                        FramesFO(end+1) = NaN;
                    else
                        FramesFO = [NaN FramesFO];
                    end
                elseif size(FramesMarkFO,1) < length(FramesFO)
                    for zz = 1 : length(FramesFO)
                        if ~isnan(FramesFO(zz))
                            break
                        end
                    end
                    if (FenetrageEssai(zz) < FramesFO(zz)) && (FramesFO(zz) < FenetrageEssai(zz+1))
                        FramesFO = FramesFO(1:end-1);
                    else
                        FramesFO = FramesFO(2:end);
                    end
                end
            end
            
            FramesMarkFO(:,mark) = FramesFO;
            
            
        end
        
        
        if min(FramesMarkFS(1,:)) < 4
            FramesMarkFS = FramesMarkFS(2:end,:) ;
        end
        if max(FramesMarkFO(end,:)) > length(Vy)-4
            FramesMarkFO = FramesMarkFO(1:end-1,:);
        end
        
        %% Choix du bon marker et comptage
        WriteFSFrames = [];
        for winS = 1 : size(FramesMarkFS,1)
            [WriteFSFrames(winS), idxS] = min(FramesMarkFS(winS,:));
            
            switch idxS
                case 1
                    if ~isnan(WriteFSFrames(winS))
                        ComptANKs = ComptANKs + 1;
                    end
                case 2
                    if ~isnan(WriteFSFrames(winS))
                        ComptHEEs = ComptHEEs + 1;
                    end
                case 3
                    if ~isnan(WriteFSFrames(winS))
                        ComptTOEs = ComptTOEs + 1;
                    end
                case 4
                    if ~isnan(WriteFSFrames(winS))
                        ComptFMHs = ComptFMHs + 1;
                    end
                case 5
                    if ~isnan(WriteFSFrames(winS))
                        ComptSMHs = ComptSMHs + 1;
                    end
                case 6
                    if ~isnan(WriteFSFrames(winS))
                        ComptVMHs = ComptVMHs + 1;
                    end
                case 7
                    if ~isnan(WriteFSFrames(winS))
                        ComptHAXs = ComptHAXs + 1;
                    end
            end
        end
        
        WriteFOFrames = [];
        for winO = 1 : size(FramesMarkFO,1)
            [WriteFOFrames(winO), idxO] = max(FramesMarkFO(winO,:));
            
            switch idxO
                case 1
                    if ~isnan(WriteFOFrames(winO))
                        ComptANKo = ComptANKo + 1;
                    end
                case 2
                    if ~isnan(WriteFOFrames(winO))
                        ComptHEEo = ComptHEEo + 1;
                    end
                case 3
                    if ~isnan(WriteFOFrames(winO))
                        ComptTOEo = ComptTOEo + 1;
                    end
                case 4
                    if ~isnan(WriteFOFrames(winO))
                        ComptFMHo = ComptFMHo + 1;
                    end
                case 5
                    if ~isnan(WriteFOFrames(winO))
                        ComptSMHo = ComptSMHo + 1;
                    end
                case 6
                    if ~isnan(WriteFOFrames(winO))
                        ComptVMHo = ComptVMHo + 1;
                    end
                case 7
                    if ~isnan(WriteFOFrames(winO))
                        ComptHAXo = ComptHAXo + 1;
                    end
            end
        end
        
        %% Write Events on Files
        WriteFSFrames = WriteFSFrames(~isnan(WriteFSFrames)) ; % suppression des NaN
        WriteFOFrames = WriteFOFrames(~isnan(WriteFOFrames)) ; % suppression des NaN
        
        WriteFSFrames = WriteFSFrames - 1 + Crop ;     % FrameC3D = FrameAlgo + 1 to compensate first frame of c3d file = 1
        WriteFOFrames = WriteFOFrames - 1 + Crop + 1 ; % +1 frame coming from optimisation of 6534 iterations on 819 c3d files
        
        if jj == 1
            LHS_loc = WriteFSFrames ;
            LTO_loc = WriteFOFrames ;
        elseif jj == 2
            RHS_loc = WriteFSFrames ;
            RTO_loc = WriteFOFrames ;
        end
        
        
    end
    % FO is deleted if begins C3D since gait cycle starts with FS
    if LTO_loc(1)<RHS_loc(1)
        LTO_loc(1) = [] ;
    end
    if RTO_loc(1)<LHS_loc(1)
        RTO_loc(1) = [] ;
    end
    % FO is deleted if ends C3D
    if LTO_loc(end)>LHS_loc(end)
        LTO_loc(end) = [] ;
    end
    if RTO_loc(end)>RHS_loc(end)
        RTO_loc(end) = [] ;
    end
    
    for jj = 1:2
        if jj == 1
            for gg = 1:length(LHS_loc)
                btkAppendEvent(btkData, 'Foot Strike', LHS_loc(gg)/f, 'Left','' ,'' ,1);
            end
            for ff = 1:length(LTO_loc)
                btkAppendEvent(btkData, 'Foot Off', LTO_loc(ff)/f, 'Left','' ,'' ,2);
            end
        elseif jj == 2
            for gg = 1:length(RHS_loc)
                btkAppendEvent(btkData, 'Foot Strike', RHS_loc(gg)/f, 'Right','' ,'' ,1);
            end
            for ff = 1:length(RTO_loc)
                btkAppendEvent(btkData, 'Foot Off', RTO_loc(ff)/f, 'Right','' ,'' ,2);
            end
        end
    end
    
    btkWriteAcquisition(btkData, strcat(path_c3d,'\',Filename));
    btkCloseAcquisition(btkData)
end



%% Display

clc
sz = [1 length(MarkerList)];
MarkTypes = {'double'};
for tab = 1 : length(MarkerList)-1
    MarkTypes{end+1} = 'double';
end
TabloFS = table('Size',sz,'VariableTypes',MarkTypes,'VariableNames',MarkerList);
TabloFS.ANK = ComptANKs; TabloFS.HEE = ComptHEEs; TabloFS.TOE = ComptTOEs; TabloFS.FMH = ComptFMHs; TabloFS.SMH = ComptSMHs; TabloFS.VMH = ComptVMHs; TabloFS.HAX = ComptHAXs;
disp('Foot Strike :')
disp(TabloFS)
TabloFO = table('Size',sz,'VariableTypes',MarkTypes,'VariableNames',MarkerList);
TabloFO.ANK = ComptANKo; TabloFO.HEE = ComptHEEo; TabloFO.TOE = ComptTOEo; TabloFO.FMH = ComptFMHo; TabloFO.SMH = ComptSMHo; TabloFO.VMH = ComptVMHo; TabloFO.HAX = ComptHAXo;
disp('Foot Off :')
disp(TabloFO)
