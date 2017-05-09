function beh_model_experiment(sj_id, run_id)

% This function presents the gridworld search paradigm.
%
%   Inputs
%           sj_id  : participant ID, string
%           run_id : run ID, scalar
%
%   Outputs
%           None, saves results file to disc
%
% Copyright (C) Lilla Horvath, Dirk Ostwald
% -------------------------------------------------------------------------

% -------------------------------------------------------------------------
% ---------------------------- Initialization -----------------------------
% -------------------------------------------------------------------------
clc
close all

% reset the state of the random nubmer generator based on computer clock
rng('shuffle');

% -------------------------------------------------------------------------
% ------------------------ Paradigm Parameters ----------------------------
% -------------------------------------------------------------------------
% experimental parameters
dim             = 5                                                         ; % grid dimension 
n_tgt           = 2                                                         ; % number of targets
stim_dir        = [pwd '\Stimuli']                                          ; % stimulus directory   
res_dir         = [pwd '\Results\' sj_id]                                   ; % participant specific result directory
min_p           = 50                                                        ; % 
b_max           = 3                                                         ; % maximal number of blocks per task
p_lim           = [ 0.1  0.2  0.4  0.2  0.1]                                ; % +- two steps with probabilities given by p
poss_lim        = [-2 -1  0  1  2]                                          ; % steps with probabilities given by p

% short vs long stimuli presentation dependent parameters of the experiment
% between subject condition

bsc = 0                                                                     ;  % between subject condition, 0 corresponds to the long and 1 to the short condition

if bsc == 0
    
    % experimental parameter
    n_task          = 4                         ; % number of tasks per run
    % timing parameters
    l_s             = 3                         ; % lower endpoint for the short duration (state and obs pres)  
    u_s             = 5                         ; % upper endpoint for the short duration (state and obs pres)  
    l_l             = 6                         ; % lower endpoint for the long duration  (fixation point)      
    u_l             = 8                         ; % upper endpoint for the long duration  (fixation point)      
    
elseif bsc == 1
    
    n_task          = 8                         ; % number of tasks per run
    % timing parameters
    l_s             = 1.5                       ; % lower endpoint for the short duration (state and obs pres)  
    u_s             = 2.5                       ; % upper endpoint for the short duration (state and obs pres)  
    l_l             = 3                         ; % lower endpoint for the long duration  (fixation point)      
    u_l             = 4                         ; % upper endpoint for the long duration  (fixation point)      
    
end
    
% cue display coordinates (up, down, left, right) 
cue_coords      = [[0 177];[0 -177]; [-177 0];[177 0]];
        
% cogent display configuration parameters
cog_ws          = 1                         ; % window size
cog_rs          = 3                         ; % display resolution (mapped with laptop screen)
cog_bc          = [0 0 0]                   ; % background color 
cog_fc          = [1 1 1]                   ; % foreground color
cog_fn          = 'Helvetia'                ; % fontname 
cog_fs          = 45                        ; % fontsize
cog_bu          = 11                        ; % number of offscreen-buffers
cog_bd          =  0                        ; % bits per pixel

% cogent keyboard configuration parameters
cog_ql          = 100                       ; % quelength
cog_kr          = 5                         ; % keyboard resolution
cog_km          = 'nonexclusive'            ; % keyboard mode

% scanner or behavioural testing flag
isscan          = 1                        ; % 0: behavioural lab, 1: MR lab

% response keys differentiation
if isscan
    
    % button box "left" keys
    % ---------------------------------------------------------------------
    % blue  , yellow, green , red
    % e = 5, w = 23 , n = 14, d =  4
    % use cursor keys
    k_left      = 5;
    k_up        = 23;
    k_down      = 14;
    k_right     = 4;
    
else
    
    % use cursor keys
    k_left      = 97;
    k_up        = 95;
    k_down      = 100;
    k_right     = 98;
    
end

% -------------------------------------------------------------------------
% -------------------------- Stimulus Loading -----------------------------
% -------------------------------------------------------------------------
% grid background stimulus filenames
grid_fname  =   {  'upper_corner_left_s.jpg'    , ...
                    'upper_corner_right_s.jpg'  , ...
                    'lower_corner_left_s.jpg'   , ...
                    'lower_corner_right_s.jpg'  , ...
                    'side_upper_s.jpg'          , ...
                    'side_left_s.jpg'           , ...
                    'side_right_s.jpg'          , ...
                    'side_lower_s.jpg'          , ...
                    'middle_s.jpg'                   };

% cue stimulus filenames
cue_fname   =   {   'dark_hor_s.jpg'            ,...
                    'light_hor_s.jpg'           ,...
                    'dark_ver_s.jpg'            ,...
                    'light_ver_s.jpg'                };
                
% decision prompt filenames
arrow_fname =   {   'up_s.jpg'                  ,...
                    'down_s.jpg'                ,...
                    'left_s.jpg'                ,...                    
                    'right_s.jpg'                   };
                
                
% load grid stimuli into workspace and normalize to cogent image format
grid_pict      = cell(1,numel(grid_fname));
for i = 1:numel(grid_fname)
    grid_pict{i} = double(imread(fullfile(stim_dir,grid_fname{i})))./255;
end
                                               
% load cue stimuli into workspace and normalize to cogent image format
cue_pict      = cell(1,numel(cue_fname));
for i = 1:numel(cue_fname)
    cue_pict{i} = double(imread(fullfile(stim_dir,cue_fname{i})))./255;
end

% load arrow stimuli into workspace and normalize to cogent image format
arrow_pict      = cell(1,numel(cue_fname));
for i = 1:numel(cue_fname)
    arrow_pict{i} = double(imread(fullfile(stim_dir,arrow_fname{i})))./255;
end

% load target stimulus into workspace and normalize to cogent image format
target = double(imread(fullfile(stim_dir,'treasure2_s.jpg')))./255;

% -------------------------------------------------------------------------
% ----------------------- Paradigm Presentation ---------------------------
% -------------------------------------------------------------------------
% initialize cogent
config_display(cog_ws, cog_rs, cog_bc, cog_fc, cog_fn, cog_fs, cog_bu, cog_bd);
config_keyboard(cog_ql, cog_kr, cog_km);
start_cogent;

% display starting screen
clearpict(1);
preparestring('Ready',1)
drawpict(1);  

% differentiate presentation start conditions 
% -------------------------------------------------------------------------
if isscan
    
    % wait for and read scanner trigger 
    outportb(890,32) 
    startState  = inportb(888);
    oldValue    = startState;
    triggerNum  = 0;
    
    while triggerNum < 4
    
        val = inportb(888);
        if val ~= oldValue
            triggerNum = triggerNum + 1;
            triggertime(triggerNum) = time;
        end
    
        oldValue = val;
    
    end
else
    wait(2000);
end

% -------------------------------------------------------------------------   
% ------------------------ Initial Fixation Block -------------------------
% -------------------------------------------------------------------------
% clear buffer 2 to default background color
clearpict(2);

% write fixation cross into buffer 2
preparestring('+',2);

% present buffer 2, t_expstart is the time of the experiment start, i.e. 
% the scanner trigger synchronized zero time point.
t_expstart  = drawpict(2); 

% sample fixation duration and wait
dur_fp_i = unifrnd(l_l,u_l)*1000; 
wait(dur_fp_i);

% -------------------------------------------------------------------------
% --------------------------- Cycle Over Tasks ----------------------------
% -------------------------------------------------------------------------
% initialize global block counter over tasks
blockcount = 1;

% cycle over tasks
for task = 1:n_task
    
    % sample target position from independent uniform discrete distribution
    % ---------------------------------------------------------------------
    % initialize target coordinate array
    tgt = ones(n_tgt,2);

    while tgt(1,:) == 1 || tgt(2,:) == 1
        
        for j = 1 : n_tgt
            tgt(j,:) = unidrnd(5,1,2);
        end
    
        while isequal(tgt(1,:),tgt(2,:))
            
            for i = n_tgt
                tgt(i,:) = unidrnd(5,1,2);
            end
            
        end
    end

    % Observation probability definition
    % ---------------------------------------------------------------------
    % initialize distance and probability arrays
    l1_map_t    = NaN(dim, dim, n_tgt);
    p_map_t     = NaN(dim,dim,n_tgt);
    max_l1      = NaN(1,n_tgt);
    lin_acc_g   = NaN(1,n_tgt);
        
    % cycle over targets
    for t = 1:n_tgt
    
        % cycle over matrix rows
        for i = 1:dim

            % cyle over matrix columns
            for j = 1:dim
                % evaluate l1 distance to the target for each grid cell and each target
                l1_map_t(i,j,t) = pdist([i,j;tgt(t,1),tgt(t,2)], 'cityblock'); % parametrize 1 and 2 later!
            end
        end

        % maximum l1 distance for each target
        max_l1(t)          = max(max(l1_map_t(:,:,t))); 

        % linear increase in cue accuracy 
        lin_acc_g(t)       = min_p/(max_l1(t) - 1); 

        % cue accuracy probability (right direction)
        p_map_t(:,:,t)           = (100-((l1_map_t(:,:,t) - 1)*lin_acc_g(t)))/100;

    end

    % join probability maps
    p_12 = NaN(dim,dim);
    for i = 1:dim
        for j = 1:dim
            p_12(i,j) =  join_p(p_map_t(i,j,1),p_map_t(i,j,2));
        end
    end

    % Optimal number of steps for the current problem
    % ---------------------------------------------------------------------
    % evaluate the optimal path according to Dijkstra's algorithm
    % CAVE: OPTIMAL FORAY USES TGT IN TRANSPOSED FORM - HARMONIZE!!!
    optglobalpath   = optimal_foray(dim,tgt');
    
    % evaluate the number optimally visited nodes, including the start node
    n_opt           = size(optglobalpath,2); 
    
    % evaluate the number of optimal performed steps
    s_opt           = n_opt - 1;

    % ---------------------------------------------------------------------
    % ------------------------ Cycle over Task Blocks ---------------------
    % ---------------------------------------------------------------------
    for block = 1:b_max
        
        % reset agent position to row 1, column 1 of grid matrix
        pos = [1 1];  
        
        step_lim = zeros(1,1); % initialize step_lim array
        
        while step_lim < 1 % sample from poss_lim until the step_limit (poss_lim + s_opt) is at least 1
        
            % evaluate the step number limit
            step_lim =  s_opt + poss_lim(find(logical(mnrnd(1,p_lim))));
            
        end
        
        % initialize block log for current block
        blocklog = NaN(step_lim+1,34);
        
        % initialize target found flag
        % -----------------------------------------------------------------
        % [0 0] : no target found
        % [1 0] : target at tgt_coord(1,:) found
        % [0 1] : target at tgt_coord(2,:) found
        % [1 1] : target at tgt_coord(1,:) and tgt_coord(2,:) found
        % -----------------------------------------------------------------
        found_t = zeros(1,n_tgt);
        f_within_tl = zeros(1,1);
    
        % -----------------------------------------------------------------
        % ------------------ Cycle over Block Trials ----------------------
        % -----------------------------------------------------------------
        for trial = 1:step_lim
            
            % record agent position prior to move
            pos_trial = pos;
 
            % Grid Cell State Preparation
            % -------------------------------------------------------------
            % clear buffer 3 to default background color
            clearpict(3);
                       
            % write the respective background grid image into buffer 3
            preparepict(grid_pict{get_pictidx(pos, dim)},3);
                        
            % additionally write the target to buffer 3, if appropriate,
            for t = 1:n_tgt
                
                % show image only if the target has not been found previously
                if isequal(pos,tgt(t,:)) && found_t(t) == 0
                    
                    % write targt into buffer 3
                    preparepict(target,3,-70,-70);
                    
                    % set target found flag to 1
                    found_t(t) = 1;
                end
            end

            % additionally write fixation cross into buffer 3 
            preparestring('+',3);
                        
            % additionally write string of current grid index into buffer 3
            preparestring(['(',num2str(pos(1)),' , ', num2str(pos(2)),')'],3,70,-70);
            
            % additionally evaluate and write string of number of trials
            % left and targets found into buffer 3
            preparestring([num2str(step_lim - trial + 1), ' (' num2str(sum(found_t)),')'],3,-190,190);
            
            % Grid Cell State Presentation 
            % -------------------------------------------------------------
            % present the grid cell state
            t_s         = drawpict(3);
            
            % evaluate state presentation time wrt experiment start
            t_s         = t_s - t_expstart;
            
            % evaluate state presentation duration
            dur_s       = unifrnd(l_s,u_s)*1000;      
            
            % present for dur_s
            wait(dur_s)
            
                            
            % Observation Preparation
            % -------------------------------------------------------------
            % clear buffer 4 to default background color
            clearpict(4);
            
            % differentiate l1 distance and probability maps
            if ~ismember(1, found_t) % no target found so far
            
                p_map   = p_12;
                l1_map  = l1_map_t;
           
            elseif found_t(1) == 1 && found_t(2) == 0  % target 1 found, target 2 not found
                              
                 p_map = p_map_t(:,:,2);
                l1_map  = l1_map_t(:,:,2);
            
            elseif found_t(1) == 0 && found_t(2) == 1  % target 1 not found, target 2 found
                
                p_map   = p_map_t(:,:,1);
                l1_map  = l1_map_t(:,:,1);
            
            elseif found_t(1) == 1 && found_t(2) == 1 % both targets found
                
                % end the loop over trials, go to final state presentation
                f_within_tl = 1;
                
                % block-log array for trial loop - variables of trial where 
                % goal reached before steps limit reached 
                % (state presentation only)
                % ---------------------------------------------------------
                blocklog(trial,1:2)     = tgt(1,:)      ; % (1:2)   target 1 position
                blocklog(trial,3:4)     = tgt(2,:)      ; % (3:4)   target 2 position
                blocklog(trial,5:6)     = found_t       ; % (5:6)   targets found flag 
                blocklog(trial,7:8)     = pos_trial     ; % (7:8)   current trial matrix row position
                blocklog(trial,9)       = block         ; % (9)     block number (= attempt on this task)  
                blocklog(trial,10)      = f_within_tl   ; % (10)    both targets found within trial loop flag
                blocklog(trial,11)      = t_s           ; % (11)    onset time of trial start/state presentation
                
                break
                
            end
            
            % get cue picture information
            [cueidx,d_cue,s_cue,o_cue] = get_cue(pos, l1_map, p_map);
            
            % write cue pictures into buffer 4 (up, down, left, right)
            for c = 1:length(cueidx)
                
                % if cueidx is zero, do write image into array
                if ~(cueidx(c) == 0)
                    preparepict(cue_pict{cueidx(c)},4, cue_coords(c,1),cue_coords(c,2));
                end
                
            end
               
            % additionally write fixation cross into buffer 4 
            preparestring('+',4);
            
            % additionally evaluate and write string of number of trials and targets found into buffer 4          
            preparestring([num2str(step_lim - trial + 1), ' (' num2str(sum(found_t)),')'],4,-190,190);
            
            % Observation Presentation
            % -------------------------------------------------------------
            % present the observation
            t_o = drawpict(4);       
            
            % evaluate observation presentation time wrt experiment start
            t_o = t_o - t_expstart;
            
            % evaluate observation presentation duration
            dur_o = unifrnd(l_s,u_s)*1000;
            
            % present for dur_o
            wait(dur_o);
       
            % Decision Prompt Preparation
            % -------------------------------------------------------------
            % clear buffer 5 to default background color
            clearpict(5);
 
            % WRITE GET_ARROWIDX SUBFUNCTION FOR CONSITENCY
            
            if pos(1) ~= 1 % not uppermost row
                preparepict(arrow_pict{1},5,0,90);      % arrow up
            end
            
            if pos(1) ~= dim % not lowermost row
                preparepict(arrow_pict{2},5,0,-90);     % arrow down
            end
                        
            if pos(2) ~= 1 % not leftmost column
                preparepict(arrow_pict{3},5,-90,0);      % arrow left
            end
            
            if pos(2) ~= dim % not rightmost column
                preparepict(arrow_pict{4},5,90,0);      % arrow right
            end
   
            % additionally write fixation cross into buffer 5
            preparestring('+',5);
            
            % additionally evaluate and write string of number of trials and targets found into buffer 5          
            preparestring([num2str(step_lim - trial + 1), ' (' num2str(sum(found_t)),')'],5,-190,190);

            % Decision Prompt Presentation 
            % -------------------------------------------------------------
            % present decision prompt 
            t_c = drawpict(5);
            
            % evaluate decision prompt presentation time
            t_c = t_c - t_expstart;
            
            % evaluate maximal decision prompt presentation time
            dur_r       = unifrnd(l_s,u_s)*1000;
                    
            % Read participant responses
            % -------------------------------------------------------------
            % clears all prior keyboard events  
            clearkeys;                    
            
            % read all key events since the last readkeyes
            readkeys; 
        
            % evaluate keyboard responses
            [key_d,t_d] = waitkeydown(dur_r, [k_up,k_left,k_right,k_down]);
                       
            % record keyboard response time            
            t_d         = t_d - t_expstart;
            
            % case that the participants responds faster than dur_r
            if t_d - t_c < dur_r
                
                % Post-Decision Fixation Preparation
                % ---------------------------------------------------------
                % clear buffer 6 to default background color
                clearpict(6);
               
                % additionally evaluate and write string of number of trials and targets found into buffer 6          
                preparestring([num2str(step_lim - trial), ' (' num2str(sum(found_t)),')'],6,-190,190);
                
                %  write fixation cross to buffer 6
                preparestring('+',6);    
                
                % present fixation 
                t_pd_f = drawpict(6); 
         
                % present fixation for the remaining time of the response interval 
                wait(dur_r - (t_d - t_c));    
                                
            end
                 
            % Decision Consequence - Agent movement
            % -------------------------------------------------------------
            % evaluate decision, if there is one
            if ~isempty(key_d)

                % evaluate new position based on key press 
                if key_d == k_up
                    pos(1) = pos(1) - 1;
                elseif key_d == k_down
                    pos(1) = pos(1) + 1;
                elseif key_d == k_left
                    pos(2) = pos(2) - 1;
                elseif key_d == k_right
                    pos(2) = pos(2) + 1;
                end
                    
                % make inappropriate steps impossible 
                for i = 1:length(pos)
                    if pos(i) < 1
                        pos(i) = 1 ;
                    elseif pos(i) > dim
                        pos(i) = dim;
                    end
                end

            % if there is no response the agent remains at its location    
            else
                key_d   = NaN;
                t_d     = NaN;
            end
            
            % fill in block-log array
            % -------------------------------------------------------------
            blocklog(trial,1:2)     = tgt(1,:)      ; % (1:2)   target 1 position
            blocklog(trial,3:4)     = tgt(2,:)      ; % (3:4)   target 2 position
            blocklog(trial,5:6)     = found_t       ; % (5:6)   targets found flag 
            blocklog(trial,7:8)     = pos_trial     ; % (7:8)   current trial matrix row position
            blocklog(trial,9)       = block         ; % (9)     block number (= attempt on this task)  
            blocklog(trial,10)      = f_within_tl   ; % (10)    both targets found within trial loop flag
            blocklog(trial,11)      = t_s           ; % (11)    onset time of trial start/state presentation
            blocklog(trial,12)      = t_o           ; % (12)    onset time of observation presentation
            blocklog(trial,13)      = t_c           ; % (13)    onset time of decision cue (arrows) presentation
            blocklog(trial,14)      = t_d           ; % (14)    onset time of participant's button press
            blocklog(trial,15)      = key_d         ; % (15)    participant/agent decision (NaN = no key pres)
            blocklog(trial,16)      = t_d - t_c     ; % (16)    reaction time (button press onset - decision cue onset, NaN: no key press - no reaction time)
            blocklog(trial,17:20)   = d_cue         ; % (17:20) display cue flags
            blocklog(trial,21:24)   = s_cue         ; % (21:24) sample cue flags
            blocklog(trial,25:28)   = o_cue         ; % (25:28) cue sampling outcomes
            blocklog(trial,29:32)   = cueidx        ; % (29:32) cue sampling outcomes
            
               
        end % trial loop

        % -----------------------------------------------------------------
        % ------------- Final State Presentation Preparation --------------
        % -----------------------------------------------------------------
        % clear buffer 7 to default background color
        clearpict(7);
        
        % Evaluate target presence
        % -----------------------------------------------------------------
        
        if f_within_tl == 0
        
            % write the respective background grid image into buffer 7
            preparepict(grid_pict{get_pictidx(pos, dim)},7);

            % additionally write the target to buffer 7, if appropriate,
            for t = 1:n_tgt

                % show image only if the target has not been found previously
                if isequal(pos,tgt(t,:)) && found_t(t) == 0

                    % write targt into buffer 3
                    preparepict(target,7,-70,-70);

                    % set target found flag to 1
                    found_t(t) = 1;
                end
            end
        
            % additionally write fixation cross to buffer 7
            preparestring('+',7);
        
            % additionally write string of current grid index into buffer 7
            preparestring(['(',num2str(pos(1)),' , ', num2str(pos(2)),')'],7,70,-70);

            % additionally evaluate and write string of number of trials and targets found into buffer 7          
            preparestring([num2str(step_lim - trial), ' (' num2str(sum(found_t)),')'],7,-190,190);
        
            % Final State Presentation Preparation
            % -----------------------------------------------------------------
            % present the final state 
            t_fin_s = drawpict(7); 
         
            % evaluate final state presentation time
            t_fin_s = t_fin_s - t_expstart;

            % evaluate duration of the final state presentation
            dur_fin_s = unifrnd(l_s,u_s)*1000;
      
            % present for dur_fin_s
            wait(dur_fin_s);  
            
            % fill in block-log array for final state 
            % ---------------------------------------------------------      
            blocklog(trial+1,1:2)     = tgt(1,:)      ; % (1:2)   target 1 position
            blocklog(trial+1,3:4)     = tgt(2,:)      ; % (3:4)   target 2 position
            blocklog(trial+1,5:6)     = found_t       ; % (5:6)   targets found flag 
            blocklog(trial+1,7:8)     = pos           ; % (7:8)   current trial matrix row position
            blocklog(trial+1,9)       = block         ; % (9)     block number (= attempt on this task)
            blocklog(trial+1,10)      = f_within_tl   ; % (10)    both targets found within trial loop flag (if 0: endstate pres)
            blocklog(trial+1,11)      = t_fin_s       ; % (11)    onset time of trial start/state presentation
            
            
        end
        
        % End of block - case distinction and information presentation
        %------------------------------------------------------------------
        
        % Reset Information Preparation
        % -----------------------------------------------------------------   
        clearpict(8);
        
        if sum(found_t) == 2
            preparestring('Both targets found',8);
            info_block = 1;
        else
            preparestring('Step Limit Reached',8);
            if block < b_max
                preparestring('Resetting Position',8, 0,-40);
                info_block = 2;
            else
                preparestring('Attempt Limit Reached',8,0,-40);
                info_block = 3;
            end
        end
        
        t_info = drawpict(8);
        
        % evaluate information presentation time wrt experiment start
        t_info = t_info - t_expstart;

        % evaluateinformation presentation duration
        dur_info = unifrnd(l_s,u_s)*1000;

        % present for dur_info
        wait(dur_info);
        
        clearpict(9);
        
        t_fix_b = NaN;
        dur_fix_b = NaN;
        
        if info_block == 2 && rem(blockcount,4) == 0
            preparestring('+',9)
            t_fix_b = drawpict(9);
            t_fix_b = t_fix_b - t_expstart;
            dur_fix_b = unifrnd(l_l,u_l)*1000;
            wait(dur_fix_b);
        end
            
        
        % fill in block-log array for block information
        % -----------------------------------------------------------------
        blocklog(trial+1,12)      = info_block ; % (12)    type of info at the end of the block
        blocklog(trial+1,13)      = t_info     ; % (13)    onset time of information presentation (last state only!)
        blocklog(trial+1,14)      = step_lim   ; % (14)    block specific step limit
        blocklog(trial+1,15)      = s_opt      ; % (15)    optimal number of steps for task
        blocklog(trial+1,16)      = t_fix_b    ; % (16)    onset time of fixation presentation after 4th block
        
        % allocate variable values to runlog array
        runlog{blockcount} = blocklog;
            
        % increase blockcounter
        blockcount = blockcount + 1;
        
        if block < b_max && sum(found_t) == 2
            break
        end
        
    end
    
    % End of task - case distinction and information presentation
    %------------------------------------------------------------------
    
    % Reset Information Preparation
    % -----------------------------------------------------------------
    clearpict(10);
    
    if task < n_task
        preparestring('Creating New Task', 10);
        info_task = 1;
    else 
        preparestring('Task Limit Reached',10);
        preparestring('Ending Program',10,0,-40);
        info_task = 2;
    end
    
    t_info_tl = drawpict(10);
    
    % evaluate information presentation time wrt experiment start
    t_info_tl = t_info_tl - t_expstart;

    % evaluateinformation presentation duration
    dur_info_tl = unifrnd(l_s,u_s)*1000;

    % present for dur_info
    wait(dur_info_tl);
    
    clearpict(11);
    
    if rem((blockcount-1),4) == 0
        preparestring('+',11)
        t_fix_b = drawpict(11);
        t_fix_b = t_fix_b - t_expstart;
        dur_fix_b = unifrnd(l_l,u_l)*1000;
        wait(dur_fix_b);
    end
    
     % fill in block-log array for end of task 
     % --------------------------------------------------------------------
     blocklog(trial+1,17)       = info_task     ; % (17)    type of info after a task
     blocklog(trial+1,18)       = t_info_tl     ; % (18)    onset time of information presentation - end of task
     blocklog(trial+1,19)       = dur_info_tl   ; % (19)    duration of info presentation - end of task
     blocklog(trial+1,20)       = t_fix_b       ; % (20)    onset time of fixation presentation after 4th block
     blocklog(trial+1,21)       = dur_fix_b     ; % (21)    duration of fixation presentation
       
%     % allocate variable values to runlog array
    runlog{blockcount-1} = blocklog;
    
    % task specific l1 map
    l1_maps_t{task} = l1_map_t;
    
    % task specific p map
    p_maps_t{task} = p_map_t;
    
    % task specific joint p map
    p_12_maps{task} = p_12;
    
    % task specific optimal steps
    s_opt_all{task} = s_opt;

end % task loop

% save  
% ---------------------------------------------------------
% save to disc 
save(fullfile(res_dir, [sj_id '_Run_' num2str(run_id), '.mat']), 'runlog', 'l1_maps_t', 'p_maps_t', 'p_12_maps', 's_opt_all', 'bsc')      

% end cogent
stop_cogent;

end

function [pictidx] = get_pictidx(pos_coord, dim)

% This function returns the indices of pictures in the grid_pict array
% based on the current coordinates of the agent
%
%   Inputs
%       pos_coord   : 2 x 1 array of agent row and column coordinatee
%       dim         : grid world dimension
%
%  Outpus
%       pictidx     : scalar index for grid_pict array
%
% Copyright (C) Dirk Ostwald
% -------------------------------------------------------------------------

% determine the index
% -------------------------------------------------------------------------
% upper left corner
if pos_coord(1) == 1 && pos_coord(2) == 1 
    pictidx = 1;

% upper right corner
elseif pos_coord(1) == 1 && pos_coord(2) == dim
    pictidx = 2;

% lower left corner
elseif pos_coord(1) == dim && pos_coord(2) == 1
    pictidx = 3;

% lower  right corner
elseif pos_coord(1) == dim && pos_coord(2) == dim
    pictidx = 4;

% uppermost row
elseif pos_coord(1) == 1 && (pos_coord(2) > 1 && pos_coord(2) < dim)
    pictidx = 5;

% leftmost row
elseif (pos_coord(1) > 1 && pos_coord(1) < dim) &&  pos_coord(2) == 1
    pictidx = 6;

% rightmost row
elseif (pos_coord(1) > 1 && pos_coord(1) < dim) &&  pos_coord(2) == dim
    pictidx = 7 ;       

% lowermost row
elseif pos_coord(1) == dim && (pos_coord(2) > 1 && pos_coord(2) < dim)
    pictidx = 8;     

% interior of the grid
else
    pictidx = 9;    
end % grid cell image selection
end

function [cueidx,d_cue,s_cue,o_cue] = get_cue(pos, l1_map, p_map)

% This function specifies four variables, each with two possible values (
% filenames of images), where the current value of the variable depends on 
% the function of  Manhattan-distance (l1 norm), and the outcome of 
% sampling from a Bernoulli distribution
%
% Input
%       pos         : 1 x 2 array if agent gridworld matrix coordinate 
%       l1_map      : dim x dim x n_tgt array of l1 distances to target
%       p_map       : dim x dim array of probabilities (joint probability
%                     or probability map of target1/target2 
%
% Output
%       cueidx      : 4 x 1 (up, down, left, right) array of cue image array
%                     indices with the semantic 
%                               (0) = no bar
%                               (1) = dark vertical 
%                               (2) = light vertical 
%                               (3) = dark  horizontal
%                               (4) = light horizontal
%       d_cue       : 4 x 1 (up, down, left, right) array of display cue flags
%       s_cue       : 4 x 1 (up, down, left, right) array of sample  cue flags
%       o_cue       : 4 x 1 (up, down, left, right) array of cue sample outcomes
%
% Copyright (C) Lilla Horvath, Dirk Ostwald
% -------------------------------------------------------------------------
%
% Initialization
% -------------------------------------------------------------------------
% recover problem cardinalities
n_act       = 4                     ; % largest possible action set
n_dim       = size(l1_map,1)        ; % assuming a square grid world

% evaluate possible future locations("act"ion "o"ut"c"omes)
act_oc      = repmat(pos,n_act,1)   ;
act_oc(1,1) = act_oc(1,1) - 1       ; % step up
act_oc(2,1) = act_oc(2,1) + 1       ; % step down
act_oc(3,2) = act_oc(3,2) - 1       ; % step left
act_oc(4,2) = act_oc(4,2) + 1       ; % step right

% cue display, sample, and sample outcome flag initialization
% order of flags => up,down,left,right cues
d_cue       = NaN(n_act,1)          ; % display cue     ? (0/1)
s_cue       = NaN(n_act,1)          ; % sample  cue     ? (0/1)
o_cue       = NaN(n_act,1)          ; % sample outcome  ? (0/1/NaN = not sampled)

% Evaluation of display, cue, and outcome flags
% -------------------------------------------------------------------------
% cycle over possible action outcomes
for a = 1:n_act

   
    
    % if action outcome takes the agent of grid, do not display cue and do not sample
    if any(act_oc(a,:) < 1) || any(act_oc(a,:) > n_dim)

        d_cue(a) = 0;
        s_cue(a) = 0;

    % if action outcome does not take agent of the grid 
    else

        % display cue
        d_cue(a)  = 1;

        % set sampling flag to zero as default
        s_cue(a)  = 0;

        % if, in addition action outcome decreases l1 distance to any
        % target, set sampling flag to 1
        for t = 1:size(l1_map,3)

             % action outcome decreases distance to target t
            if l1_map(act_oc(a,1),act_oc(a,2),t) < l1_map(pos(1),pos(2),t)
  
                % sample cue
                s_cue(a) = 1;
            end
        end
    end
end

% sample the relevant cues from a Bernoulli distribution with parameter 
% determined from the probability map
for a = 1:n_act
   
    % determine whether to sample or not
    if s_cue(a)
        
        % sample cue outcome based on Binomial distribution with n = 1 and
        % p as function of probability map and agent position
        o_cue(a) = binornd(1,p_map(pos(1), pos(2)));
       
    end
end

% determine cue array index
% -------------------------------------------------------------------------
% cue_array order: (1) dark_hor,(2) light_hor,(3) dark_ver,(4) light_ver           

% initialize cueidx to no cue depiction
cueidx = zeros(n_act,1);

% cycle over action set
for a = 1:n_act

    % update cueidx, if a cue is to be displayed
    if d_cue(a) 
       
        % up/down action outcomes -> horizontal bars
        if a < 3
        
            % a cue is to be displayed, but it was not sampled -> dark bar
            if s_cue(a) == 0
                cueidx(a) = 1;
            % a cue is to be displyed and it was sampled
            else
                % determine dark or light horizontal bar
                switch o_cue(a)
                
                    % 0 was sampled
                    case 0
                        % dark horizontal bar is displayed
                        cueidx(a) = 1;
                    % 1 was sampled
                    case 1
                        % light horizontal bar is displayed
                        cueidx(a) = 2;
                end % switch
            end % if
                   
        % left/right action outcomes -> vertical bars
        else
            
            % a cue is to be displayed, but it was not sampled -> dark bar
            if s_cue(a) == 0
                cueidx(a) = 3;
            % a cue is to be displyed and it was sampled
            else       
                % determine dark or light vertical bar
                switch o_cue(a)
                    % 0 was sampled
                    case 0
                        % dark vertical bar is displayed
                        cueidx(a) = 3;
                    % 1 was sampled
                    case 1
                        % light vertical bar is displayed
                        cueidx(a) = 4;
                end % switch
            end% if
        end % if
    end % if 
end % for
end % function

function p_12               = join_p(p_1,p_2)

% This function evaluates weights w_1, w_2 for two Bernoulli distribution
% parameters p_1, p_2 in [0.5,1] such that the weighted sum of p_1 and p_2
% using w_1, w_2, i.e.
%
%                      p_12 = w_1*p_2 + w_2*p_2
%
% falls into the interval [.5,1]. 
%
%       Inputs
%           p_1,p_2 : Bernoulli distribution parameters in [0.5,1]
%          
%       Outputs
%           w_1,w_2 : Bernoulli distribution parameter weights in [0,1]
%
% Copyright (C) Dirk Ostwald
% -------------------------------------------------------------------------
    if p_1 >= p_2
         p_12 = p_2 + (p_1 - p_2)*p_1;
    else
        p_12 = p_1 + (p_2 - p_1)*p_2;
    end
end

function [optglpathmidxs]   = optimal_foray(n,tgt_coord)

% This function is a training script that uses Dijkstra's algorithm to
% determine the optimal foray path to recover size(tgt,2) items in a n x n
% grid-world starting in the upper left corner. The implemenation of
% Dijkstra's algorithm capitalizes on graph notation.
%
%   Inputs
%           n               : scalar indicating the gridworld dimension
%           tgt_coords      : 2 x number of target array of target coordinates
%
%  Output
%           optglpathmidxs  : 2 x n array of the optimal path in matrix
%                             coordinates, n corresponds to the number of
%                             visited nodes
%
% Copyright (C) Dirk Ostwald
% -------------------------------------------------------------------------
clc
close all

% convert matrix target coordinates into graph nodes
% -------------------------------------------------------------------------
% CAVE: Matlab's linear indices increase row-wise, not column-wise,,
% while the optimal foray graph nodes do so vice versa. Hence the 
% transposition of the input target coordinates
tgt = NaN(1,size(tgt_coord,2));
for i = 1:size(tgt_coord,2)
    tgt(i) = sub2ind([n,n],tgt_coord(2,i),tgt_coord(1,i));
end

% tgt will be constructively destroyed below, thus save in additionally
tgt_mark = tgt;


% graph nodes: a grid world
% -------------------------------------------------------------------------
% specify matrix coordinates = graph nodes
M =      NaN(n^2,2);
m_idx   = 1;
for i = 1:n
    for j = 1:n
        M(m_idx,1) = i;
        M(m_idx,2) = j;
        m_idx      = m_idx + 1;
    end
end

% graph edges: a grid world adjacency matrix 
% -------------------------------------------------------------------------
% initialize as completely unconnected graph
A = zeros(size(M,1));

% specify square grid-world edges row wise
for i = 1:(n^2-1)
    
    % right-most column graph nodes
    if mod(i,n) == 0              
        
        A(i,i+n) = 1;
    
    % last-row column graph nodes
    elseif i >= (n^2 - n)
         
        A(i,i+1) = 1;
    
    % standard case
    else
        
        A(i,i+1) = 1;
        A(i,i+n) = 1;
        
    end
end

% specify symmetric connections
A = A + A';

% cost matrix for steps between nodes i and j
% -------------------------------------------------------------------------
C = A;

% Determine optimal path
% -------------------------------------------------------------------------
% initialize start node to upper left corner
start_node  = 1;

% initialize the optimal global path covering all targets
optglobalpath   = 1;

% iteratively determine the optimal global path
while ~isempty(tgt)

    % evaluate optimal paths from start_node to remaining target nodes
    pathcost = NaN(1,length(tgt));
    optpath  = cell(1,length(tgt));

    for i = 1:length(tgt)

        % determine costs and paths using Dijkstra's algorithm
        [localcost,localpath] = dijkstras_algorithm(A,C,start_node,tgt(i));
        pathcost(i)           = localcost;
        optpath{i}            = localpath;

    end

    % determine the cheapest path 
    [minpathcost, minidx]   = min(pathcost);

    % concatenate the global path, remove starting node
    optglobalpath              = [optglobalpath optpath{minidx}(2:end)];

    % redefine the starting point for the next search
    start_node              = tgt(minidx);
        
    % remove the target from the target list
    tgt(minidx)             = [];
    
end


% reconvert the optimal global path to matrix indices
% -------------------------------------------------------------------------
optglpathmidxs = NaN(2,length(optglobalpath));
for i = 1:length(optglobalpath)
    [col,row]               = ind2sub([n,n],optglobalpath(i));
    optglpathmidxs(1,i)    = row;
    optglpathmidxs(2,i)    = col;
end



% convert graph space into matrix coordinates into Cartesian coordinates
% ------------------------------------------------------------------------- 
% convert the matrix graph node coordinates to Cartesian coordinates
K   = NaN(n^2,2);
for m = 1:size(K,1)
   K(m,:)   = mat2cart(M(m,:),n);
end

% initialize connection array in matrix index space
src_node_mat = NaN(sum(sum(A)),2);
tgt_node_mat = NaN(sum(sum(A)),2);
src_node_crt = NaN(sum(sum(A)),2);
tgt_node_crt = NaN(sum(sum(A)),2);

% cycle over adjacency matrix rows and colummn
minidx = 1;
for i = 1:size(A,1)
    for j = 1:size(A,2)
        
        % there exists a connection between M coordinates n and m
        if A(i,j) == 1
            src_node_mat(minidx,:) = M(i,:);
            tgt_node_mat(minidx,:) = M(j,:);
            minidx                 = minidx + 1;
                
       end     
    end
end

% convert to list of Cartesian coordinates
for i = 1:size(src_node_mat,1)
    src_node_crt(i,:) = mat2cart(src_node_mat(i,:),n);
    tgt_node_crt(i,:) = mat2cart(tgt_node_mat(i,:),n);
end

% convert the path defined as edge nodes into matrix coordinates
% -------------------------------------------------------------------------
% initialize path coordinated arrays in matrix and Cartesian space
path_src_node_mat = NaN(length(optglobalpath)-1,2);
path_tgt_node_mat = NaN(length(optglobalpath)-1,2);
path_src_node_crt = NaN(length(optglobalpath)-1,2);
path_tgt_node_crt = NaN(length(optglobalpath)-1,2);

minidx = 1;
for i = 2:length(optglobalpath)
    
    % determine source and target nodes in matrix coordinates
    path_src_node_mat(minidx,:) = M(optglobalpath(i-1),:);
    path_tgt_node_mat(minidx,:) = M(optglobalpath(i)  ,:);
    minidx                      = minidx + 1;
                
end     

% convert matrix coordinate path into list of Cartesian coordinate path
for i = 1:size(path_src_node_mat,1)
    path_src_node_crt(i,:) = mat2cart(path_src_node_mat(i,:),n);
    path_tgt_node_crt(i,:) = mat2cart(path_tgt_node_mat(i,:),n);
end



% visualize the problem and its solution
% -------------------------------------------------------------------------
% visualize the grid world
do_plot = 0;

if do_plot
h = figure;
set(h, 'Color', [1 1 1])
hold on
plot(K(:,1), K(:,2), 'ko', 'MarkerFaceColor', 'k')
for i = 1:size(src_node_mat,1)
    plot([src_node_crt(i,1) tgt_node_crt(i,1)],[src_node_crt(i,2) tgt_node_crt(i,2)])
end

% visualize the target nodes
for i = 1:length(tgt_mark)
    plot(K(tgt_mark(i),1), K(tgt_mark(i),2), 'go', 'MarkerFaceColor', 'g')
end

% visualize the optimal path
for i = 1:size(path_src_node_crt,1)
    plot([path_src_node_crt(i,1) path_tgt_node_crt(i,1)],[path_src_node_crt(i,2) path_tgt_node_crt(i,2)], 'r', 'LineWidth', 3)
end

xlim([min(K(:,1))-1, max(K(:,1))+1])
ylim([min(K(:,2))-1, max(K(:,2))+1])
axis square
axis off
end


end

function [xy]               = mat2cart(ij,n)

% This function transforms matrix coordinates(i,j)to Cartesian coordinates
% (x,y) to allow for straightforward plotting of matrix space properties.
% It implements the mapping
%
%  f:N_n^2 -> R^2, (i,j) |-> f(i,j) := (x,y) := (j - (n+1)/2,-i +(n+1)/2 ) 
%
% where (i,j) refer to matrix row and column indices and (x,y) to Cartesian
% coordinates for a given square matrix size of size  n x n. The Cartesian 
% coordinates are zero centered at the center of the matrix
%
% Inputs
%          
%       ij  : (2 x 1) array of matrix row index i and column index j
%       n   :  matrix size 
%
% Outputs
%      xy   : (2 x 1) array of Cartesian x-ordinate x and y-ordinate y 
%
% Copyright (C) Dirk Ostwald
% -------------------------------------------------------------------------

% implement the transform
% -------------------------------------------------------------------------
xy(1) =   ij(2) - ((n + 1)/2);
xy(2) =  -ij(1) + ((n + 1)/2);

end
  
function [costs, paths]     = dijkstras_algorithm(A,C,s,f)

% This function implements Dijkstra's algorithm, a dynamic programming
% routine, to evaluate the optimal (lowest cost) path from a specified 
% start node to a specified target node on a graph specified by the 
% adjacency matrix A with an associated cost matrix for each edge given by
% C. It is based on Joseph Kirk's Matlab implementation of Dijkstra's 
% algorithm as available from http://www.mathworks.com/matlabcentral/
% fileexchange/20025-advanced-dijkstras-minimum-path-algorithm
%
%   Inputs
%
%       A : (n^2 x n^2) adjaency matrix specifying the edges of a graph
%            comprising n nodes
%       C : (n^2 x n^2) cost matrix specifying the cost associated with
%            taking a specific edge
%       s : scalar start node
%       t : scalar target node
%
%   Outputs
%
%       costs:
%       paths:
%
% Copyright (C) Dirk Ostwald
% -------------------------------------------------------------------------
% recover adjacency matrix size
[n, nc] = size(A);

% recover cost matrix size
[m, mc] = size(C);

% convert adjacency matrix to edge list
[I,J]   = find(A);
E       = [I J];
cost    = C;

L       = length(s);
M       = length(f);
costs   = zeros(L,M);
paths   = num2cell(NaN(L,M));

% Find the Minimum Costs and Paths using Dijkstra's Algorithm
for k = 1:L
    
    % Initializations
    TBL         = sparse(1,n); 
    min_cost    = Inf(1,n);
    settled     = zeros(1,n);
    path        = num2cell(nan(1,n));
    I           = s(k);
    min_cost(I) = 0;
    TBL(I)      = 0;
    settled(I)  = 1;
    path(I)     = {I};

    while any(~settled(f))
        
        % Update the Table
        TAB     = TBL;
        TBL(I)  = 0;
        nids    = find(E(:,1) == I);
        
        % Calculate the Costs to the Neighbor Points and Record Paths
        for kk = 1:length(nids)
        
            J = E(nids(kk),2);
            
            if ~settled(J)
                
                c = cost(I,J);
                empty = ~TAB(J); 
            
                if empty || (TAB(J) > (TAB(I) + c))
                    TBL(J) = TAB(I) + c;
                    path{J} = [path{I} J];
                else
                    TBL(J) = TAB(J);
                end
            end
        end
        
        K = find(TBL);
        
        % Find the Minimum Value in the Table
        N = find(TBL(K) == min(TBL(K)));
        
        if isempty(N)
            break
        else
            % Settle the Minimum Value
            I           = K(N(1));
            min_cost(I) = TBL(I);
            settled(I)  = 1;
        end
    end
    
    % store Costs and Paths
    costs(k,:) = min_cost(f);
    paths(k,:) = path(f);
    
end

if L == 1 && M == 1
    paths = paths{1};
end
end
