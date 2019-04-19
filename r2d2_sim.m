%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andrew Payne
% 14 April 2019
%
% r2d2_sim.m - A simulator to show platoon moving through battlespace with
% algorithms programmed to dodge R2D2 units.
%
%
% Basic Instructions:
%   Run in Octave.  Opens up a split GUI figure with battlespace on top and
%   platoon detection time values along the bottom.  Platoon automatically starts 
%   from (1,1), and goes until it finishes at (1000,1000).
%   In the bottom of the GUI window are the count and the R2D2 Detections count.
%   The Count is simply the main loop counter.  The R2D2 Detections count is the
%   number of times the R2D2 (enemy) units spot the platoon.
% It is required to have config.txt present in the same directory.
%
% input arguments:
%   - GRID_SIZE, <#> (default 1000)
%   - NUM_R2D2, <#>  (Try 5 - 50)
%   - RD2D_RANGE, <#> (default 5)
%   - MAX_RANGE, <#>  (Try 50, soldier sensor distance of r2d2 detection)
%
% Example Calls:
%   r2d2_sim
%   r2d2_sim('GRID_SIZE', 500, 'NUM_R2D2', 20)
%   r2d2_sim('MAX_RANGE', 50, 'GRID_SIZE', 1000)
%   r2d2_sim('GRID_SIZE', 1000, 'NUM_R2D2', 40, 'R2D2_RANGE', 6, 'MAX_RANGE', 20)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ret_val = r2d2_sim(varargin)
  GRID_SIZE = 1000;
  NUM_R2D2 = 40;
  R2D2_RANGE = 10;
  MAX_RANGE = 40;

  % Parsing input arguments.
  for i = 1:2:length(varargin)
    switch upper(varargin{i})
      case 'GRID_SIZE'
        GRID_SIZE = varargin{i+1};
      case 'NUM_R2D2'
        NUM_R2D2 = varargin{i+1};
      case 'R2D2_RANGE'
        R2D2_RANGE = varargin{i+1};
      case 'MAX_RANGE'
        MAX_RANGE = varargin{i+1};
      otherwise
        msgbox(['Command ', varargin{i}, ' is not recognized. Exiting.']);
        ret_val = -1;
        return
    end
  end
  
  % Maneuvering step counts:
  DIRECTIVE_FORWARD_STEPS = GRID_SIZE/50;   % Seems to work the best.
  DIRECTIVE_RETRACE_STEPS = GRID_SIZE/50;
  
  config_file_contents = [];
  
  Read_Config_File;
  
  ret_val = 0;
  
  grid = zeros(GRID_SIZE);

  % Generate R2D2.
  r2d2.num = NUM_R2D2;
  % Positions (columns wrap top-to-bottom).
  r2d2.raw_positions = fix(rand(1,NUM_R2D2)*GRID_SIZE^2);
  % Set ordered pair positions via R2d2_Positions:
  r2d2.num_detections = 0;
  R2d2_Positions;

  % Place R2D2s into grid.
  grid(r2d2.raw_positions) = 1;

  % Place platoon at the starting mark.
  platoon.position = [1,1];
  platoon.directive = {'N/A', 0;
                       'N/A', 0;};
  platoon.last_move = [0 0 0 0]; % n s e w.
  platoon.move_count = 0;
  platoon.detections = zeros(DIRECTIVE_FORWARD_STEPS, 1);
  platoon.maneuver = 0;
  
  sensor.num_hits = 0;
  sensor.heatmap = zeros(2);  % Give a fresh 2x2 sensor relative heatmap grid.
  sensor.old_heatmap = zeros(2);
  sensor.change = zeros(2);
  splot=[];
  
  Setup_Plot;
  
  done = 0;

  % Main time-step loop.
  while (done == 0)
    % Turn on sensor. Collect Data about R2D2.
    % Gets a 2x2 grid of NW, SW, NE, SE data correlating to severity of sensor data.
    Platoon_Sense;

    % Decides how to proceed to the finish line based on sensor heatmap.
    Platoon_Decide;

    % R2d2 detection:
    R2d2_Detection;
    
    % Updates the plot for this simulation.
    Update_Plot;

    % Loop ending when finish line reached.
    done = (platoon.position(1) > GRID_SIZE-1 && platoon.position(2) > GRID_SIZE-1);
    
    %disp([platoon.directive{1}, ' ', num2str(platoon.directive{2})]);
    
    %disp([num2str(platoon.move_count), ' ', num2str(platoon.position)]);
  end
  
    disp([num2str(platoon.position(1)), ' ', num2str(platoon.position(2))]);
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% Read_Config_File - A function to read the config.txt file into a 
%                    cell array called config_file_contents.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function read_success = Read_Config_File
    read_success = 0;
    i = 1;
    fileID = fopen('config.txt');
    while ~feof(fileID)
      f_line = strtrim(fgetl(fileID));
      if length(f_line)
        % needs to have a : to denote key:value pair, or drop it.
        if strfind(f_line,':')  
          % Remove comments.
          k = min(strfind(f_line, '//'));  % min: find the first occurance.
          if k == 1
            % don't include, the whole line contains a comment.
          elseif k > 1
            % Partially-commented. Key-value pair.
            config_file_contents{i} = strsplit(f_line(1:k-1),':');
            i = i + 1;
          else  % include whole line. Key-value pair.
            config_file_contents{i} = strsplit(f_line,':');
            i = i + 1;
          end
        end
      end
    end
    fclose(fileID);
    
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% Get_Value -  A function to read the configuration from config.txt file, get to
%              the key specified, and return the cell array value(s) to the
%              right of the '<key>:' label in config.txt. Returns a cell array
%              of value(s).
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  function value = Get_Value(key)
    % Search for key
    if length(config_file_contents)
      for i=1:length(config_file_contents)
        if strfind(config_file_contents{i}{1}, key)
          % found the key, now take the values and break loop.
          value = strtrim(strsplit(config_file_contents{i}{2},','));  % split up by ',' delimeter.
          % Deletes white-space.
          break;
        end
      end
    else  
      msgbox('Unsucessful read of config.txt');
      return
    end
  end
  
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% Setup_Plot - A function to create a figure window and display the battlespace 
%              scatter plot and heatmap spikes.
%
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  function Setup_Plot
    for (i = 1:NUM_R2D2)
      splot.xdata(i) = r2d2.positions{i}(1);
      splot.ydata(i) = r2d2.positions{i}(2);
    end
    % Add the platoon position:
    splot.xdata(NUM_R2D2+1) = platoon.position(1);
    splot.ydata(NUM_R2D2+1) = platoon.position(2);
    
    f=figure('units', 'normalized', 'position', [.2 .2 .6 .6]);
    subplot(2,1,1);
    
    % Start the graph
    splot.graph = scatter(splot.xdata, splot.ydata);
    
    axis([0 GRID_SIZE*1.1 0 GRID_SIZE*1.1]);

    subplot(2,1,2);
    
    splot.xheat = linspace(1,1000,GRID_SIZE*5);
    splot.yheat = zeros(1,GRID_SIZE*5);
    splot.heatmap_plot = plot(splot.xheat, splot.yheat);
    splot.move_counter = text(1,4,'Count = 0');
    splot.r2d2_detections = text(1,3,'R2D2 Detections = 0');
    axis([0 GRID_SIZE*5 0 5]);
    
  end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% Update_Plot - A function to edit the platoon position, the counter text, the
%               heatmap indicator, and the R2D2 detection count.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  function Update_Plot
    xdata_tmp = get(splot.graph, 'xdata');
    ydata_tmp = get(splot.graph, 'ydata');
    
    xdata_tmp(NUM_R2D2+1) = platoon.position(1);
    ydata_tmp(NUM_R2D2+1) = platoon.position(2);
    
    set(splot.graph, 'xdata', xdata_tmp);
    set(splot.graph, 'ydata', ydata_tmp);
    
    yheat_tmp = get(splot.heatmap_plot,'ydata');
    yheat_tmp(platoon.move_count) = sum(sensor.heatmap(:));
    
    set(splot.heatmap_plot, 'ydata', yheat_tmp);
    set(splot.move_counter, 'string', ['Count = ',num2str(platoon.move_count)]);
    set(splot.r2d2_detections, 'string', ['R2D2 Detections = ',num2str(r2d2.num_detections)]);
    
    drawnow
     
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% R2d2_Positions - A function to install the R2D2s on the battlespace as deemed
%                  by r2d2.raw_positions array.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  function R2d2_Positions
    for i=1:NUM_R2D2
      % X position
      x_pos = ceil(r2d2.raw_positions(i)/GRID_SIZE);
      %x_pos = 50*i + 20;
      
      % Y position
      y_pos = mod(r2d2.raw_positions(i), GRID_SIZE);
      %y_pos = 50*i - 5;
      if (y_pos==0)
        y_pos = GRID_SIZE;
      end 
      
      % Set ordered pair variable.
      r2d2.positions{i} = [x_pos, y_pos]; 

    end
  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% R2d2_Detection - A function to check if R2D2 detected a platoon in accordance
%                  with the R2D2_RANGE variable.  It will increment the count
%                  within the R2D2 structure if it does spot a platoon.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  function R2d2_Detection
     % Go through list of r2d2s, calculating Manhatten Distance.
    for i = 1:NUM_R2D2
      raw_offsets{i} = r2d2.positions{i} - platoon.position;
      raw_man_dist{i} = sum(abs(raw_offsets{i}));
    end
    
    % Reorder from closest to furthest, based on Manhatten Distance:
    [man_dist, sort_order] = sort(cell2mat(raw_man_dist));
    man_dist = num2cell(man_dist);
    offsets = raw_offsets(sort_order); 
    
    R2d2_Measure;
    
    function R2d2_Measure
      r2d2.range = R2D2_RANGE;
      
      % Go through list of r2d2s.
      for r2d2_ctr = 1:NUM_R2D2
        % Manhatten Distance comparison with sensor.range.
        if (man_dist{r2d2_ctr} <= r2d2.range)
          % Add to list of hit(s)
          r2d2.num_detections = r2d2.num_detections + 1;
        else
          break; % Keeps getting further away, break.
        end
      end  
    end
  end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% Platoon_Decide - A function to process the beginning through the end of the
%                  platoon's movements. It is very high level. It starts with
%                  deciding whether it's the first step, if it is, just head
%                  north. If there's a number >0 of directives in the queue, then
%                  it calls Platoon_Manuever. Else it checks the sensor's heatmap
%                  and calls a new directive if it is registering an R2D2.
%                  Otherwise, a routine movement northeast is all that it calls.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
  function Platoon_Decide
    if (platoon.move_count == 0)
      Platoon_Move('North');
      return
    elseif (platoon.directive{2}>0)
      
      % Call the directives.
      Platoon_Manuever

  % Heatmap analysis.
  
  % -------------
  % | NW  | NE  |
  % |     |     |
  % -------------
  % | SW  | SE  |
  % |     |     |
  % -------------
  
  %
    elseif (sensor.num_hits)
      last_move = platoon.directive{1};
      dir_finder = zeros(2);
    
      % Find the direction of the r2d2 being sensed.
      [max_heat,~] = max(sensor.heatmap(:));
      
      [x, y] = find(sensor.heatmap==max_heat);
      
      % Find the direction of the r2d2s that are most increasing.
      %[~,max_chg_dir] = max(sensor.change(:));
      
      dir_finder(x,y) = 1;
      disp(dir_finder)
      %dir_finder(max_chg_dir) = 1;
      
      % Form a string that Get_Value function can read.
      % e.g., turns [0 1;0 0] into '[0 1;0 0'
      dir_finder_str = ['[',num2str(dir_finder(1)),' ',num2str(dir_finder(3)),...
                        ';',num2str(dir_finder(2)),' ',num2str(dir_finder(4)),...
                        ']'];
      
      % Get the directives from the config.txt file.
      directives = Get_Value(dir_finder_str);
      % Assign the directives as deemed from config.txt.
      if (length(directives) == 2)
        platoon.directive(1:2) = {directives{1}, str2double(directives{2})};
      elseif (length(directives) == 4)
        platoon.directive(1:2) = {directives{1}, str2double(directives{2})};
        platoon.directive(3:4) = {directives{3}, str2double(directives{4})};
      end
      
      % Starts the move. The next loop around, this will not run because there
      % will be queued up directives.    
      Platoon_Move(platoon.directive{1});

      % Just a routine move north or east to reach the finish line.
      else
        if platoon.last_move(1)
          Platoon_Move('East');
         % disp('routine move east.');
        % North
        else
          Platoon_Move('North');
         % disp('routine move north.');
        end
      end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% Platoon_Manuever - A function to handle calling Platoon_Move function from the
%                    directives queue, which is a multidimensional cell array
%                    holding the direction to travel and the number of steps to
%                    take on that direction.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%      
      % Maneuvering function.
      function Platoon_Manuever
      
        Platoon_Move(platoon.directive{1});
        % Decrement the count.
        platoon.directive{2} = platoon.directive{2} - 1;
        
        % When the first row of directives is exhausted, move second row up.
        if (platoon.directive{2} == 0 && platoon.directive{4}>0)
          platoon.directive{1} = platoon.directive{3};
          platoon.directive{2} = platoon.directive{4};
          platoon.directive{4} = 0;
        end
      
    end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% Platoon_Move - A function to move the Platoon in a direction (required input
%                argument string).  It records the last move before exiting the
%                function. It also keeps the platoon in-bounds.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
      function Platoon_Move(direction)
        direction = upper(direction);
        switch direction
          case 'NORTH'
            if (platoon.position(2) < GRID_SIZE)
              platoon.position(2) = platoon.position(2) + 1;
              platoon.last_move = [1 0 0 0];
              platoon.move_count = platoon.move_count + 1; 
            elseif (platoon.position(1) < GRID_SIZE)
              Platoon_Move('East');
            else
              Platoon_Move('South');
            end
          case 'SOUTH'
            platoon.position(2) = platoon.position(2) - 1;
            platoon.last_move = [0 1 0 0];
            platoon.move_count = platoon.move_count + 1; 
          case 'EAST'
            if (platoon.position(1) < GRID_SIZE)
              platoon.position(1) = platoon.position(1) + 1;
              platoon.last_move = [0 0 1 0];
              platoon.move_count = platoon.move_count + 1; 
            elseif (platoon.position(2) < GRID_SIZE)
              Platoon_Move('North');
            else
              Platoon_Move('South');
            end
          case 'WEST'
            platoon.position(1) = platoon.position(1) - 1;
            platoon.last_move = [0 0 0 1];
            platoon.move_count = platoon.move_count + 1; 
          case 'NORTHEAST'
            if (platoon.last_move == [1 0 0 0])
              Platoon_Move('East');
            else
              Platoon_Move('North');
            end
          case 'NORTHWEST'
            if (platoon.last_move == [1 0 0 0])
              Platoon_Move('West');
            else
              Platoon_Move('North');
            end
          case 'SOUTHEAST'
            if (platoon.last_move == [0 1 0 0])
              Platoon_Move('East');
            else
              Platoon_Move('South');
            end
          case 'SOUTHWEST'
            if (platoon.last_move == [0 1 0 0])
              Platoon_Move('West');
            else
              Platoon_Move('South');
            end
        end  
      end

  end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% Platoon_Sense - A function to use the sensor that the Platoon has to
%                 sense R2D2s.  It calls measure and calculation routines that
%                 convert the sensor data into usable information.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
  function Platoon_Sense
    % Go through list of r2d2s, calculating Manhatten Distance.
    % This part of the code is from the know-all computer, which knows every
    % r2d2.position.  Later it will be filtered by what the sensor actually sees
    % in Sensor_Measure.
    for i = 1:NUM_R2D2
      raw_offsets{i} = r2d2.positions{i} - platoon.position;
      raw_man_dist{i} = sum(abs(raw_offsets{i}));
    end
    
    % Reorder from closest to furthest, based on Manhatten Distance:
    [man_dist, sort_order] = sort(cell2mat(raw_man_dist));
    man_dist = num2cell(man_dist);
    offsets = raw_offsets(sort_order);
    
    % Run Sensor_Measure to pick up any R2D2(s).
    Sensor_Measure;
    
    % Run Sensor_Calc to calculate severity of measurements and mark the change from last detection.
    Sensor_Calc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% Sensor_Measure - A function to actually compare the omniscent computer
%                  generated offsets with the Platoon's inferior sensor
%                  (sensor.range is derived from a probability).  It will rack
%                  up the num_hits and catalog the manhatten distances.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    function Sensor_Measure
      sensor.capability = .2 + .5*rand; % luck of the draw.
      sensor.range = sensor.capability * MAX_RANGE;
      sensor.old_heatmap = sensor.heatmap;
      sensor.heatmap = zeros(2);  % Give a fresh 2x2 sensor relative heatmap grid.
      sensor.num_hits = 0; % reset.
      
      % Go through list of r2d2s.
      for r2d2_ctr = 1:NUM_R2D2
        % Manhatten Distance comparison with sensor.range.
        if (man_dist{r2d2_ctr} <= sensor.range)
          % Add to list of hit(s)
          sensor.num_hits = sensor.num_hits + 1;
          sensor.hits{r2d2_ctr} = offsets{r2d2_ctr};
          sensor.man_dist{r2d2_ctr} = man_dist{r2d2_ctr};
        else
          break; % Keeps getting further away, break.
        end
      end  
    end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% FUNCTION:
% Sensor_Calc - A function to take the data from the Sensor_Measure routine and
%               create the heatmap, the information that the Platoon uses in
%               Platoon_Decide to decide where to maneuver.
%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
    function Sensor_Calc
      for hit_ctr = 1:sensor.num_hits
        % Find the N/S/E/W direction via the sign.
        if (sensor.hits{hit_ctr}(1) > 0)
          x_sign = 1;
        else
          x_sign = -1;
        end
        if (sensor.hits{hit_ctr}(2) > 0)
          y_sign = 1;
        else
          y_sign = -1;
        end
        
        norm_vector = sensor.hits{hit_ctr}/norm(abs(sensor.hits{hit_ctr}));
        
        % If X direction is severe:
        if (norm_vector(2) ~= 0)
          if (abs(norm_vector(1)/norm_vector(2)) > 5)
            norm_vector(1) = 1;
            norm_vector(2) = 0;
          end
        end
        if (norm_vector(1) ~= 0)
          if (abs(norm_vector(2)/norm_vector(1)) > 5)
            norm_vector(1) = 0;
            norm_vector(2) = 1;
          end
        end
        
        % Gives the strength of the signals:
        x_gauge = norm_vector(1)*x_sign*(1 - abs(sensor.hits{hit_ctr}(1) / sensor.range));
        y_gauge = norm_vector(2)*y_sign*(1 - abs(sensor.hits{hit_ctr}(2) / sensor.range));
   
        if (x_gauge > 0) % Target is at right.
          sensor.heatmap(1,2) = sensor.heatmap(1,2) + x_gauge;
          sensor.heatmap(2,2) = sensor.heatmap(2,2) + x_gauge;
        elseif (x_gauge < 0) % Target is at left.
          sensor.heatmap(1,1) = sensor.heatmap(1,1) - x_gauge;
          sensor.heatmap(2,1) = sensor.heatmap(2,1) - x_gauge;
        end
        if (y_gauge > 0) % Target is above.
          sensor.heatmap(1,1) = sensor.heatmap(1,1) + y_gauge;
          sensor.heatmap(1,2) = sensor.heatmap(1,2) + y_gauge;
        elseif (y_gauge < 0) % Target is below.
          sensor.heatmap(2,1) = sensor.heatmap(2,1) - y_gauge;
          sensor.heatmap(2,2) = sensor.heatmap(2,2) - y_gauge;
        end
        
      end
      
      % Get the change from last detection:
      sensor.change = sensor.heatmap - sensor.old_heatmap;
      
    end
end

end
