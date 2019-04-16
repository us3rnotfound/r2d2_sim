%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Andrew Payne
% 14 April 2019
%
% r2d2_sim.m - A simulator to show platoon moving through battlespace with
% algorithms programmed to dodge R2D2 units.
%
% Basic Instructions:
%   Run in Octave.  Opens up a split GUI figure with battlespace on top and
%   platoon detection time values along the bottom.  Platoon automatically starts 
%   from (1,1), and goes until it finishes at (1000,1000).
%   In the bottom of the GUI window are the count and the R2D2 Detections count.
%   The Count is simply the main loop counter.  The R2D2 Detections count is the
%   number of times the R2D2 (enemy) units spot the platoon.
%
%   Modify the following parameters to test different scenarios:
%    - GRID_SIZE (default 1000)
%    - NUM_R2D2  
%    - R2D2_RANGE  (default 5)
%    - MAX_RANGE  (soldier sensor distance of r2d2 detection)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function ret_val = r2d2_sim
  GRID_SIZE = 1000;
  NUM_R2D2 = 40;
  R2D2_RANGE = 10;
  MAX_RANGE = 40;
  % Maneuvering step counts:
  DIRECTIVE_FORWARD_STEPS = 5;
  DIRECTIVE_RETRACE_STEPS = 20; 
  
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
  endwhile
  
    disp([num2str(platoon.position(1)), ' ', num2str(platoon.position(2))]);
 
 
  function Setup_Plot
    for (i = 1:NUM_R2D2)
      splot.xdata(i) = r2d2.positions{i}(1);
      splot.ydata(i) = r2d2.positions{i}(2);
    endfor
    % Add the platoon position:
    splot.xdata(NUM_R2D2+1) = platoon.position(1);
    splot.ydata(NUM_R2D2+1) = platoon.position(2);
    
    f=figure;
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
    
  endfunction
  
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
     
  endfunction

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
      endif 
      
      % Set ordered pair variable.
      r2d2.positions{i} = [x_pos, y_pos]; 

    endfor
  endfunction
  
  function R2d2_Detection
     % Go through list of r2d2s, calculating Manhatten Distance.
    for i = 1:NUM_R2D2
      raw_offsets{i} = r2d2.positions{i} - platoon.position;
      raw_man_dist{i} = sum(abs(raw_offsets{i}));
    endfor
    
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
        endif
      endfor  
    endfunction
  endfunction
  
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
      
      switch (dir_finder)
        case [0 0;
              0 1]
          platoon.directive(1:2)={'NorthWest', DIRECTIVE_FORWARD_STEPS};

        case [0 0;
              1 0]
          platoon.directive(1:2)={'NorthEast', DIRECTIVE_FORWARD_STEPS};
        
        case [0 0;
              1 1]
          platoon.directive(1:2)={'North', DIRECTIVE_FORWARD_STEPS};
       
        case [0 1;
              0 0]
          platoon.directive(1:2)={'South', DIRECTIVE_RETRACE_STEPS};
          platoon.directive(3:4)={'East', DIRECTIVE_RETRACE_STEPS};

        case [0 1;
              0 1]
          platoon.directive(1:2)={'NorthWest', DIRECTIVE_RETRACE_STEPS};
          platoon.directive(1:2)={'East', DIRECTIVE_RETRACE_STEPS};
        case [0 1;
              1 0]
          platoon.directive(1:2)={'SouthEast', DIRECTIVE_FORWARD_STEPS};
        
        case [0 1;
              1 1]
          platoon.directive(1:2)={'NorthWest', DIRECTIVE_FORWARD_STEPS};
        
        case [1 0;
              0 0]
          platoon.directive(1:2)={'SouthEast', DIRECTIVE_FORWARD_STEPS};
        
        case [1 0;
              0 1]
          platoon.directive(1:2)={'NorthEast', DIRECTIVE_FORWARD_STEPS};
        
        case [1 0;
              1 0]
          platoon.directive(1:2)={'East', DIRECTIVE_FORWARD_STEPS};
        
        case [1 0;
              1 1]
          platoon.directive(1:2)={'NorthEast', DIRECTIVE_FORWARD_STEPS};
        
        case [1 1;
              0 0]
          platoon.directive(1:2)={'South', DIRECTIVE_FORWARD_STEPS};
        
        case [1 1;
              0 1]
          platoon.directive(1:2)={'SouthWest', DIRECTIVE_RETRACE_STEPS};
        otherwise
          disp('should not see this');
      endswitch
      
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
        endif
      endif
      
      % Maneuvering function.
      function Platoon_Manuever
      
        Platoon_Move(platoon.directive{1});
        % Decrement the count.
        platoon.directive{2} = platoon.directive{2} - 1;
        
        % When the first row of directives is exhausted, move second row up.
        if (platoon.directive{2} == 0 && platoon.directive{4}>0)
          disp('Changing Directives');
          platoon.directive{1} = platoon.directive{3};
          platoon.directive{2} = platoon.directive{4};
          platoon.directive{4} = 0;
        endif
      
      endfunction
  
      % Basic moving function.
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
            endif
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
            endif
          case 'WEST'
            platoon.position(1) = platoon.position(1) - 1;
            platoon.last_move = [0 0 0 1];
            platoon.move_count = platoon.move_count + 1; 
          case 'NORTHEAST'
            if (platoon.last_move == [1 0 0 0])
              Platoon_Move('East');
            else
              Platoon_Move('North');
            endif
          case 'NORTHWEST'
            if (platoon.last_move == [1 0 0 0])
              Platoon_Move('West');
            else
              Platoon_Move('North');
            endif
          case 'SOUTHEAST'
            if (platoon.last_move == [0 1 0 0])
              Platoon_Move('East');
            else
              Platoon_Move('South');
            endif
          case 'SOUTHWEST'
            if (platoon.last_move == [0 1 0 0])
              Platoon_Move('West');
            else
              Platoon_Move('South');
            endif
        endswitch  
      endfunction

  endfunction
  
  function Platoon_Sense
    % Go through list of r2d2s, calculating Manhatten Distance.
    for i = 1:NUM_R2D2
      raw_offsets{i} = r2d2.positions{i} - platoon.position;
      raw_man_dist{i} = sum(abs(raw_offsets{i}));
    endfor
    
    % Reorder from closest to furthest, based on Manhatten Distance:
    [man_dist, sort_order] = sort(cell2mat(raw_man_dist));
    man_dist = num2cell(man_dist);
    offsets = raw_offsets(sort_order);
    
    % Run Sensor_Measure to pick up any R2D2(s).
    Sensor_Measure;
    
    % Run Sensor_Calc to calculate severity of measurements and mark the change from last detection.
    Sensor_Calc;

    function Sensor_Measure
      %sensor.capability = 0.8; % luck of the draw.
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
        endif
      endfor  
    endfunction
    
    function Sensor_Calc
      for hit_ctr = 1:sensor.num_hits
        % Find the N/S/E/W direction via the sign.
        if (sensor.hits{hit_ctr}(1) > 0)
          x_sign = 1;
        else
          x_sign = -1;
        endif
        if (sensor.hits{hit_ctr}(2) > 0)
          y_sign = 1;
        else
          y_sign = -1;
        endif
        
        norm_vector = sensor.hits{hit_ctr}/norm(abs(sensor.hits{hit_ctr}));
        
        % If X direction is severe:
        if (abs(norm_vector(1)/norm_vector(2)) > 5)
          norm_vector(1) = 1;
          norm_vector(2) = 0;
        elseif (abs(norm_vector(2)/norm_vector(1)) > 5)
          norm_vector(1) = 0;
          norm_vector(2) = 1;
        endif
        
        % Gives the strength of the signals:
        x_gauge = norm_vector(1)*x_sign*(1 - abs(sensor.hits{hit_ctr}(1) / sensor.range));
        y_gauge = norm_vector(2)*y_sign*(1 - abs(sensor.hits{hit_ctr}(2) / sensor.range));
   
        if (x_gauge > 0) % Target is at right.
          sensor.heatmap(1,2) = sensor.heatmap(1,2) + x_gauge;
          sensor.heatmap(2,2) = sensor.heatmap(2,2) + x_gauge;
        elseif (x_gauge < 0) % Target is at left.
          sensor.heatmap(1,1) = sensor.heatmap(1,1) - x_gauge;
          sensor.heatmap(2,1) = sensor.heatmap(2,1) - x_gauge;
        endif
        if (y_gauge > 0) % Target is above.
          sensor.heatmap(1,1) = sensor.heatmap(1,1) + y_gauge;
          sensor.heatmap(1,2) = sensor.heatmap(1,2) + y_gauge;
        elseif (y_gauge < 0) % Target is below.
          sensor.heatmap(2,1) = sensor.heatmap(2,1) - y_gauge;
          sensor.heatmap(2,2) = sensor.heatmap(2,2) - y_gauge;
        endif
        
      endfor
      
      % Get the change from last detection:
      sensor.change = sensor.heatmap - sensor.old_heatmap;
      
    endfunction
endfunction

endfunction
