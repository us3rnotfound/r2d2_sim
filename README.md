# Deliverables:

### r2d2_sim.m
Matlab/Octave script showing platoon maneuvering around Remote Radio Detection Devices (R2D2s) from (0,0) to the maximal (x,y) point.

### config.txt
Required configuration file to be in same directory as r2d2_sim.m.  It defines platoon movement for each of the possible sensor "heatmap" scenarios.

### jyn_prompt2_battlespace.pptx
MS Powerpoint presentation detailing the plan to utilize a computer simulation (r2d2_sim.m) to simulate the maneuvering of Platoon through battlespace to avoid detection from the R2D2s.

# r2d2_sim.m Basic Instructions:

   Run in Octave or Matlab.  I changed the way the platoon.directive cell array was setup to make it compatible with Matlab.
   It opens up a split GUI figure with battlespace on top and platoon detection time values along the bottom.
   Platoon automatically starts from (1,1), and goes until it finishes at (1000,1000).
   In the bottom of the GUI window are the count and the R2D2 Detections count.
   The Count is simply the main loop counter.  The R2D2 Detections count is the
   number of times the R2D2 (enemy) units spot the platoon (hint: this is rare unless you decrease the MAX_RANGE variable).
  
   **It is required to have config.txt present in the same directory.**
   

### input arguments:
   - 'GRID_SIZE', <#> (default 1000)
   - 'NUM_R2D2', <#>  (Try 5 - 50)
   - 'RD2D_RANGE', <#> (default 5)
   - 'MAX_RANGE', <#>  (Try 50, soldier sensor distance of r2d2 detection)
   
### Example Calls (Ran from Octave/Matlab Main Console):
   - r2d2_sim
   - r2d2_sim('GRID_SIZE', 500, 'NUM_R2D2', 20)
   - r2d2_sim('MAX_RANGE', 50, 'GRID_SIZE', 1000)
   - r2d2_sim('GRID_SIZE', 1000, 'NUM_R2D2', 40, 'R2D2_RANGE', 6, 'MAX_RANGE', 20)
