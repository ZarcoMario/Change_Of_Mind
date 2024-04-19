# Change_Of_Mind
 Change of Mind 

The preprocessing step used to calculate changes of mind is 
to align trajectories to have a common initial point.
This is done with the method  align_trial_start in align_trajectories.py

You can find two examples for two and four targets in 
test_chom_two_targets.py and test_chom_four_targets.py, respectively.

The method get_changes_of_mind_two_targets in change_of_mind.py uses another method 
called samples_outside_region to define a small circle 
and check that the number of points (N) defined to consider a trial a change of mind
is higher outside that region. This is not a necessary step but we used it 
because we noticed what we considered "false positives" as 
there might be N points in the alternative region as a result of the a slow movement initiation.