function [PE] = PEPF( est_group, group0 )




    PE = 1 - mean(true_class == est_group);
end
