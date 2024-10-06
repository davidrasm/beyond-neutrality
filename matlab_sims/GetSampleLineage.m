function this_lineage = GetSampleLineage(focal_sample, parent)

this_lineage = focal_sample; cntr = 1;

while 1 
    lineage_parent = parent(this_lineage(cntr));
    this_lineage(cntr+1) = lineage_parent;
    if (isnan(lineage_parent)) % got to the first infected individual  
        return;
    end
    cntr = cntr + 1;
end