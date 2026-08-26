fix_issue_3

Add OFF-mode workflow for cohort-level immunopeptidome analyses

When measuring selection across multiple patients, patient-specific OFF-target regions can differ substantially. This workflow provides a common OFF-target regions as the complement of the immunopeptidome intersections across all the patients in the cohort. 
In the OFF-mode, SOPRANO ON values correspond biologicall to OFF-target selection estimates

How to run the ON mode 
./run_localSSBselection_vLOCAL.sh \
-i data/MISSONI/MISSONI_ON_clonal_cohort_SOPRANO.anno \
-b immunopeptidomes/MISSONI/IP_MISSONI_clonal_merged.bed \
-m ssb192 \
-e false \
-n MISSONI_ON

How to run the OFF mode
./run_localSSBselection_vLOCAL_MOD4OFF.sh \
-i data/MISSONI/MISSONI_OFF_clonal_cohort_SOPRANO.anno \
-b immunopeptidomes/MISSONI/IP_MISSONI_clonal_intersection_complement_modified.bed \
-m ssb192 \
-e false \
-n MISSONI_OFF \
-o results
