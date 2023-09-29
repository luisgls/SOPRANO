# Output dN/dS file

The output of SOPRANO is a TSV file in the following format:

| coverage    | ON_dnds | ON_lowci | ON_highci | ON_muts | OFF_dnds | OFF_lowci | OFF_highci | OFF_muts | Pvalue | ON_na | ON_NA | ON_ns | ON_NS | OFF_na | OFF_NA | OFF_ns | OFF_NS |
|-------------| ------- | -------- | --------- | ------- | -------- | --------- | ---------- | -------- |--------| ----- | ----- | ----- | ----- | ------ | ------ | ------ | ------ |
| Exonic_Only |
| Exonic_Intronic |

Only first row (`Exonic_Only`) of the table data will be generated if there are 
no intronic mutations. The SOPRANO algorithm uses intronic mutations to improve
the background counts of silent mutations.

### Contents descriptions

- ```ON_dnds``` dN/dS of the target region provided in the bed file

- ```ON_lowci``` lower value for the 95% CI of the target

- ```ON_highci``` upper value for the 95% CI of the target

- ```ON_muts``` number of mutations observed inside the target region

- ```OFF_dnds``` dN/dS of the OFF-target region provided in the bed file

- ```OFF_lowci``` lower value for the 95% CI of the OFF-target

- ```OFF_highci``` upper value for the 95% CI of the OFF-target

- ```OFF_muts``` number of mutations observed outside the target region

- ```P-val``` P-value estimated from the comparison of the confidence intervals from ON and OFF dN/dS values

- ```ON_na``` Observed number of nonsilent mutations ON target

- ```ON_NA``` Number of nonsilent sites (corrected) ON target

- ```ON_ns``` Observed number of silent mutations ON target

- ```ON_NS``` Number of silent sites (corrected) ON target

- ```OFF_na``` Observed number of nonsilent mutations OFF target

- ```OFF_NA``` Number of nonsilent sites (corrected) OFF target

- ```OFF_ns``` Number of silent sites (corrected) OFF target

- ```OFF_NS``` Number of silent sites (corrected) OFF target