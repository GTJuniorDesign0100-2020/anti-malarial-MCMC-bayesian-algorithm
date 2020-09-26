# memory.limit(size=50000)
# options(java.parameters = "-Xmx4096m")

import pandas as pd

##### read in data

def onload():

  inputfile = "Angola2017_example.xlsx";

  excel_file = pd.ExcelFile(inputfile);
  genotypedata_latefailures = excel_file.parse("Late Treatment Failures", skiprows=3);


#genotypedata_latefailures[genotypedata_latefailures == 0] = 'NA'; # missing data has to be coded as NA
#genotypedata_latefailures[genotypedata_latefailures == "0"] = 'NA'; # missing data has to be coded as NA
#genotypedata_latefailures[genotypedata_latefailures == "N/A"] = 'NA'; # missing data has to be coded as NA
#genotypedata_latefailures[genotypedata_latefailures == "-"] = 'NA'; # missing data has to be coded as NA
#genotypedata_latefailures[genotypedata_latefailures == "NA"] = 'NA'; # missing data has to be coded as NA



### recode sample names so that each pair has a " Day 0" and a " Day Failure"
  SampleID = genotypedata_latefailures["Sample ID"];
  for index, text in SampleID.iteritems():
    end = text[-3:];
    if end[0] == '_':
      SampleID[index] = text.replace(end, "_ Day 0");
    if end[0] == 'D':
      SampleID[index] = text.replace(end, " Day Failure");

    genotypedata_latefailures["Sample ID"] = genotypedata_latefailures["Sample ID"].replace([index], SampleID[index]);


# each sample in genotypedata_RR has to have day 0 and day of Failure
  patients = [];
  for index, text in SampleID.iteritems():
    patient = text[0:7];
    patients.append(patient);

  for patient in patients:
    count = patients.count(patient);
    if count % 2 != 0:
      print("Error - each sample must have day 0 and day of failure data");

### background samples (from "Additional" tab)


  additional_genotypedata = excel_file.parse("Additional", skiprows=3);
  return genotypedata_latefailures, additional_genotypedata


#  if (dim(additional_genotypedata)[1] > 0) {
#    additional_genotypedata[additional_genotypedata == 0] = NA # missing data has to be coded as NA
#    additional_genotypedata[additional_genotypedata == "0"] = NA
#    additional_genotypedata[additional_genotypedata == "N/A"] = NA
#    additional_genotypedata[additional_genotypedata == "-"] = NA
#    additional_genotypedata[additional_genotypedata == "NA"] = NA
#    additional_genotypedata$Sample.ID = sub("_D0"," Day 0",additional_genotypedata$Sample.ID)
#    additional_genotypedata$Sample.ID = sub("_D[0-9]*"," Day Failure",additional_genotypedata$Sample.ID)
#  }


# recode as numeric


#genotypedata_latefailures[columns[-2:]] = genotypedata_latefailures[columns[-2:]].map(function (x) int(character(genotypedata_latefailures[,x])));

#additional_genotypedata[,colnames(additional_genotypedata)[-c(1,2)]] = sapply(colnames(additional_genotypedata)[-c(1,2)],function (x) as.numeric(as.character(additional_genotypedata[,x])))
