library(DataCombine)

#set path----
path=("C:/Users/mrademaker/Documents/Research projects/STCNWS/DATRAS/cpue_length_hour/")

#file with species info----
species_info=read.table(paste(path,"Data/spec_info.txt",sep=""),sep="\t",header=TRUE)
species_info$file_name =gsub(" ","_",species_info$Scientific.name)

spec=species_info$file_name[1]

data=read.csv(paste(path,sprintf("Seascapes/Q1/%s/GAM/Diff_smooth/num_diff.csv",spec),sep=""))
combis=c("1 - 2"  ,"1 - 3" ,"2 - 3" ,"2 - 4" ,"3 - 5" ,"3 - 6" ,"4 - 5" ,"5 - 10" ,"5 - 6"  ,"5 - 7" ,"5 - 8" ,"5 - 9" ,"6 - 8" ,"7 - 8" ,"7 - 9" ,"8 - 10" ,"8 - 9" ,"9 - 10")
row510=data[35,]
data <- data[data$combi %in% combis, ]
data=data[-c(12),]
data=InsertRow(data, NewRow = row510, RowNum = 8)
spec
min_ciw=min(data$ci_width)
max_ciw=max(data$ci_width)
range=max_ciw/min_ciw
print(range)
