#RNA-seq workshop _ https://hackmd.io/SnorsWTbTTyRenptpjrhww?view

install.packages("ggplot2")
install.packages("tidyr")
install.packages("dplyr")
install.packages("cowplot")
library(ggplot2)
library(tidyr)
library(dplyr)
library(scales)
library(cowplot)




favorite_genes <- c("BRCA1", "JUN",  "GNRH1", "TH", "AR")
favorite_genes

# 1. set wd and import rnaseq sample file
setwd("xyz")
getwd()
samples <- read.csv("rnaseq/data/samples.csv") #contains sample info

head(samples) #count files can be very long and wide therefore we only view part of the data
tail(samples)
dim(samples) #shows no. of samples (rows) and sample features (columns)

str(samples) #compactly displays internal structure
summary(samples) #summary() computes statistics

# 2. check number of samples per tissue; then per tissue and sex; then sex, age, and hardy scale?
dplyr::count(samples, SMTS) #counts # of sample in sample dataframe; SMTS is the tissue type column
dplyr::count(samples, SMTS, SEX)
head(dplyr::count(samples, SMTS, SEX, AGE, DTHHRDY))
#Can we test the effect of sex on gene expr. in all tissues? No, some tissues are assigned to only one sex ie uterus
#Do you have enough samples to test the effects of Sex, Age, and Hardy Scale in the Heart? 


# 3. import count data file; then check the no. of samples (ie rows?) and no. columns (no. of features per sample) 
# file contains differential gene expression analysis of heart tissue from 20-29 vs 50-59 year olds - investigates effect of age (but not sex or hardy scale) on gene expression in the heart; contains numerical values such as p-value, LogFC etc
counts <- read.csv("rnaseq/data/countData.HEART.csv.gz", row.names = 1) #use gene id (ie ensemble id) as row names; contains  numerical info for each transcript
dim(counts)
View(counts) #confirm data

results <- read.table("rnaseq/data/GTEx_Heart_20-29_vs_50-59.tsv")
head(results)

# 4. create a visual plot of how many samples there are per tissue, sex, and hardiness
samples_subset <- dplyr::count(samples, SMTS, SEX, DTHHRDY)
ggplot(samples, aes(
              x=SMTS, fill = SEX)) +
  geom_bar(stat = "count") +
  facet_wrap(~DTHHRDY, nrow=1) +
  theme(axis.text.x = element_text(angle = 90), axis.text = element_text(size = 6))

##OR easier to read; according to tissue type and age
ggplot(samples, aes(
             x=SMTS, fill = DTHHRDY)) +
  geom_bar(stat = "count") +
  facet_wrap(~SEX) + #groups data by sex
  coord_flip()+ #makes x-axis easier to read
  theme(axis.text.x = element_text(angle = 90), axis.text = element_text(size = 6))

ggplot(samples, aes(x = AGE, fill = DTHHRDY))  +
  geom_bar(stat = "count") +
  facet_wrap(~SEX) 


# 5. create V-plot 
#Volcano Plots: type of scatter plots showing log fold change (logFC) x axis vs.  inverse log of a p-value (corrected for multiple hypothesis testing (adj.P.Val)) --> -log10(adj.P.Val)
ggplot(results, aes(x = logFC, y=-log10(adj.P.Val))) +
  geom_point()

# inverse log of p<05 is '-log10(0.05)'; add horizontal line to observe no. of significant genes (or plot points)
ggplot(results, aes(x = logFC, y=-log10(adj.P.Val), colour = adj.P.Val)) + #colour points along a gradient
  geom_point() +
  geom_hline(yintercept= -log10(0.05), linetype=2, colour = "black") 


# create plot to allow colour to show significant genes:
# `aes(color = ifelse( adj.P.Val < 0.05, "p < 0.05", "NS"))`. This essentially  called   and gives each factor a color different color on the plot.
ggplot(results, aes(x = logFC, y=-log10(adj.P.Val), colour = adj.P.Val)) +
  geom_point(aes(colour = ifelse(adj.P.Val < 0.05, "p < 0.05", "NS"))) + #create two factors, "p < 0.05" and "NS", and gives each a different colour
  geom_hline(yintercept= -log10(0.05), linetype=2) #select line type

# clean plot; move legend to bottom, add titles, captions etc.
heart_20_50_plot <- ggplot(results, aes(x = logFC, y=-log10(adj.P.Val), colour = adj.P.Val)) +
              geom_point(aes(colour = ifelse(adj.P.Val < 0.05, "p < 0.05", "NS"))) + 
              geom_hline(yintercept= -log10(0.05), linetype=2) +
              geom_vline(xintercept= c(-0.5, 0.5), linetype=2) +
               # theme(legend.position = "bottom") +
               labs(color = "adjusted P.Val", 
                    title = "Heart Tissue Gene Expression", 
                    subtitle = "expr. data for 20-29 vs 50-59 year olds",
                    caption = "Volcano plot displaying log fold change (logFC) vs. inverse log of a p-value")



# 6. Are there any interactions between RIN score (SMRIN) and sequencing facility (SMCENTER)?
# create box-plots using samples file containing RIN scores (data quality)
ggplot(samples, aes(x=SMCENTER, y=SMRIN)) +
  geom_boxplot()

# then add individual points
ggplot(samples, aes(x=SMCENTER, y=SMRIN)) +
  geom_boxplot() +
  geom_jitter(aes(color = SMRIN)) #shortcut for geom_point(position="jitter") - adds small amount of random variation of each point, useful handling overplotting caused by discreteness in smaller datasets




# Wrangle - Tidy and Transform Data
# makes data more valuable and assures quality and usefulness


# 7. identify approved names and symbols of differentially expressed genes (DEGs) by filtering by adjusted p-values<0.05 then further filtered by fold-change if needed (ie too many results to handle; too many with small fold-changes)
# subset both results and genes files
#
results_subset <- filter(results, adj.P.Val < 0.05)
results_subset <- filter(results, adj.P.Val < 0.05, logFC > 1 | logFC < -1) %>% #filters results further
                arrange(adj.P.Val) %>%     #sorts order by p-value ascending
                    head()
results_subset

results_DEGs <- rownames(results_subset) # shows names of differentially expr. genes filtered using specifics
  View(results_DEGs)


# 8. CHECKPOINT: Mutating data
# most RNA-pipelines require count file formatting: sample = column and genes = rows; values = integers or doubles
# however, we then need a corresponding counts file formatted: row names = sample id matching column names

# import colData - a dataframe has metadata for each sample, containing sample identifier, relevant experimental factors (ie. treatment/control, cell type, tissue, age etc)
colData <- read.csv("rnaseq/data/colData.HEART.csv", row.names = 1)
head(colData)


# rownames(colData) == colnames(counts) should show TRUE statements, if FALSE data cannot be processed
head(rownames(colData) == colnames(counts)) #check whether row names are the same as column names
head(colnames(counts))
head(rownames(colData)) #we can see dashes in colData have been placed with dots in counts data

# 10. Identify DEGs
#DESeq2 prefers dots - so use gsub() to replace all; sub() replaces first occurrence
tidy_colData <- colData %>% 
        mutate(SAMPID = gsub("-", ".", SAMPID)) #adds/ mutate new variable, SAMPID to now be replaced by with gsub function

rownames(tidy_colData) <- tidy_colData$SAMPID #this uses SAMPID-containing .dots. instead of rownames containing -dashes-  
View(tidy_colData)

# importing and modifying original files take too long due to size therefore we create "tidy" versions for downstream analyses
# first create variable of row names and name 'mycols' of tidied data
# then use the row names of tidied colData (ie mycols) to represent them as columns in countData
mycols <- rownames(tidy_colData)
counts_tidy <- counts %>%
  select(all_of(mycols)) #all_of function - also checks for missing characters in a selective way

head(rownames(tidy_colData) == colnames(counts_tidy)) #should now show TRUE



# 11. Joining Data
# we need to identify and match gene names from the counts "genes.txt" file, in this case the EnsembleIDs in the 'id' column, to the results file which, instead, uses gene symbols

# this particular counts file uses Ensemble IDs, import /read "genes.txt" 
# combine this with results file to append gene info, in this case symbols column
genes <- read.table("rnaseq/data/ensembl_genes.tsv", sep = "\t", header = TRUE, fill = TRUE) #wout "fill=" error message shows a row doesn't have the same #of columns as the rest
View(genes)


resultsSymbol <- results %>%
      mutate(name = row.names(results))
head(resultsSymbol)

# combine two data frames to contain results, gene names and info
#`left_join` incl all records from the left table and any matching values from the right. 
#`right_join` incl all values from the right table and any matching values from the left. 
#`inner_join` incl records that have values in both tables. 
#`full_join` incl everything. 
resultsName <- left_join(resultsSymbol, genes, by = "name") #adds all of resultsSymbol and any matches from genes according to 'name' col
View(resultsName)

# filter and select columns to make a pretty table of the DEGs using final combined df
View(results_DEGs) # checks differentially expr. genes 
results_DEGs2 <- resultsName %>%
      filter(adj.P.Val < 0.05, logFC > 1 | logFC < -1) %>%
      arrange(adj.P.Val) %>%  #sort according to p-value
      select(name, description, id, logFC, AveExpr, adj.P.Val)  #selects column order
View(results_DEGs2)

# list Ensemble IDs of DEGs
DEGs_ids <- results_DEGs2 %>%
      drop_na(id) %>% # drops rows where column contains missing values
      pull(id) #similar to '$'; converts column in dataframe to list
DEGs_ids


# 12. Lengthening data
# The matrix form of countData is required for some pipelines, but a long format (each row is an observation ie sample) is better in R programs are  suited to 
# create a tidy long vers. of countData by first filtering (due to data size) using 'rowSums(.) > 0' ; removes data with values of 0
# this allows the data to be easily subset by variables or genes of interest (GOIs ie rownames) 
counts_tidy2 <- counts_tidy %>%
  mutate(id = row.names(.)) %>%   # use ensembleID column as rownames 
  filter(id %in% DEGs_ids)  #then filter id column with ensembleIDs of DEGs
View(counts_tidy2)  #for some reason order of function matters here


# we can pivot the dataframe to make rows longer and lesser columns
# create df by using all col names (ie sample id) and adding to new SAMPLID col, then add count values
count_tidy2_long <- counts_tidy2 %>%
  pivot_longer(cols = all_of(mycols),  # cols - specifies column names (variables) turned into observations (rows ie samples);
               names_to = "SAMPID",   # names_to - specifies name of new column
               values_to = "counts")   # values_to - names the column with corresponding values 
View(count_tidy2_long)


#now we have matching column names; join the genes & tidy_colData dataframes to count_tidy2_long according to specified col
counts_colData_join <- count_tidy2_long %>%
    left_join(.,tidy_colData, by = "SAMPID") %>%
  left_join(., genes, by = "id")
View(counts_colData_join)


# now plot we can plot the counts for DEGs
# 'scales' package allows scientific notation for the axes

final_counts_DEG_plot <- ggplot(counts_colData_join, aes(AGE, counts, fill = SEX)) +
   geom_boxplot() +
   geom_jitter() +
  facet_wrap(~name, scales = "free_y") +   #free_y scale; diff for each graph
  theme(axis.text.x = element_text(angle = 75, hjust  = 1),
          strip.text = element_text(face = "italic")) +
  scale_y_log10(labels = label_number_si()) #transforms the y-axis using log10 and labels the numbers with SI prefixes 
final_counts_DEG_plot










#### Challenges

# 1. How to read the following files:
#   a. GTEx results comparing *muscle* of 20-29 year old to 50-59 year olds?
#   b. csv info file describing muscle samples? 
#   c. how about GTEx results comparing *muscle* of 20-29 year old to 70-79 year olds?

# 2. How to import 'data/colData.MUSCLE.csv' and count the no. of muscles samples per sex, age? 
# 3. How many female muscles samples are there from age group 30-39?

# 4. Create a volcano plot to compare *muscle* tissue of 20-29 year olds vs 50-59 year olds? (HINT: log fold change (logFC) x axis vs. -log10(adj.P.Val))
# 5. Are there more or less DEGs in the muscle compared to the heart for this age group?
# 6. Do the same comparison for the *muscle* tissue of 20-29 vs. 70-79 year olds?

#1a.
muscle_20_50_results <- read.table("rnaseq/data/GTEx_Muscle_20-29_vs_50-59.tsv")
View(muscle_20_50_results)
#1b.
View(samples) #already saved
#1c.
muscle_20_70_results <- read.table("rnaseq/data/GTEx_Muscle_20-29_vs_70-79.tsv")
View(muscle_20_70_results)
#2.
colData_muscle <- read.csv("rnaseq/data/colData.MUSCLE.csv", row.names = 1) #differential gene expr.
View(colData_muscle)
#3.
count(colData_muscle, SEX, AGE) #shows 11 types of samples that fit 'criteria'; there are three 30-39yo female muscle samples
#4.
muscle_20_50_plot <- ggplot(muscle_20_50_results, aes(logFC, -log10(adj.P.Val), colour = ifelse(adj.P.Val < 0.05, "p < 0.05", "NS"))) +
  geom_point() +
  geom_hline(yintercept= -log10(0.05), linetype=2)
muscle_20_50_plot
#4.1 make plot nicer
#first create new 'Expression' column showing gene regulation as categories
muscle_20_50_results2 <- muscle_20_50_results %>%
  mutate(Expression = case_when(logFC >= 0.05 & adj.P.Val <= 0.05 ~ "UPREGULATED",   #if logFC is more than/equal to 0.05 and significant then it means its upregulated expr.
                                logFC <= -0.05 & adj.P.Val <= 0.05 ~ "DOWNREGULATED",
                                TRUE ~ "nosig")) #case_when vectorises if_else statements
View(muscle_20_50_results2)

#then plot
muscle_20_50_plot2 <- ggplot(muscle_20_50_results2, aes(logFC, -log10(adj.P.Val), colour = Expression)) +
  geom_point() +
  scale_color_manual(values = c("darkslategray2", "gray72", "coral1")) + #colour codes https://www.nceas.ucsb.edu/sites/default/files/2020-04/colorPaletteCheatsheet.pdf
  geom_jitter(alpha = 0.25, width = 0.2) + # adds random noise ie transparency to points and limit its width
  geom_hline(yintercept= -log10(0.05), linetype=2) +
  geom_vline(xintercept= c(-0.5, 0.5), linetype=2)
muscle_20_50_plot2

#then add other info
muscle_20_50_plot3 <- ggplot(muscle_20_50_results2, aes(logFC, -log10(adj.P.Val), colour = Expression)) +
  geom_point() +
  scale_color_manual(values = c("darkslategray2", "gray72", "coral1")) +
  geom_jitter(alpha = 0.25, width = 0.2) +
  geom_hline(yintercept= -log10(0.05), linetype=2) +
  geom_vline(xintercept= c(-0.5, 0.5), linetype=2) +
  #xlim(-3, 3) + #manually changes x-axis limit
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-0.5, 13)) +
  xlab(expression("log"[2]*"FC")) +
  ylab(expression("-log"[10]*"adj.P.Val")) + #alters y-axis-label; subscript 10
  labs(color = "adjusted P.Val", 
       title = "Muscle Tissue Gene Expression", 
       subtitle = "expr. data for 20-29 vs 50-59 year olds",
       caption = "Volcano plot displaying log fold change (logFC) vs. inverse log of adj. p-value")
muscle_20_50_plot3

# 5.
install.packages("cowplot")
library(cowplot)
#re-create heart plot to have the same axis scale as muscle
heart_20_50_plot2 <- heart_20_50_plot +
  scale_x_continuous(limits = c(-3, 3)) +
  scale_y_continuous(limits = c(-0.5, 13))
heart_20_50_plot2 

plot_grid(muscle_20_50_plot3, heart_20_50_plot2) #more DEGs in the muscle > heart






  


### Descriptions of functions
#----
  | Function         | Description |
  |------------------|-------------|
  | `log10()` | A built-in function for a log transformation  |
  | `read.csv()`  | A base R function for importing comma separated tabular data  |
  | `read_csv()`  | A tidyR function for importing .csv files as tibbles |
  | `read.table()` | A base R function for importing tabular data with any delimiter |
  | `read_tsv()`  | A tidyR function for importing .tsv files as tibbles | 
  | `as_tibble()` | Convert data frames to tibbles | 
  | `head()` and `tail()` | Print the first or last 6 lines of an object  | 
  | `dim()`  | A function that prints the dimensions of an object | 
  | `length()` | Calculate the length of an object |
  | `count()` | A dplyr function that counts number of samples per group |
  | `str()` | A function that prints the internal structure of an object  |  
  | `summary()` | A function that summarizes each variable |  
  
  | `ggplot2` | An open-source data visualization package for the statistical programming language R | 
  | `ggplot()`  | The function used to construct the initial plot object, and is almost always followed by + to add component to the plot |
  | `aes()`  | Aesthetic mappings that describe how variables in the data are mapped to visual properties (aesthetics) of geoms |
  | `geom_point()`  | A function used to create scatter plots | 
  | `geom_bar()`  | A function used to create bar plots | 
  | `coord_flip()`  | Flips the x and y axis | 
  | `geom_hline()`  | Add a horizontal line to plots | 
  
  | `filter()`       |  A function for filtering data |
  | `mutate()`       |  A function for create new columns  |
  | `select()`       |  A function for selecting/reordering columns   |
  | `arrange()`      |  A function for ordering observations  |
  | `full_join()`    |  Join 2 tables, return all observations |
  | `left_join()`    |  Join 2 tables, return all observations in the left and matching observations in the right table |
  | `inner_join()`   |  Join 2 tables, return observations with values in both tables        |
  | `pivot_wider()`  |  Widen a data frame  |
  | `pivot_longer()` |  Lengthen a data frame |
  | `drop_na()`      |  Remove missing values  |
  | `separate()`     |  Separate a column into two columns   |
  
  ----Operators & accessors
  `[]` indexes into a vector, matrix, array, list or dataframe; also used to extract specific elements from a vector or matrix

  `[[]]` extracts element in a list?
  
  `^` represents the beginning of the string
  `$` represents the end of the string
  `.` stands for any character 
  `*` defines a repetition (zero or more)

  
