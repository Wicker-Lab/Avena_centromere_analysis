## This script contains the code for calculating the centromere boundaries

# read in data
ChIP_bamcompare_500kb_MQ30 <- read.table("CRR515328_CRR515329_500kb_MQ30.bedgraph")

# initiate dataframe (only once)
df_centromere <- data.frame(chromosome = character(), start = numeric(), end = numeric())

#change (e.g., "2D")
chr <- ""

#identify the window in which 25% (threshold) is met first and last
ChIP_bamcompare_500kb_MQ30_chr <- ChIP_bamcompare_500kb_MQ30[ChIP_bamcompare_500kb_MQ30$V1 == paste0("Asat_OT3098_v2_chr", chr),]
ChIP_bamcompare_500kb_MQ30_chr[which(ChIP_bamcompare_500kb_MQ30_chr$V4 > max(ChIP_bamcompare_500kb_MQ30_chr$V4)*0.25),]

#individual threshold per chromosome
max(ChIP_bamcompare_500kb_MQ30_chr$V4)*0.25

### to calculate the boundaries of the centromere
## enter V2 (start) from first window over the threshold and the one before that and vice versa

# window start last before threshold is met (or b - window size)
a = 
  
# window start where threshold is first met 
b = 
  
delta_x = abs(a-b)
delta_y = abs((ChIP_bamcompare_500kb_MQ30_chr[ChIP_bamcompare_500kb_MQ30_chr$V2 == a,]$V4)-(ChIP_bamcompare_500kb_MQ30_chr[ChIP_bamcompare_500kb_MQ30_chr$V2 == b,]$V4))
delta_y2 = abs((max(ChIP_bamcompare_500kb_MQ30_chr$V4)*0.25)-(ChIP_bamcompare_500kb_MQ30_chr[ChIP_bamcompare_500kb_MQ30_chr$V2 == a,]$V4))
delta_x2 = (delta_x/delta_y)*delta_y2
start_cen = (a+delta_x2)
start_cen

# window start where threshold is last met
c = 
  
# window start from following window (or c + window size)
d = 
  
delta_x = abs(c-d)
delta_y = abs((ChIP_bamcompare_500kb_MQ30_chr[ChIP_bamcompare_500kb_MQ30_chr$V2 == c,]$V4)-(ChIP_bamcompare_500kb_MQ30_chr[ChIP_bamcompare_500kb_MQ30_chr$V2 == d,]$V4))
delta_y2 = abs((max(ChIP_bamcompare_500kb_MQ30_chr$V4)*0.25)-(ChIP_bamcompare_500kb_MQ30_chr[ChIP_bamcompare_500kb_MQ30_chr$V2 == d,]$V4))
delta_x2 = (delta_x/delta_y)*delta_y2
end_cen = (d-delta_x2)
end_cen

# collect data
new_data <- data.frame(chromosome = chr,  start = round(start_cen,0),  end = round(end_cen,0))

# feed data into dataframe
df_centromere <- rbind(df_centromere, new_data)