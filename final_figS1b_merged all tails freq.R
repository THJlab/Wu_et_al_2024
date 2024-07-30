library(dplyr)
library(ggplot2)

setwd("./softclip_counts")

filenames <- list.files(pattern=".counts")
for (x in 1:length(filenames)){
  file_name=filenames[[x]]
  tmp=read.table(filenames[[x]], stringsAsFactors = F)
  names(tmp)=c('motif','occurences')
  assign(file_name,tmp)
}
rm(tmp,file_name,filenames,x)

#makes one list of all dfs
all_df<- Filter(function(x) is(x, "data.frame"), mget(ls()))

rm(list=setdiff(ls(), "all_df"))

setwd("..")

setwd("./spike_softclip_counts")

filenames <- list.files(pattern=".counts")
for (x in 1:length(filenames)){
  file_name=filenames[[x]]
  tmp=read.table(filenames[[x]], stringsAsFactors = F)
  names(tmp)=c('motif','spike_occurences')
  assign(file_name,tmp)
}
rm(tmp,file_name,filenames,x)

#makes one list of all dfs
all_spikes<- Filter(function(x) is(x, "data.frame"), mget(ls()))

rm(list=setdiff(ls(), c("all_df",'all_spikes')))

setwd("..")






setwd("./Plots")

#make one df of all others
Main_df=bind_rows(all_df, .id = "column_label")

simplified=dplyr::select(Main_df, -c('column_label'))

simplified_deduplicated <- simplified %>%                                    
  group_by(motif) %>%
  dplyr::summarise(occurences = sum(occurences))



spike_df=bind_rows(all_spikes, .id = "column_label")

spike_simplified=dplyr::select(spike_df, -c('column_label'))

spike_simplified_deduplicated <- spike_simplified %>%                                    
  group_by(motif) %>%
  dplyr::summarise(spike_occurences = sum(spike_occurences))



merged_genomic_spikes=merge(simplified_deduplicated,spike_simplified_deduplicated, by='motif', all.x = TRUE)

merged_genomic_spikes$spike_occurences[is.na(merged_genomic_spikes$spike_occurences)] <- 0

merged_genomic_spikes=merged_genomic_spikes %>%
  mutate(genomic_cnt=occurences-spike_occurences)

#start to work on the data (legth, freq,..)

merged_genomic_spikes$length=nchar(merged_genomic_spikes$motif)

tails_1_to_9=merged_genomic_spikes %>%
  filter(.$length <10)


tails_1_to_9_no_pure_A=tails_1_to_9[10:nrow(tails_1_to_9),]

rm(list=setdiff(ls(), 'tails_1_to_9_no_pure_A'))


#bellow code could use a proper rewriting either using functions and apply or better sorting methods
no_A_tails1=tails_1_to_9_no_pure_A %>%
  filter(.$length==1)

no_A_tails2=tails_1_to_9_no_pure_A %>%
  filter(.$length==2)

no_A_tails3=tails_1_to_9_no_pure_A %>%
  filter(.$length==3)

no_A_tails4=tails_1_to_9_no_pure_A %>%
  filter(.$length==4)

no_A_tails5=tails_1_to_9_no_pure_A %>%
  filter(.$length==5)

no_A_tails6=tails_1_to_9_no_pure_A %>%
  filter(.$length==6)

no_A_tails7=tails_1_to_9_no_pure_A %>%
  filter(.$length==7)

no_A_tails8=tails_1_to_9_no_pure_A %>%
  filter(.$length==8)

no_A_tails9=tails_1_to_9_no_pure_A %>%
  filter(.$length==9)



ranked_no_A_tails1=no_A_tails1[order(no_A_tails1$genomic_cnt,decreasing = TRUE),]
ranked_no_A_tails2=no_A_tails2[order(no_A_tails2$genomic_cnt,decreasing = TRUE),]
ranked_no_A_tails3=no_A_tails3[order(no_A_tails3$genomic_cnt,decreasing = TRUE),]
ranked_no_A_tails4=no_A_tails4[order(no_A_tails4$genomic_cnt,decreasing = TRUE),]
ranked_no_A_tails5=no_A_tails5[order(no_A_tails5$genomic_cnt,decreasing = TRUE),]
ranked_no_A_tails6=no_A_tails6[order(no_A_tails6$genomic_cnt,decreasing = TRUE),]
ranked_no_A_tails7=no_A_tails7[order(no_A_tails7$genomic_cnt,decreasing = TRUE),]
ranked_no_A_tails8=no_A_tails8[order(no_A_tails8$genomic_cnt,decreasing = TRUE),]
ranked_no_A_tails9=no_A_tails9[order(no_A_tails9$genomic_cnt,decreasing = TRUE),]

top9_ranked_no_A_tails2=ranked_no_A_tails2[1:9,]
top9_ranked_no_A_tails2$rank=c(1:9)
top9_ranked_no_A_tails3=ranked_no_A_tails3[1:9,]
top9_ranked_no_A_tails3$rank=c(1:9)
top9_ranked_no_A_tails4=ranked_no_A_tails4[1:9,]
top9_ranked_no_A_tails4$rank=c(1:9)
top9_ranked_no_A_tails5=ranked_no_A_tails5[1:9,]
top9_ranked_no_A_tails5$rank=c(1:9)
top9_ranked_no_A_tails6=ranked_no_A_tails6[1:9,]
top9_ranked_no_A_tails6$rank=c(1:9)
top9_ranked_no_A_tails7=ranked_no_A_tails7[1:9,]
top9_ranked_no_A_tails7$rank=c(1:9)
top9_ranked_no_A_tails8=ranked_no_A_tails8[1:9,]
top9_ranked_no_A_tails8$rank=c(1:9)
top9_ranked_no_A_tails9=ranked_no_A_tails9[1:9,]
top9_ranked_no_A_tails9$rank=c(1:9)

ranked_no_A_tails1$rank=c(1:3)

dfs_for_plot=list(ranked_no_A_tails1,
                  top9_ranked_no_A_tails2,
                  top9_ranked_no_A_tails3,
                  top9_ranked_no_A_tails4,
                  top9_ranked_no_A_tails5,
                  top9_ranked_no_A_tails6,
                  top9_ranked_no_A_tails7,
                  top9_ranked_no_A_tails8,
                  top9_ranked_no_A_tails9)




plotting_df=bind_rows(dfs_for_plot, .id = "column_label")

simplified_plotting_df=data.frame(motif=plotting_df$motif, true_cnt=plotting_df$genomic_cnt, length=plotting_df$length, rank=plotting_df$rank)

rm(list=setdiff(ls(), 'simplified_plotting_df'))

#replace Ts by Us
simplified_plotting_df$motif <- gsub('T', 'U', simplified_plotting_df$motif)
simplified_plotting_df$Utail=0

simplified_plotting_df$Utail[c(2,9,20,27,34,44,52,62,71)]=c('TRUE')

simplified_plotting_df=simplified_plotting_df %>%
  mutate(Utail=ifelse(Utail==0,'FALSE',Utail))

ggplot(simplified_plotting_df,
       aes(x=factor(length),
           y=factor(rank),
           fill=log10(true_cnt),
           label=motif)) +
  geom_tile() +
  geom_text(aes(color=Utail), size=5) +
  scale_color_manual(values=c('lightgray', 'red')) +
  xlab('softclip length') +
  ylab('ordered by abundance of softclip')
ggsave('motif_count.pdf')


















