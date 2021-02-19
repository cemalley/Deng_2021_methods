# continuation of noci-deseq2.R plotting percent of genes expressed over normalized count of 100----

length(expressed.D28.100ct[Category=='GPCR', GeneId]) /length(genelist.all[Category=='GPCR', GeneId]) *100 #22.55639

length(expressed.D28.100ct[Category=='Neuropeptide', GeneId]) /length(genelist.all[Category=='Neuropeptide', GeneId]) *100 #25

length(expressed.D28.100ct[Category=='Nuclear receptor', GeneId]) /length(genelist.all[Category=='Nuclear receptor', GeneId]) *100 #60.41667

length(expressed.D28.100ct[Category=='Neurotransmitter transporter', GeneId]) /length(genelist.all[Category=='Neurotransmitter transporter', GeneId]) *100 #56

length(expressed.D28.100ct[Category=='Protein kinase', GeneId]) /length(genelist.all[Category=='Protein kinase', GeneId]) *100 #75.88933
length(expressed.D28.100ct[Category=='Ion channel', GeneId]) /length(genelist.all[Category=='Ion channel', GeneId]) *100 #75.46296 old list -- > 52.38095 new list

# graph-----
perc.data <- as.data.table(readxl::read_xlsx('/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/Percent_genes_expressed.xlsx', sheet='Combined'))
perc.data[,threshold:=as.character(threshold)]
perc.data[,category:=factor(category, levels=c('GPCR','Neuropeptide','Neurotransmitter transporter','Ion channel','Nuclear receptor','Protein kinase'))]

perc.data[,count_n := c(399, 84, 48, 25, 506, 411)]

perc.data[,category:= paste0(category, ' (',count_n,')')]

# ggplot(perc.data, aes(x=threshold, y=percent, fill=category))+ geom_bar(position='dodge', stat='identity')+
#   labs(x='Minimum normalized expression, mean', y='Percent expressed', title='Genes expressed in D28 Nociceptors by gene family category',
#        fill='Gene family')+
#   theme_bw()+
#   geom_text(aes(label=paste0(round(percent), '%')), position=position_dodge(width=0.9), vjust=-0.25, size=5)+
#   theme(axis.text=element_text(size=15), axis.title = element_text(size=15), title=element_text(size=16),
#         legend.text = element_text(size=12), legend.key.size = unit(2,"line"))+ scale_fill_viridis_d(option = 'D')

perc.data.100 <- perc.data[threshold==100 & category != 'Overall mean',]
perc.data.100[,category:=factor(category, levels=c('GPCR (399)','Neuropeptide (84)','Ion channel (411)','Neurotransmitter transporter (25)','Nuclear receptor (48)','Protein kinase (506)'))]



ggplot(perc.data.100, aes(x=category, y=percent, fill=category))+ geom_bar(position='dodge', stat='identity')+
  labs(x='Gene family category', y='Percent expressed', title='Genes expressed in D28 Nociceptors by gene family',
       fill='Gene family')+
  theme_bw()+
  geom_text(aes(label=paste0(round(percent), '%')), position=position_dodge(width=0.9), size=5, hjust=0)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15), title=element_text(size=15))+
  scale_fill_viridis_d(option = 'D', direction=-1)+guides(fill=FALSE) + coord_flip() + ylim(0,100)

ggplot(perc.data.100, aes(x=category, y=percent, fill=category))+ geom_bar(position='dodge', stat='identity')+
  labs(x='Gene family category', y='Percent expressed', title='Genes expressed in D28 Nociceptors by gene family',
       fill='Gene family')+
  theme_bw()+
  geom_text(aes(label=paste0(round(percent), '%')), position=position_dodge(width=0.9), size=5, hjust=0)+
  theme(axis.text.y=element_blank(), axis.text.x=element_text(size=15), title=element_text(size=15))+
  scale_fill_viridis_d(option = 'D', direction=-1)+guides(fill=FALSE) + coord_flip() + ylim(0,100)


# ions only graph

ion.channels.perc[,proportion := as.numeric(proportion)]

ggplot(ion.channels.perc[proportion >0,], aes(y=proportion, x=reorder(Subcategory, proportion), fill=reorder(Subcategory, proportion)))+ geom_bar(position='dodge',stat='identity') + coord_flip()+
  labs(y='Subcategory', x='Percent expressed', title='Ion channel genes expressed in D28 Nociceptors\n by gene family subcategory',
       fill='Gene family')+
  theme_bw()+
  geom_text(aes(label=paste0(round(proportion), '%')), position=position_dodge(width=0.9), size=5, hjust=0)+
  theme(axis.text=element_text(size=15), axis.title = element_text(size=15), title=element_text(size=16))+
  scale_fill_viridis_d(option = 'D', direction = -1)+guides(fill=FALSE) + ylim(0, 100)

ggplot(ion.channels.perc[proportion >0,], aes(y=proportion, x=reorder(Subcategory, proportion), fill=reorder(Subcategory, proportion)))+ geom_bar(position='dodge',stat='identity') + coord_flip()+
  labs(y='Subcategory', x='Percent expressed', title='Ion channel genes expressed in D28 Nociceptors\n by gene family subcategory',
       fill='Gene family')+
  theme_bw()+
  geom_text(aes(label=paste0(round(proportion), '%')), position=position_dodge(width=0.9), size=5, hjust=0)+
  theme(axis.text.y=element_blank(), axis.text.x=element_text(size=15), title=element_text(size=15))+
  scale_fill_viridis_d(option = 'D', direction = -1)+guides(fill=FALSE) + ylim(0, 100)

save(ion.channels.perc, file='/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/Ionchannels.RData')
save(perc.data.100, file='/Volumes/ncatssctl/NGS_related/BulkRNA/ISB012/Tao/Barplotdata.RData')






