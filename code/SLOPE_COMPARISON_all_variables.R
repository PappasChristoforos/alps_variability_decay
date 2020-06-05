# load libraries
library(dplyr)


#####
# 1 # --------------------------------------------------------------------------
##### Combine slopes of SAP, DEN, and TRW

# SAPFLOW data -----------------------------------------------------------------   
load("D:/postpro/SAP_Data_GS_Decaystdev.RData")

lm_SAP_h_diff = subset(SAP_h_diff_stdev, Hours >= 24) %>% 
    group_by(ID) %>% 
    do({
      mod = lm(log10(Stdev) ~ log10(Hours), data = .)
      data.frame(Slope = coef(mod)[2],
                 StdError = summary(mod)$coefficients[2,2],      
                 R2adj = summary(mod)$adj.r.squared)
    })
lm_SAP_h_diff$Type  = "W_SFD_H"

lm_SAP_h_diff$Elevation = NA
lm_SAP_h_diff[substr(lm_SAP_h_diff$ID, 2,3) == "08", "Elevation"] = 800
lm_SAP_h_diff[substr(lm_SAP_h_diff$ID, 2,3) == "13", "Elevation"] = 1300
lm_SAP_h_diff[substr(lm_SAP_h_diff$ID, 2,3) == "16", "Elevation"] = 1600
lm_SAP_h_diff[substr(lm_SAP_h_diff$ID, 2,3) == "19", "Elevation"] = 1900
lm_SAP_h_diff[substr(lm_SAP_h_diff$ID, 2,3) == "22", "Elevation"] = 2200

lm_SAP_h_diff$Species = "Picea abies"
lm_SAP_h_diff[substr(lm_SAP_h_diff$ID, 7,7) == "L", "Species"] = "Larix decidua"


# DENDROMETER data -------------------------------------------------------------
load("D:/postpro/Dendro_Data_GS_clipped_Decaystdev.RData")

R_h_diff_stdev$Variable = "deltaR_H"
GRO_zwf_diff_stdev$Variable = "G_H"
DR_h_stdev$Variable = "W_R_H"
                                 
lm_W_R_H = subset(DR_h_stdev, Hours >= 24) %>% 
    group_by(ID) %>% 
    do({
      mod = lm(log10(Stdev) ~ log10(Hours), data = .)
      data.frame(Slope = coef(mod)[2],
                 StdError = summary(mod)$coefficients[2,2],      
                 R2adj = summary(mod)$adj.r.squared)
    })
lm_W_R_H$Type  = "W_R_H"
 
lm_GRO_zwf_diff = subset(GRO_zwf_diff_stdev, Hours >= 24) %>% 
    group_by(ID) %>% 
    do({
      mod = lm(log10(Stdev) ~ log10(Hours), data = .)
      data.frame(Slope = coef(mod)[2],
                 StdError = summary(mod)$coefficients[2,2],      
                 R2adj = summary(mod)$adj.r.squared)
    })
lm_GRO_zwf_diff$Type  = "G_H"
 
lm_R_h_diff = subset(R_h_diff_stdev, Hours >= 24) %>% 
    group_by(ID) %>% 
    do({
      mod = lm(log10(Stdev) ~ log10(Hours), data = .)
      data.frame(Slope = coef(mod)[2],
                 StdError = summary(mod)$coefficients[2,2],      
                 R2adj = summary(mod)$adj.r.squared)
    })
lm_R_h_diff$Type  = "R_h_diff"

# make a data.frame
dendro_climacs_slopes = data.frame(rbind(lm_W_R_H, lm_GRO_zwf_diff, lm_R_h_diff))

dendro_climacs_slopes$Elevation = NA
dendro_climacs_slopes[substr(dendro_climacs_slopes$ID, 2,3) == "08", "Elevation"] = 800
dendro_climacs_slopes[substr(dendro_climacs_slopes$ID, 2,3) == "13", "Elevation"] = 1300
dendro_climacs_slopes[substr(dendro_climacs_slopes$ID, 2,3) == "16", "Elevation"] = 1600
dendro_climacs_slopes[substr(dendro_climacs_slopes$ID, 2,3) == "19", "Elevation"] = 1900
dendro_climacs_slopes[substr(dendro_climacs_slopes$ID, 2,3) == "22", "Elevation"] = 2200

dendro_climacs_slopes$Species = "Picea abies"
dendro_climacs_slopes[substr(dendro_climacs_slopes$ID, 7,7) == "L", "Species"] = "Larix decidua"


# TREE RING WIDTH data ---------------------------------------------------------
load("D:/postpro/slopes_df_full_TRWs_112yr.RData")

slopes_TRWs = data.frame(ID=slopes_df_full_TRWs_112yr$.id, 
                         Slope=slopes_df_full_TRWs_112yr$Slope,
                         StdError=slopes_df_full_TRWs_112yr$StdError,
                         R2adj=slopes_df_full_TRWs_112yr$R2adj,
                         Type="TRW",
                         Elevation=slopes_df_full_TRWs_112yr$Elevation,
                         Species=slopes_df_full_TRWs_112yr$Species)


# Put all variables together ---------------------------------------------------
slopes_df = rbind(dendro_climacs_slopes, 
                  data.frame(lm_SAP_h_diff),
                  data.frame(slopes_TRWs))

slopes_df_selc = slopes_df[slopes_df$Type %in% c("G_H", "TRW", "W_R_H", "W_SFD_H"), ]

slopes_df_selc$Type = factor(slopes_df_selc$Type)


DEN_slopes = dendro_climacs_slopes
SAP_slopes = lm_SAP_h_diff
TRW_slopes = slopes_TRWs
Ltal_slopes = slopes_df_selc


# save data.frames
save(Ltal_slopes, DEN_slopes, SAP_slopes, TRW_slopes, 
     file = "D:/postpro/Ltal_slopes.RData")


#####
# 2 # --------------------------------------------------------------------------
##### Boxplots

# load libraries
library(ggplot2)
#library(ggsignif)
library(grid)
library(gtable)
library(gridExtra)
library(scales)

# load Ltal data
load("D:/postpro/Ltal_slopes.RData")

# Growth Boxplot -----------------------------------------------------------------
growth_slopes_Ltal = Ltal_slopes[Ltal_slopes$Type %in% c("G_H", "TRW"), ]
       
# load EU TRW data 
load("D:/postpro/TRW_slopes_EU.Rdata")     
TRW_slopes_EU$Type = "TRW_EUR"
TRW_slopes_EU$Elevation = TRW_slopes_EU$elev 

growth_gg = data.frame(
             rbind(growth_slopes_Ltal[,c("Slope", "Type", "Elevation", "Species")],
                   TRW_slopes_EU[,c("Slope", "Type", "Elevation", "Species")]
                   )
             )
         
cols = c("chartreuse3","mediumseagreen", "NA") 

GG_growth = ggplot(data= growth_gg, aes(x=Species, y=Slope)) + 
  geom_boxplot(aes(fill=Type), outlier.shape = NA) +   

#  geom_signif(comparisons = list(c("Larix decidua", "Picea abies")),
#  test="t.test", test.args=list(alternative = "less"), map_signif_level=T)+
  scale_y_continuous("Slope [-]", limits=c(-0.5, 0))+
  xlab("")+
  scale_fill_manual("", values=cols)+
  theme_bw()+
             theme(aspect.ratio = 1.5,
                   panel.grid.minor.x = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   strip.background = element_blank())
  


# Hydro Boxplot -----------------------------------------------------------------
hydro_slopes_Ltal = Ltal_slopes[Ltal_slopes $Type %in% c("W_R_H", "W_SFD_H"), ]
       
cols = c("royalblue1", "skyblue2") 

GG_hydro = ggplot(data= hydro_slopes_Ltal, aes(x=Species, y=Slope, fill=Type)) + 
  geom_boxplot() +    

#  geom_signif(comparisons = list(c("Larix decidua", "Picea abies")),
#  test="t.test", test.args=list(alternative = "less"), map_signif_level=T)+
  scale_y_continuous("Slope [-]", limits=c(-1.1, -0.9))+
  xlab("")+
  scale_fill_manual("", values=cols)+
  theme_bw()+
             theme(aspect.ratio = 1.5,
                   panel.grid.minor.x = element_blank(),
                   panel.grid.minor.y = element_blank(),
                   panel.grid.major.x = element_blank(),
                   panel.grid.major.y = element_blank(),
                   strip.background = element_blank())


g2 <- ggplotGrob(GG_hydro)
g3 <- ggplotGrob(GG_growth)
g <- cbind(g2, g3, size = "first")
g$heights <- unit.pmax(g2$heights, g3$heights)
grid.newpage()
grid.draw(g)

ggsave("D:/Figs/Slopes_bpolt.eps", plot=g, device="eps", width=20, units="cm")



#####
# 3 # --------------------------------------------------------------------------
##### Gradient plot

cols = c("yellowgreen","mediumseagreen", "royalblue1", "skyblue2", "black") 
                                      
# gradient plot
grad_plot = 
 ggplot(Ltal_slopes)+
  
  geom_hline(yintercept = c(-1, -0.5, 0), linetype = "dashed", colour="gray80")+ 
  
  geom_jitter(data=TRW_slopes_EU , mapping=aes(group=siteID, x=Elevation, y=Slope),
              shape=21, fill="gray50", colour="white", show.legend = FALSE, 
              size=1, width = 50)+  
  
  geom_smooth(data=TRW_slopes_EU, mapping=aes(x=Elevation, y=Slope),
              colour="black", method="loess", se=T, show.legend = FALSE, alpha=1)+
             
  geom_smooth(mapping=aes(x=Elevation, y=Slope, colour=Type, 
                          group=interaction(Species, Type)),
              method="loess", se=T, show.legend = FALSE, alpha=1)+

  geom_jitter(mapping=aes(group=ID, x=Elevation, y=Slope, fill=Type),
              shape=21, colour="white", show.legend = FALSE, size=2, width = 50)+
                                                             
  geom_tile(mapping=aes(group=ID, x=Elevation, y=Slope, fill=Type))+
                                        
  scale_x_continuous("Elevation [m asl.]", 
                     limits=c(0,2500), 
                     breaks=c(0, 800, 1300, 1600, 1900, 2200))+
             
  scale_y_continuous("Slope [-]", 
                     limits=c(-1.2,0.1), 
                     breaks=c(-1, -0.5, 0))+
            
  scale_colour_manual(values=cols)+
  scale_fill_manual("", values=cols)+
            
  facet_grid(.~Species)+
  
  theme_bw()+
  
  theme(aspect.ratio = 1.4,
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_blank(),
        strip.background = element_blank())+
  
  coord_flip()
       
ggsave("D:/Figs/slopes_gradient.eps", plot=grad_plot, device="eps", width=16, units="cm")


