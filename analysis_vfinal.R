library(ggplot2)
library(cowplot)
library(viridis)
library(dplyr)
library(tidyr)
theme_set(theme_cowplot())
# change this to the your directory containing the scripts
setwd("YOUR_WORKING_DIRECTORY")

#### Load data ####

# simulation varying epimutation rate
df_u = read.csv("RI_sim2.6.1_2024-04-02_20-56-22.csv")

# recalculate epimutaiton rate from i and c parameters
df_u$u = (df_u$i_1+df_u$c_1)/2
df_u$u = factor(df_u$u, levels = sort(unique(df_u$u),decreasing = F)) # make factor for plotting

# recalculate skew from i and c parameters
df_u$skew = (df_u$i_1-df_u$c_1)/(df_u$i_1+df_u$c_1)
# division by small numbers causes slight differences in skew values
# that should be the same, rounding removes these small differences
df_u$skew = round(df_u$skew, digits = 1) 
# fill in skew values when mu = 0, even though skew doesn't affect
# the results of these simulations this makes plotting easier later
df_u$skew[is.nan(df_u$skew)] = rep(c(-1,0,0.5,1),each = 144)
df_u$skew = factor(df_u$skew, levels = sort(unique(df_u$skew),decreasing = F)) # make factor for plotting

# create vector marking when mu = 0 for plotting the genetic locus
df_u$u0 = F
df_u$u0[df_u$u == 0] = T

# simulation varying epimutation rate
# with pre-selection epimutation
df_u_pre = read.csv("RI_sim2.6.1_2024-04-02_21-09-35.csv")
# same as above
df_u_pre$u = (df_u_pre$i_1+df_u_pre$c_1)/2
df_u_pre$u = factor(df_u_pre$u, levels = sort(unique(df_u_pre$u),decreasing = F))
df_u_pre$skew = (df_u_pre$i_1-df_u_pre$c_1)/(df_u_pre$i_1+df_u_pre$c_1)
df_u_pre$skew = round(df_u_pre$skew, digits = 1)
df_u_pre$skew[is.nan(df_u_pre$skew)] = rep(c(-1,0,0.5,1),each = 144)
df_u_pre$skew = factor(df_u_pre$skew, levels = sort(unique(df_u_pre$skew),decreasing = F))
df_u_pre$u0 = F
df_u_pre$u0[df_u_pre$u == 0] = T

# simulation varying skew
df_skew = read.csv("RI_sim2.6.1_2024-04-02_21-22-14.csv")
df_skew$u = (df_skew$i_1+df_skew$c_1)/2
df_skew$u = factor(df_skew$u, levels = sort(unique(df_skew$u),decreasing = F))
df_skew$skew = (df_skew$i_1-df_skew$c_1)/(df_skew$i_1+df_skew$c_1)
df_skew$skew = round(df_skew$skew, digits = 1)
df_skew$skew = factor(df_skew$skew, levels = sort(unique(df_skew$skew),decreasing = F))
df_skew$u0 = F
df_skew$u0[df_skew$u == 0] = T

# simulation measuring RI over time
df_t = read.csv("RI_sim2.6.1_2024-04-02_21-22-28.csv")
df_t$u = (df_t$i_1+df_t$c_1)/2
df_t$u = factor(df_t$u, levels = sort(unique(df_t$u),decreasing = F))
df_t$skew = (df_t$i_1-df_t$c_1)/(df_t$i_1+df_t$c_1)
df_t$skew = round(df_t$skew, digits = 1)
df_t$skew[is.nan(df_t$skew)] = 1
df_t$skew = factor(df_t$skew, levels = sort(unique(df_t$skew),decreasing = F))
df_t$u0 = F
df_t$u0[df_t$u == 0] = T

# create vectors to select frequency columns that are homozygous E
monoE1 = colnames(df_u)[grepl("EE",colnames(df_u)) &
                        grepl("freq",colnames(df_u)) &
                        grepl("1",colnames(df_u))]
monoE2 = colnames(df_u)[grepl("EE",colnames(df_u)) &
                        grepl("freq",colnames(df_u)) &
                        grepl("2",colnames(df_u))]

# create vectors to select frequency columns that are heterozygous E
heteroE1 = colnames(df_u)[(grepl("Ee",colnames(df_u)) | grepl("eE",colnames(df_u))) &
                          grepl("freq",colnames(df_u)) &
                          grepl("1",colnames(df_u))]
heteroE2 = colnames(df_u)[(grepl("Ee",colnames(df_u)) | grepl("eE",colnames(df_u))) &
                          grepl("freq",colnames(df_u)) &
                          grepl("2",colnames(df_u))]

# calculate the proportion of the E allele in both demes
# and the difference between demes for Figure S2A
df_u$prop_E_1 = apply(df_u[,monoE1],1,sum) + 0.5 * apply(df_u[,heteroE1],1,sum)
df_u$prop_E_2 = apply(df_u[,monoE2],1,sum) + 0.5 * apply(df_u[,heteroE2],1,sum)
df_u$prop_E_dif = df_u$prop_E_1 - df_u$prop_E_2

# calculate the proportion of the E allele in both demes
# and the difference between demes for Figure S2B
df_skew$prop_E_1 = apply(df_skew[,monoE1],1,sum) + 0.5 * apply(df_skew[,heteroE1],1,sum)
df_skew$prop_E_2 = apply(df_skew[,monoE2],1,sum) + 0.5 * apply(df_skew[,heteroE2],1,sum)
df_skew$prop_E_dif = df_skew$prop_E_1 - df_skew$prop_E_2

#### Figure 2 & 3 ####

# create a custom color palette adding red to a viridis palette
viridis_pal = c("red",viridis(13, begin = 0.2, end = 0.8))

# create a base plot (p stands for plot) to build figures off of 
# by changing the data used and adding other settings
pX_base = ggplot(subset(df_u, s.e == 0.5 & r2 == 0.5 & m != 0 & skew == 1),
                 aes(x = m,
                     y = RI,
                     color = u,
                     group = u,
                     lty = u0)) +
                 scale_color_manual(values = viridis_pal) +
                 scale_linetype_manual(values = c(1,2)) +
                 geom_line() +
                 geom_point(pch = 1, size = 1) +
                 labs(color = expression(mu),
                      y = "RI",
                      x = "Migration (m)") +
                 guides(lty = "none",
                        color = guide_legend(override.aes = list(size = 1.5))) +
                 theme(legend.text = element_text(size = 10),
                       legend.key.height = unit(3.5,"mm"),
                       legend.key.width = unit(5,"mm"))

# Figure 2
p2 = pX_base %+% subset(df_u, s.e == 0.5 & r2 == 0.5 & m != 0 & skew == 1) +
        # redrawing the red points and lines so they sit on top of the other colors
        geom_line(data = subset(df_u, u == 0 & s.e == 0.5 & r2 == 0.5 & m != 0 & skew == 1)) +
        geom_point(data = subset(df_u, u == 0 & s.e == 0.5 & r2 == 0.5 & m != 0 & skew == 1), pch = 1, size = 1)
# Figure 3
p3 = pX_base %+% subset(df_u %>% arrange(u), s.e == 0.01 & r2 == 0.1 & m != 0 & skew == 1) +
        geom_line(data = subset(df_u, u == 0 & s.e == 0.01 & r2 == 0.1 & m != 0 & skew == 1)) +
        geom_point(data = subset(df_u, u == 0 & s.e == 0.01 & r2 == 0.1 & m != 0 & skew == 1), pch = 1, size = 1)

ggsave("fig2.pdf",p2,width = 130,height = 100, units = "mm")
ggsave("fig3.pdf",p3,width = 130,height = 100, units = "mm")

#### Figure 4 ####

p4 = ggplot(df_skew,
            aes(x = m, y = RI,
                color = skew,
                group = skew,
                lty = as.factor(r2))) +
            scale_linetype_manual(values = c(1,2)) + 
            scale_color_viridis_d(begin = 0.2, end = 0.8) +
            geom_point(pch = 1, size = 1) +
            geom_line() +
            geom_point(data = subset(df_u, u == 0 & s.e == 0.5 & r2 == 0.5 & skew == 0), pch = 1, color = "red", size = 1) +
            geom_line(data = subset(df_u, u == 0 & s.e == 0.5 & r2 == 0.5 & skew == 0), color = "red") +
            geom_point(data = subset(df_u, u == 0 & s.e == 0.5 & r2 == 0.001 & skew == 0), pch = 1, color = "red", size = 1) +
            geom_line(data = subset(df_u, u == 0 & s.e == 0.5 & r2 == 0.001 & skew == 0), color = "red") +
            labs(color = expression(varphi),
                 y = "RI",
                 x = "Migration (m)",
                 lty = "r") +
            guides(color = guide_legend(override.aes = list(size = 1.5)),
                   linetype = guide_legend(override.aes = list(color = "black"))) +
            theme(legend.text = element_text(size = 10),
                  legend.key.height = unit(3.5,"mm"),
                  legend.key.width = unit(5,"mm"))

ggsave("fig4.pdf",p4,width = 130,height = 100, units = "mm")

#### Figure 5 ####

# create vectors to select columns recording  
# the frequency of B over time
select_B1 = colnames(df_t)[grepl("prop_B_1_gen_",colnames(df_t))]
select_B2 = colnames(df_t)[grepl("prop_B_2_gen_",colnames(df_t))]

# save total frequency of B in the columns for deme 1
# this line cannot be run multiple times
df_t[,select_B1] = df_t[,select_B1] + df_t[,select_B2]
# create vector to select all columns but those recording deme 2
# this is important to avoid duplicating thousands of columns
# thousands of times when gathering the dataframe
filt = colnames(df_t) %in% select_B2

# convert df_t to long form
df_t_long = gather(df_t[,!filt],gen,prop_B,all_of(select_B1))
# convert generation column names into a numeric
df_t_long$gen = gsub("prop_B_1_gen_","",df_t_long$gen)
df_t_long$gen = gsub("e.","e",df_t_long$gen)
df_t_long$gen = as.numeric(df_t_long$gen)
df_t_long$gen = df_t_long$gen - min(df_t_long$gen) + 1 # remove the burn-in generations from the count

p5 = ggplot(subset(df_t_long, m == 0.1 & r2 == 0.001 & s.e == 0.5 & skew == 1),
            aes(x = log(gen),
                y = 1-prop_B/m,
                color = u,
                group = u,
                lty = u0)) +
    scale_color_manual(values = viridis_pal) +
    scale_linetype_manual(values = c(1,2)) +
    geom_line() +
    geom_point(pch = 1, size = 1) +
    geom_line(data = subset(df_t_long, u == 0 & m == 0.1 & r2 == 0.001 & s.e == 0.5 & skew == 1)) +
    geom_point(data = subset(df_t_long, u == 0 & m == 0.1 & r2 == 0.001 & s.e == 0.5 & skew == 1), pch = 1, size = 1) +
    xlim(c(0,5)) + # just show the generations where significant RI occurs
    labs(x = "log(Generation)",
         y = "RI",
         color = expression(mu)) +
    guides(lty = "none",
           color = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.text = element_text(size = 10),
          legend.key.height = unit(3.5,"mm"),
          legend.key.width = unit(5,"mm"))

ggsave("fig5.pdf",p5,width = 130,height = 100, units = "mm")

#### Supplemental Figures ####

#### Figure S2 ####

ps2.a = ggplot(subset(df_u, s.e == 0.5 & skew %in% c(0,1) & r2 == 0.5),
               aes(x = m,
                   y = prop_E_dif,
                   color = u,
                   group = u:skew,
                   lty = skew)) +
               scale_linetype_manual(values = c(1,2)) +
               scale_color_manual(values = viridis_pal) +
               geom_line() +
               geom_point(pch = 1, size = 1) +
               geom_line(data = subset(df_u, u == 0 & s.e == 0.5 & skew == 0 & r2 == 0.5)) +
               geom_point(data = subset(df_u, u == 0 & s.e == 0.5 & skew == 0 & r2 == 0.5), pch = 1, size = 1) +
               labs(color = expression(mu),
                    y = "Epigenetic divergence\nat equilibrium",
                    x = "Migration (m)",
                    lty = expression(varphi)) +
               guides(linetype = guide_legend(override.aes = list(color = "black")))

ps2.b = ggplot(subset(df_skew, s.e == 0.5 & u == 0.5),
               aes(x = m,
                   y = prop_E_dif,
                   color = skew,
                   group = skew,
                   lty = u)) +
               scale_linetype_manual(values = c(1,2)) +
               scale_color_viridis_d(begin = 0.2, end = 0.8) +
               geom_point(pch = 1, size = 1) +
               geom_line() +
               geom_line(data = subset(df_u, u == 0 & s.e == 0.5 & skew == 0 & r2 == 0.5), color = "red") +
               geom_point(data = subset(df_u, u == 0 & s.e == 0.5 & skew == 0 & r2 == 0.5), pch = 1, color = "red", size = 1) +
               labs(color = expression(varphi),
                    lty = expression(mu),
                    y = "Epigenetic divergence\nat equilibrium",
                    x = "Migration (m)") +
               guides(linetype = guide_legend(override.aes = list(color = "black")))

ps2 = plot_grid(ps2.a,ps2.b, ncol = 2, labels = "auto")

ggsave("figS2.pdf",ps2,width = 280,height = 130, units = "mm")

#### Figure S3 ####

# create title variables for reuse in plots
title_psX.a = "r = 0.001; s = 0.01"
title_psX.b = "r = 0.001; s = 0.1"
title_psX.c = "r = 0.001; s = 0.5"
title_psX.d = "r = 0.1; s = 0.01"
title_psX.e = "r = 0.1; s = 0.1"
title_psX.f = "r = 0.1; s = 0.5"
title_psX.g = "r = 0.5; s = 0.01"
title_psX.h = "r = 0.5; s = 0.1"
title_psX.i = "r = 0.5; s = 0.5"

# data subset for sub-plots
data_ps3.a = subset(df_u, skew == 1 & r2 == 0.001 & s.e == 0.01)
data_ps3.b = subset(df_u, skew == 1 & r2 == 0.001 & s.e == 0.1)
data_ps3.c = subset(df_u, skew == 1 & r2 == 0.001 & s.e == 0.5)
data_ps3.d = subset(df_u, skew == 1 & r2 == 0.1 & s.e == 0.01)
data_ps3.e = subset(df_u, skew == 1 & r2 == 0.1 & s.e == 0.1)
data_ps3.f = subset(df_u, skew == 1 & r2 == 0.1 & s.e == 0.5)
data_ps3.g = subset(df_u, skew == 1 & r2 == 0.5 & s.e == 0.01)
data_ps3.h = subset(df_u, skew == 1 & r2 == 0.5 & s.e == 0.1)
data_ps3.i = subset(df_u, skew == 1 & r2 == 0.5 & s.e == 0.5)

# create a base plot to build plots off of
psX_base = ggplot(data_ps3.a,
                  aes(x = m,
                      y = RI,
                      color = u,
                      group = u,
                      lty = u0)) +
    scale_color_manual(values = viridis_pal) +
    scale_linetype_manual(values = c(1,2)) +
    geom_line() +
    geom_point(pch = 1, size = 1) +
    geom_hline(yintercept = 0, lty = 2) +
    guides(lty = "none") +
    theme(legend.position = "none",
          axis.title.y = element_blank(),
          axis.title.x = element_blank())

# make individual sub-plots
ps3.a = psX_base %+% data_ps3.a +
    # adding red lines and points on top
    geom_line(data = subset(data_ps3.a, u == 0)) +
    geom_point(data = subset(data_ps3.a, u == 0), pch = 1, size = 1) +
    labs(title = title_psX.a)
ps3.b = psX_base %+% data_ps3.b +
    geom_line(data = subset(data_ps3.b, u == 0)) +
    geom_point(data = subset(data_ps3.b, u == 0), pch = 1, size = 1) +
    labs(title = title_psX.b)
ps3.c = psX_base %+% data_ps3.c +
    geom_line(data = subset(data_ps3.c, u == 0)) +
    geom_point(data = subset(data_ps3.c, u == 0), pch = 1, size = 1) +
    labs(title = title_psX.c)
ps3.d = psX_base %+% data_ps3.d +
    geom_line(data = subset(data_ps3.d, u == 0)) +
    geom_point(data = subset(data_ps3.d, u == 0), pch = 1, size = 1) +
    labs(title = title_psX.d) + theme(axis.title.y = element_text()) +
    labs(y = "RI") # adding y axis label to middle left plot
ps3.e = psX_base %+% data_ps3.e +
    geom_line(data = subset(data_ps3.e, u == 0)) +
    geom_point(data = subset(data_ps3.e, u == 0), pch = 1, size = 1) +
    labs(title = title_psX.e)
ps3.f = psX_base %+% data_ps3.f +
    geom_line(data = subset(data_ps3.f, u == 0)) +
    geom_point(data = subset(data_ps3.f, u == 0), pch = 1, size = 1) +
    labs(title = title_psX.f)
ps3.g = psX_base %+% data_ps3.g +
    geom_line(data = subset(data_ps3.g, u == 0)) +
    geom_point(data = subset(data_ps3.g, u == 0), pch = 1, size = 1) +
    labs(title = title_psX.g)
ps3.h = psX_base %+% data_ps3.h +
    geom_line(data = subset(data_ps3.h, u == 0)) +
    geom_point(data = subset(data_ps3.h, u == 0), pch = 1, size = 1) +
    labs(title = title_psX.h) + theme(axis.title.x = element_text()) +
    labs(x = "Migraiton (m)") # adding x axis label to middle bottom plot
ps3.i = psX_base %+% data_ps3.i +
    geom_line(data = subset(data_ps3.i, u == 0)) +
    geom_point(data = subset(data_ps3.i, u == 0), pch = 1, size = 1) +
    labs(title = title_psX.i)

# plotting subplots together
ps3_noleg = plot_grid(ps3.a,ps3.b,ps3.c,
                      ps3.d,ps3.e,ps3.f,
                      ps3.g,ps3.h,ps3.i,
                      ncol = 3, labels = "auto", align = "hv")

# getting legend from the base plot
ps3_leg <- get_legend(psX_base +
                      theme(legend.position = "right") +
                      labs(color = expression(mu)))

# plotting the subplots together with the legend
ps3_notitle = plot_grid(ps3_noleg,ps3_leg, rel_widths = c(1,0.1))

# making a title
ps3_title = ggdraw() + 
                draw_label(expression("Adaptive Epimutation ("*varphi==1*")"),
                           fontface="bold",
                           x = 0.01,
                           hjust = 0,
                           size = 20)

# plotting the title together with the rest of the plot
ps3 = plot_grid(ps3_title, ps3_notitle, ncol = 1, rel_heights = c(0.1,1))

ggsave("figS3.pdf",ps3,width = 280,height = 240, units = "mm")

#### Figure S4 ####

data_ps4.a = subset(df_u, skew == 0.5 & r2 == 0.001 & s.e == 0.01)
data_ps4.b = subset(df_u, skew == 0.5 & r2 == 0.001 & s.e == 0.1)
data_ps4.c = subset(df_u, skew == 0.5 & r2 == 0.001 & s.e == 0.5)
data_ps4.d = subset(df_u, skew == 0.5 & r2 == 0.1 & s.e == 0.01)
data_ps4.e = subset(df_u, skew == 0.5 & r2 == 0.1 & s.e == 0.1)
data_ps4.f = subset(df_u, skew == 0.5 & r2 == 0.1 & s.e == 0.5)
data_ps4.g = subset(df_u, skew == 0.5 & r2 == 0.5 & s.e == 0.01)
data_ps4.h = subset(df_u, skew == 0.5 & r2 == 0.5 & s.e == 0.1)
data_ps4.i = subset(df_u, skew == 0.5 & r2 == 0.5 & s.e == 0.5)

# changing out data from ps3 to make ps4
# since skew does not affect the genetic locus,
# the red lines and points don't have to be replotted
# with the new data
ps4.a = ps3.a %+% data_ps4.a
ps4.a = ps3.a %+% data_ps4.a
ps4.b = ps3.b %+% data_ps4.b
ps4.c = ps3.c %+% data_ps4.c
ps4.d = ps3.d %+% data_ps4.d
ps4.e = ps3.e %+% data_ps4.e
ps4.f = ps3.f %+% data_ps4.f
ps4.g = ps3.g %+% data_ps4.g
ps4.h = ps3.h %+% data_ps4.h
ps4.i = ps3.i %+% data_ps4.i

# same as before
ps4_noleg = plot_grid(ps4.a,ps4.b,ps4.c,
                      ps4.d,ps4.e,ps4.f,
                      ps4.g,ps4.h,ps4.i,
                      ncol = 3, labels = "auto", align = "hv")

ps4_leg <- get_legend(psX_base +
                          theme(legend.position = "right") +
                          labs(color = expression(mu)))

ps4_notitle = plot_grid(ps4_noleg,ps4_leg, rel_widths = c(1,0.1))

ps4_title = ggdraw() + 
    draw_label(expression("Partially-Adaptive Epimutation ("*varphi==0.5*")"),
               fontface="bold",
               x = 0.01,
               hjust = 0,
               size = 20)

ps4 = plot_grid(ps4_title, ps4_notitle, ncol = 1, rel_heights = c(0.1,1))

ggsave("figS4.pdf",ps4,width = 280,height = 240, units = "mm")


#### Figure S5 ####

data_ps5.a = subset(df_u, skew == 0 & r2 == 0.001 & s.e == 0.01)
data_ps5.b = subset(df_u, skew == 0 & r2 == 0.001 & s.e == 0.1)
data_ps5.c = subset(df_u, skew == 0 & r2 == 0.001 & s.e == 0.5)
data_ps5.d = subset(df_u, skew == 0 & r2 == 0.1 & s.e == 0.01)
data_ps5.e = subset(df_u, skew == 0 & r2 == 0.1 & s.e == 0.1)
data_ps5.f = subset(df_u, skew == 0 & r2 == 0.1 & s.e == 0.5)
data_ps5.g = subset(df_u, skew == 0 & r2 == 0.5 & s.e == 0.01)
data_ps5.h = subset(df_u, skew == 0 & r2 == 0.5 & s.e == 0.1)
data_ps5.i = subset(df_u, skew == 0 & r2 == 0.5 & s.e == 0.5)

ps5.a = ps3.a %+% data_ps5.a
ps5.a = ps3.a %+% data_ps5.a
ps5.b = ps3.b %+% data_ps5.b
ps5.c = ps3.c %+% data_ps5.c
ps5.d = ps3.d %+% data_ps5.d
ps5.e = ps3.e %+% data_ps5.e
ps5.f = ps3.f %+% data_ps5.f
ps5.g = ps3.g %+% data_ps5.g
ps5.h = ps3.h %+% data_ps5.h
ps5.i = ps3.i %+% data_ps5.i

ps5_noleg = plot_grid(ps5.a,ps5.b,ps5.c,
                      ps5.d,ps5.e,ps5.f,
                      ps5.g,ps5.h,ps5.i,
                      ncol = 3, labels = "auto", align = "hv")

ps5_leg <- get_legend(psX_base +
                          theme(legend.position = "right") +
                          labs(color = expression(mu)))

ps5_notitle = plot_grid(ps5_noleg,ps5_leg, rel_widths = c(1,0.1))

ps5_title = ggdraw() + 
    draw_label(expression("Neutral Epimutation ("*varphi==0*")"),
               fontface="bold",
               x = 0.01,
               hjust = 0,
               size = 20)

ps5 = plot_grid(ps5_title, ps5_notitle, ncol = 1, rel_heights = c(0.1,1))

ggsave("figS5.pdf",ps5,width = 280,height = 240, units = "mm")

#### Figure S6 ####

data_ps6.a = subset(df_u, skew == -1 & r2 == 0.001 & s.e == 0.01)
data_ps6.b = subset(df_u, skew == -1 & r2 == 0.001 & s.e == 0.1)
data_ps6.c = subset(df_u, skew == -1 & r2 == 0.001 & s.e == 0.5)
data_ps6.d = subset(df_u, skew == -1 & r2 == 0.1 & s.e == 0.01)
data_ps6.e = subset(df_u, skew == -1 & r2 == 0.1 & s.e == 0.1)
data_ps6.f = subset(df_u, skew == -1 & r2 == 0.1 & s.e == 0.5)
data_ps6.g = subset(df_u, skew == -1 & r2 == 0.5 & s.e == 0.01)
data_ps6.h = subset(df_u, skew == -1 & r2 == 0.5 & s.e == 0.1)
data_ps6.i = subset(df_u, skew == -1 & r2 == 0.5 & s.e == 0.5)

ps6.a = ps3.a %+% data_ps6.a
ps6.a = ps3.a %+% data_ps6.a
ps6.b = ps3.b %+% data_ps6.b
ps6.c = ps3.c %+% data_ps6.c
ps6.d = ps3.d %+% data_ps6.d
ps6.e = ps3.e %+% data_ps6.e
ps6.f = ps3.f %+% data_ps6.f
ps6.g = ps3.g %+% data_ps6.g
ps6.h = ps3.h %+% data_ps6.h
ps6.i = ps3.i %+% data_ps6.i

ps6_noleg = plot_grid(ps6.a,ps6.b,ps6.c,
                      ps6.d,ps6.e,ps6.f,
                      ps6.g,ps6.h,ps6.i,
                      ncol = 3, labels = "auto", align = "hv")

ps6_leg <- get_legend(psX_base +
                          theme(legend.position = "right") +
                          labs(color = expression(mu)))

ps6_notitle = plot_grid(ps6_noleg,ps6_leg, rel_widths = c(1,0.1))

ps6_title = ggdraw() + 
    draw_label(expression("Maladaptive Epimutation ("*varphi==-1*")"),
               fontface="bold",
               x = 0.01,
               hjust = 0,
               size = 20)

ps6 = plot_grid(ps6_title, ps6_notitle, ncol = 1, rel_heights = c(0.1,1))

ggsave("figS6.pdf",ps6,width = 280,height = 240, units = "mm")

## Figure S7 ##

data_ps7.a = subset(df_u_pre, skew == 1 & r2 == 0.001 & s.e == 0.01)
data_ps7.b = subset(df_u_pre, skew == 1 & r2 == 0.001 & s.e == 0.1)
data_ps7.c = subset(df_u_pre, skew == 1 & r2 == 0.001 & s.e == 0.5)
data_ps7.d = subset(df_u_pre, skew == 1 & r2 == 0.1 & s.e == 0.01)
data_ps7.e = subset(df_u_pre, skew == 1 & r2 == 0.1 & s.e == 0.1)
data_ps7.f = subset(df_u_pre, skew == 1 & r2 == 0.1 & s.e == 0.5)
data_ps7.g = subset(df_u_pre, skew == 1 & r2 == 0.5 & s.e == 0.01)
data_ps7.h = subset(df_u_pre, skew == 1 & r2 == 0.5 & s.e == 0.1)
data_ps7.i = subset(df_u_pre, skew == 1 & r2 == 0.5 & s.e == 0.5)

ps7.a = ps3.a %+% data_ps7.a
ps7.a = ps3.a %+% data_ps7.a
ps7.b = ps3.b %+% data_ps7.b
ps7.c = ps3.c %+% data_ps7.c
ps7.d = ps3.d %+% data_ps7.d
ps7.e = ps3.e %+% data_ps7.e
ps7.f = ps3.f %+% data_ps7.f
ps7.g = ps3.g %+% data_ps7.g
ps7.h = ps3.h %+% data_ps7.h
ps7.i = ps3.i %+% data_ps7.i

ps7_noleg = plot_grid(ps7.a,ps7.b,ps7.c,
                      ps7.d,ps7.e,ps7.f,
                      ps7.g,ps7.h,ps7.i,
                      ncol = 3, labels = "auto", align = "hv")

ps7_leg <- get_legend(psX_base +
                          theme(legend.position = "right") +
                          labs(color = expression(mu)))

ps7_notitle = plot_grid(ps7_noleg,ps7_leg, rel_widths = c(1,0.1))

ps7_title = ggdraw() + 
    draw_label(expression("Pre-selection Epimutation ("*varphi==1*")"),
               fontface="bold",
               x = 0.01,
               hjust = 0,
               size = 20)

ps7 = plot_grid(ps7_title, ps7_notitle, ncol = 1, rel_heights = c(0.1,1))

ggsave("figS7.pdf",ps7,width = 280,height = 240, units = "mm")

#### Figure S8 ####

# get the maximum RI for each m, s and r value combination
# get the epimutation rate which caused that maximal RI
subset(df_u,skew == 1 & m != 0) %>%
    group_by(m,s.e,r2) %>%
    summarise(max_RI = max(RI),
              max_u = as.numeric(as.character(u))[RI == max_RI]) -> df_u_sum

# get the greatest RI within each subplot for y-axis scaling
# of subplots
df_u_sum %>%
    group_by(s.e,r2) %>%
    summarise(max_max_RI = max(max_RI),
              max_max_u = max(max_u)) -> df_u_sum_sum

# add the max_max_RI to the original summary data.frame
for(i in unique(df_u_sum_sum$s.e)) {
    for(j in unique(df_u_sum_sum$r2)) {
        df_u_sum$max_max_RI[df_u_sum$s.e == i & df_u_sum$r2 == j] = df_u_sum_sum$max_max_RI[df_u_sum_sum$s.e == i & df_u_sum_sum$r2 == j]
    }
}

# rescale max_mu values by the max_max_RI within a subplot
# this rescaling will be undone by rescaling the y-axis
# with the opposite transformation
df_u_sum$max_u_trans = 2*df_u_sum$max_u*df_u_sum$max_max_RI
df_u_sum$max_mu = df_u_sum$max_u_trans # just to rename it

# get max_RI and max_mu into a single column
df_u_sum_long = gather(df_u_sum, measure, value, max_RI, max_mu)

# subset for subplots
data_ps8.a = subset(df_u_sum_long, r2 == 0.001 & s.e == 0.01)
data_ps8.b = subset(df_u_sum_long, r2 == 0.001 & s.e == 0.1)
data_ps8.c = subset(df_u_sum_long, r2 == 0.001 & s.e == 0.5)
data_ps8.d = subset(df_u_sum_long, r2 == 0.1 & s.e == 0.01)
data_ps8.e = subset(df_u_sum_long, r2 == 0.1 & s.e == 0.1)
data_ps8.f = subset(df_u_sum_long, r2 == 0.1 & s.e == 0.5)
data_ps8.g = subset(df_u_sum_long, r2 == 0.5 & s.e == 0.01)
data_ps8.h = subset(df_u_sum_long, r2 == 0.5 & s.e == 0.1)
data_ps8.i = subset(df_u_sum_long, r2 == 0.5 & s.e == 0.5)

# get a subplots max_max_RI as a variable
max_ps8.a = data_ps8.a$max_max_RI[1]
max_ps8.b = data_ps8.b$max_max_RI[1]
max_ps8.c = data_ps8.c$max_max_RI[1]
max_ps8.d = data_ps8.d$max_max_RI[1]
max_ps8.e = data_ps8.e$max_max_RI[1]
max_ps8.f = data_ps8.f$max_max_RI[1]
max_ps8.g = data_ps8.g$max_max_RI[1]
max_ps8.h = data_ps8.h$max_max_RI[1]
max_ps8.i = data_ps8.i$max_max_RI[1]

# create base plot to build other plots on
ps8_base = ggplot(data_ps8.a,
                   aes(x = m,
                       y = value,
                       group = measure,
                       pch = measure,
                       lty = measure)) +
                   geom_point() +
                   geom_line() +
                   theme(axis.title.x=element_blank(),
                         axis.title.y=element_blank())

ps8.a = ps8_base %+% data_ps8.a +
    theme(legend.position = "none") +
    labs(title = title_psX.a) +
    # this is where we resize the secondary y-axis to get back
    # the original max_mu values, this is done so that we can
    # scale each subplots y axis independently while having
    # two y-axes
    scale_y_continuous(sec.axis = sec_axis(~0.5*./max_ps8.a))
ps8.b = ps8_base %+% data_ps8.b +
    theme(legend.position = "none") +
    labs(title = title_psX.b) +
    scale_y_continuous(sec.axis = sec_axis(~0.5*./max_ps8.b))
ps8.c = ps8_base %+% data_ps8.c +
    theme(legend.position = "none") +
    labs(title = title_psX.c) +
    scale_y_continuous(sec.axis = sec_axis(~0.5*./max_ps8.c))
ps8.d = ps8_base %+% data_ps8.d +
    theme(legend.position = "none") +
    labs(title = title_psX.d) +
    theme(axis.title.y = element_text()) +
    scale_y_continuous(name = "Maximum possible RI", sec.axis = sec_axis(~0.5*./max_ps8.d))
ps8.e = ps8_base %+% data_ps8.e +
    theme(legend.position = "none") +
    labs(title = title_psX.e) +
    scale_y_continuous(sec.axis = sec_axis(~0.5*./max_ps8.e))
ps8.f = ps8_base %+% data_ps8.f +
    theme(legend.position = "none") +
    labs(title = title_psX.f) +
    theme(axis.title.y = element_text()) +
    scale_y_continuous(name = "", sec.axis = sec_axis(~0.5*./max_ps8.f, name = "RI maximizing epimutation rate"))
ps8.g = ps8_base %+% data_ps8.g +
    theme(legend.position = "none") +
    labs(title = title_psX.g) +
    scale_y_continuous(sec.axis = sec_axis(~0.5*./max_ps8.g))
ps8.h = ps8_base %+% data_ps8.h +
    theme(legend.position = "none") +
    theme(axis.title.x = element_text()) +
    labs(title = title_psX.h, x = "Migration (m)") +
    scale_y_continuous(sec.axis = sec_axis(~0.5*./max_ps8.h))
ps8.i = ps8_base %+% data_ps8.i +
    theme(legend.position = "none") +
    labs(title = title_psX.i) +
    scale_y_continuous(sec.axis = sec_axis(~0.5*./max_ps8.i))

ps8_noleg = plot_grid(ps8.a,ps8.b,ps8.c,
                      ps8.d,ps8.e,ps8.f,
                      ps8.g,ps8.h,ps8.i,
                      ncol = 3, labels = "auto", align = "hv")

ps8_leg <- get_legend(ps8_base)

ps8_notitle = plot_grid(ps8_noleg,ps8_leg, rel_widths = c(1,0.1))

ps8_title = ggdraw() + 
    draw_label(expression("RI-Maximizing Epimutation ("*varphi==1*")"),
               fontface="bold",
               x = 0.01,
               hjust = 0,
               size = 20)

ps8 = plot_grid(ps8_title, ps8_notitle, ncol = 1, rel_heights = c(0.1,1))

ggsave("figS8.pdf",ps8,width = 280,height = 220, units = "mm")

