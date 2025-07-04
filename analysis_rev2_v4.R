library(ggplot2)
library(cowplot)
library(viridis)
library(dplyr)
library(tidyr)
library(EigenR)
library(Matrix)
theme_set(theme_cowplot())

# change this to the your directory containing the scripts
setwd("/Users/anon/Desktop/MS5/drafts/proc_B_rev2/rev2_code")

#### Data prepartation function ####

add_calculated_cols <- function(df) {
    
    # recalculate epimutaiton rate from i and c parameters
    df$u = (df$i_1 + df$c_1) / 2
    # division by small numbers causes slight differences in u values
    # that should be the same, rounding removes these small differences
    df$u = round(df$u, digits = 8)
    df$u_fact = factor(df$u, levels = sort(unique(df$u), decreasing = F)[c(2:length(unique(df$u)),1)]) # make factor for plotting
    # df$u_fact2 = factor(df$u_fact, levels = levels(df$u_fact)[c(2:11,1,12)])
    
    # recalculate skew from i and c parameters
    df$skew = (df$i_1 - df$c_1) / (df$i_1 + df$c_1)
    # division by small numbers causes slight differences in skew values
    # that should be the same, rounding removes these small differences
    df$skew = round(df$skew, digits = 2) 
    
    # fill in skew values when mu = 0, even though skew doesn't affect
    # the results of these simulations this makes plotting easier later
    uniq_skew = unique(df$skew)[!is.nan(unique(df$skew))]
    nan_tot = sum(is.nan(df$skew))
    df$skew[is.nan(df$skew)] = rep(uniq_skew, each = nan_tot/length(uniq_skew))
    df$skew_fact = factor(df$skew, levels = sort(unique(df$skew), decreasing = F)) # make factor for plotting
    
    # create vector marking when mu = 0 for plotting the genetic locus
    df$u0 = F
    df$u0[df$u == 0] = T
    df$u0[df$u == 0.5] = T
    
    # create vectors to select frequency columns that are homozygous E
    monoE1 = colnames(df)[grepl("EE", colnames(df)) &
                              grepl("freq", colnames(df)) &
                              grepl("1", colnames(df))]
    monoE2 = colnames(df)[grepl("EE", colnames(df)) &
                              grepl("freq", colnames(df)) &
                              grepl("2", colnames(df))]
    
    # create vectors to select frequency columns that are heterozygous E
    heteroE1 = colnames(df)[(grepl("Ee", colnames(df)) | grepl("eE", colnames(df))) &
                                grepl("freq", colnames(df)) &
                                grepl("1", colnames(df))]
    heteroE2 = colnames(df)[(grepl("Ee", colnames(df)) | grepl("eE", colnames(df))) &
                                grepl("freq", colnames(df)) &
                                grepl("2", colnames(df))]
    
    # calculate the proportion of the E allele in both demes
    # and the difference between demes for Figure S2A
    df$prop_E_1 = apply(df[, monoE1], 1, sum) + 0.5 * apply(df[, heteroE1], 1, sum)
    df$prop_E_2 = apply(df[, monoE2], 1, sum) + 0.5 * apply(df[, heteroE2], 1, sum)
    df$prop_E_dif = df$prop_E_1 - df$prop_E_2
    
    # create vectors to select frequency columns that are homozygous A
    monoA1 = colnames(df)[grepl("AA", colnames(df)) &
                              grepl("freq", colnames(df)) &
                              grepl("1", colnames(df))]
    monoA2 = colnames(df)[grepl("AA", colnames(df)) &
                              grepl("freq", colnames(df)) &
                              grepl("2", colnames(df))]
    
    # create vectors to select frequency columns that are heterozygous A
    heteroA1 = colnames(df)[(grepl("Aa", colnames(df)) | grepl("aA", colnames(df))) &
                                grepl("freq", colnames(df)) &
                                grepl("1", colnames(df))]
    heteroA2 = colnames(df)[(grepl("Aa", colnames(df)) | grepl("aA", colnames(df))) &
                                grepl("freq", colnames(df)) &
                                grepl("2", colnames(df))]
    
    # calculate the proportion of the A allele in both demes
    # and the difference between demes for Figure S2A
    df$prop_A_1 = apply(df[, monoA1], 1, sum) + 0.5 * apply(df[, heteroA1], 1, sum)
    df$prop_A_2 = apply(df[, monoA2], 1, sum) + 0.5 * apply(df[, heteroA2], 1, sum)
    df$prop_A_dif = df$prop_A_1 - df$prop_A_2
    
    return(df)
}


#### Load data ####

# simulation varying epimutation rate
# default
df_u = read.csv("RI_sim2.6.3_2025-03-17_12-11-22.csv")
df_u = add_calculated_cols(df_u)

# simulation varying epimutation rate
# with secondary contact
df_u_2nd = read.csv("RI_sim2.6.3_2025-03-17_12-15-25.csv")
df_u_2nd = add_calculated_cols(df_u_2nd)

# simulation varying epimutation rate
# with pre-selection epimutation
df_u_pre = read.csv("RI_sim2.6.3_2025-03-17_12-20-03.csv")
df_u_pre = add_calculated_cols(df_u_pre)

# simulation varying skew
df_skew = read.csv("RI_sim2.6.3_2025-03-17_12-23-42.csv")
df_skew = add_calculated_cols(df_skew)

# simulation varying selection
df_s = read.csv("RI_sim2.6.3_2025-03-17_12-24-21.csv")
df_s = add_calculated_cols(df_s)

# simulation fine gradations of recombination
df_r = read.csv("RI_sim2.6.3_2025-03-17_15-50-38.csv")
df_r = add_calculated_cols(df_r)

# simulations of LD over time
df_LD = read.csv("RI_sim2.6.3_2025-03-17_15-45-04.csv")
colnames(df_LD)[colnames(df_LD) == "Dnorm_eB_3_gen_2e.06"] = "Dnorm_eB_3_gen_2000000"
df_LD = add_calculated_cols(df_LD)

### simulations of different periods of allopatry

# make a list of .csv files
file_list <- c(
    "RI_sim2.6.3_2025-03-17_16-15-25.csv",
    "RI_sim2.6.3_2025-03-17_16-16-00.csv",
    "RI_sim2.6.3_2025-03-17_16-16-35.csv",
    "RI_sim2.6.3_2025-03-17_16-17-14.csv",
    "RI_sim2.6.3_2025-03-17_16-17-55.csv",
    "RI_sim2.6.3_2025-03-17_16-18-31.csv",
    "RI_sim2.6.3_2025-03-17_16-19-14.csv",
    "RI_sim2.6.3_2025-03-17_16-19-54.csv",
    "RI_sim2.6.3_2025-03-17_16-20-48.csv",
    "RI_sim2.6.3_2025-03-17_16-21-38.csv",
    "RI_sim2.6.3_2025-03-17_16-22-19.csv",
    "RI_sim2.6.3_2025-03-17_16-23-16.csv",
    "RI_sim2.6.3_2025-03-17_16-24-06.csv",
    "RI_sim2.6.3_2025-03-17_16-24-59.csv",
    "RI_sim2.6.3_2025-03-17_16-25-50.csv",
    "RI_sim2.6.3_2025-03-17_16-26-47.csv",
    "RI_sim2.6.3_2025-03-17_16-27-33.csv",
    "RI_sim2.6.3_2025-03-17_16-28-13.csv",
    "RI_sim2.6.3_2025-03-17_16-28-51.csv",
    "RI_sim2.6.3_2025-03-17_16-29-38.csv",
    "RI_sim2.6.3_2025-03-17_16-30-25.csv",
    "RI_sim2.6.3_2025-03-17_16-31-13.csv",
    "RI_sim2.6.3_2025-03-17_16-32-20.csv",
    "RI_sim2.6.3_2025-03-17_16-33-37.csv",
    "RI_sim2.6.3_2025-03-17_16-34-59.csv",
    "RI_sim2.6.3_2025-03-17_16-36-36.csv",
    "RI_sim2.6.3_2025-03-17_16-38-01.csv",
    "RI_sim2.6.3_2025-03-17_16-39-10.csv",
    "RI_sim2.6.3_2025-03-17_16-39-59.csv",
    "RI_sim2.6.3_2025-03-17_16-40-44.csv",
    "RI_sim2.6.3_2025-03-17_16-41-55.csv",
    "RI_sim2.6.3_2025-03-17_16-42-57.csv",
    "RI_sim2.6.3_2025-03-17_16-43-46.csv",
    "RI_sim2.6.3_2025-03-17_16-44-38.csv",
    "RI_sim2.6.3_2025-03-17_16-46-05.csv",
    "RI_sim2.6.3_2025-03-17_16-48-25.csv",
    "RI_sim2.6.3_2025-03-17_16-54-43.csv",
    "RI_sim2.6.3_2025-03-17_16-58-47.csv",
    "RI_sim2.6.3_2025-03-17_17-02-26.csv",
    "RI_sim2.6.3_2025-03-17_17-05-41.csv",
    "RI_sim2.6.3_2025-03-17_17-08-55.csv",
    "RI_sim2.6.3_2025-03-17_17-12-41.csv",
    "RI_sim2.6.3_2025-03-17_17-14-37.csv",
    "RI_sim2.6.3_2025-03-17_17-16-15.csv",
    "RI_sim2.6.3_2025-03-17_17-17-49.csv",
    "RI_sim2.6.3_2025-03-17_17-19-22.csv"
)

# Create the allo sequence: 1 to 20 by 1, then 25 to 100 by 10
allo_values <- c(seq(0, 20, by = 1), seq(25, 95, by = 5), seq(100, 1000, by = 100))
# allo_values <- c(seq(0, 20, by = 1), seq(25, 50, by = 5))

# Read all CSV files into a list of data frames, apply the function and add 'allo' column
df_list <- lapply(seq_along(file_list), function(i) {
    df <- read.csv(file_list[i])
    df <- add_calculated_cols(df)  # Apply the add_calculated_cols function
    df$allo <- paste0("df", allo_values[i])  # Add the 'allo' column
    df$allo_num <- allo_values[i]  # Add the 'allo_num' column
    return(df)
})

# Combine all data frames into one
df_all <- do.call(rbind, df_list)

#### Custom palettes for coloring plots ####

viridis_pal_p3 = c("black",viridis(5, begin = 0.2, end = 0.8),"red")

viridis_pal_p2 = viridis(5, begin = 0.2, end = 0.8)

viridis_pal_ps2 = c(viridis(5, begin = 0.2, end = 0.8),
                    magma(5, begin = 0.1, end = 0.8))

#### facet labels for plots ####

# Custom labeller function for facets
facet_labels <- labeller(
    `as.factor(r2):as.factor(s.e)` = function(labels) {
        sapply(labels, function(label) {
            s_e <- strsplit(label, ":")[[1]][2]  # Extract s.e
            r2 <- strsplit(label, ":")[[1]][1]  # Extract r2
            paste0("r = ", r2, "; s = ", s_e)
        })
    }
)

# ps9
# Create a custom labeller function to parse the expressions
facet_labels_ps9 <- as_labeller(c(
    `0.01:0.001` = "s == 0.01~';' ~m == 0.001",
    `0.01:0.1` = "s == 0.01~',' ~m == 0.1",
    `0.01:0.5` = "s == 0.01~';' ~m == 0.5",
    `0.1:0.001` = "s == 0.1~';' ~m == 0.001",
    `0.1:0.1` = "s == 0.1~';' ~m == 0.1",
    `0.1:0.5` = "s == 0.1~';' ~m == 0.5",
    `0.5:0.001` = "s == 0.5~';' ~m == 0.001",
    `0.5:0.1` = "s == 0.5~';' ~m == 0.1",
    `0.5:0.5` = "s == 0.5~';' ~m == 0.5"
), label_parsed)

# p4a
facet_labels_p4a <- as_labeller(c(
    `0.5` = "(A)~'   '~m == 0.5*'  '",
    `0.1` = "(B)~'   '~m == 0.1*'  '",
    `0.001` = "(C)~' '~m == 0.001"
), label_parsed)

# p4b
facet_labels_p4b <- as_labeller(c(
    `0.5` = "(D)~'   '~m == 0.5*'  '",
    `0.1` = "(E)~'   '~m == 0.1*'  '",
    `0.001` = "(F)~' '~m == 0.001"
), label_parsed)

# p5
facet_labels_p5 <- as_labeller(c(
    `0.001` = "(A)~' '~r == 0.001",
    `0.5` = "(B)~'   '~r == 0.5*'  '"
), label_parsed)

# ps2 ps3
# Custom labeller function for facets
facet_labels_ps2 <- labeller(
    `as.factor(s.e):as.factor(m)` = function(labels) {
        sapply(labels, function(label) {
            s_e <- strsplit(label, ":")[[1]][1]  # Extract s.e
            m <- strsplit(label, ":")[[1]][2]  # Extract r2
            paste0("m = ", m, "; s = ", s_e)
        })
    }
)


# ps4 ps5
# Create a custom labeller function to parse the expressions
facet_labels_ps4 <- as_labeller(c(
    `0.01:0.001` = "s == 0.01~';' ~r == 0.001",
    `0.01:0.1` = "s == 0.01~',' ~r == 0.1",
    `0.01:0.5` = "s == 0.01~';' ~r == 0.5",
    `0.1:0.001` = "s == 0.1~';' ~r == 0.001",
    `0.1:0.1` = "s == 0.1~';' ~r == 0.1",
    `0.1:0.5` = "s == 0.1~';' ~r == 0.1",
    `0.5:0.001` = "s == 0.5~';' ~r == 0.001",
    `0.5:0.1` = "s == 0.5~';' ~r == 0.1",
    `0.5:0.5` = "s == 0.5~';' ~r == 0.5"
), label_parsed)

#ps6-ps8
# Create a custom labeller function to parse the expressions
facet_labels_ps6 <- as_labeller(c(
    `0.01:0` = "s == 0.01~';' ~varphi == 0",
    `0.01:0.5` = "s == 0.01~',' ~varphi == 0.5",
    `0.01:1` = "s == 0.01~';' ~varphi == 1",
    `0.1:0` = "s == 0.1~';' ~varphi == 0",
    `0.1:0.5` = "s == 0.1~';' ~varphi == 0.5",
    `0.1:1` = "s == 0.1~';' ~varphi == 1",
    `0.5:0` = "s == 0.5~';' ~varphi == 0",
    `0.5:0.5` = "s == 0.5~';' ~varphi == 0.5",
    `0.5:1` = "s == 0.5~';' ~varphi == 1"
), label_parsed)

#### Figures S10 to S16 ####

# Base plot
psX_base <- ggplot(subset(df_u, skew == 1),
                   aes(x = m,
                       y = RI,
                       color = u_fact,
                       group = u_fact,
                       lty = u0,
                       size = u0)) +
    scale_color_manual(values = viridis_pal_p3,
                       breaks = c(0,0.02,0.1,0.2,0.3,0.4,0.5)) +
    scale_linetype_manual(values = c("solid","12")) +
    scale_size_manual(values = c(0.5,1.3), guide = "none") +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    geom_line() +
    guides(lty = "none",
           color = guide_legend(
               override.aes = list(
                   linetype = c("12", "solid", "solid", "solid", "solid", "solid", "12"), 
                   linewidth = c(1.3, 0.5, 0.5, 0.5, 0.5, 0.5, 1.3)
               )
           )) +
    facet_wrap(~as.factor(r2):as.factor(s.e), scales = "free_y", ncol = 3, labeller = facet_labels) +
    labs(x = "Migration (m)", y = "RI", color = expression(mu),
         title = expression("Environmentally Induced Epimutation ("*varphi*" = 1)"))

psX_base_small <- ggplot(subset(df_u, skew == 1),
                         aes(x = m, y = RI, color = u_fact, group = u_fact, lty = u0)) +
    scale_color_manual(values = viridis_pal_p3, breaks = levels(df_u$u_fact)) +
    scale_linetype_manual(values = c(1, 2)) +
    geom_line(size = 0.3) +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    guides(lty = "none") +
    facet_wrap(~as.factor(r2):as.factor(s.e), scales = "free_y", ncol = 3, labeller = facet_labels) +
    labs(x = "Migration (m)", y = "RI", color = expression(mu),
         title = expression("Environmental Induction ("*varphi*" = 1)"))


ps10_base <- psX_base %+% subset(df_u, skew == 0 & s.e %in% c(0.01,0.1,0.5)) +
    labs(title = expression("Neutral Epimutation ("*varphi*" = 0)"))  # Update the title to varphi = 0

inlay_ps10_1 = psX_base %+% subset(df_u, skew == 0 & r2 == 0 & s.e == 0.01) +
    facet_null() +
    theme(title = element_blank(),
          legend.position = "none") +
    coord_cartesian(ylim = c(0,0.008)) +
    theme(axis.text = element_text(size = 8))

inlay_ps10_2 = inlay_ps10_1 %+% subset(df_u, skew == 0 & r2 == 0.001 & s.e == 0.01)
inlay_ps10_3 = inlay_ps10_1 %+% subset(df_u, skew == 0 & r2 == 0.1 & s.e == 0.01)
inlay_ps10_4 = inlay_ps10_1 %+% subset(df_u, skew == 0 & r2 == 0.5 & s.e == 0.01)

ps10 = ggdraw(ps10_base) +
    draw_plot(inlay_ps10_1, x = 0.12, y = 0.75, width = 0.2, height = 0.18) +
    draw_plot(inlay_ps10_2, x = 0.12, y = 0.52, width = 0.2, height = 0.18) +
    draw_plot(inlay_ps10_3, x = 0.12, y = 0.29, width = 0.2, height = 0.18) +
    draw_plot(inlay_ps10_4, x = 0.12, y = 0.07, width = 0.2, height = 0.18)

ps11_base <- psX_base %+% subset(df_u, skew == 1 & s.e %in% c(0.01,0.1,0.5)) +
    labs(title = expression("Environmental Induction ("*varphi*" = 1)"))  # Update the title to varphi = 1

inlay_ps11_1 = inlay_ps10_1 %+% subset(df_u, skew == 1 & r2 == 0 & s.e == 0.01) +
    coord_cartesian(ylim = c(0,0.03))
inlay_ps11_2 = inlay_ps10_1 %+% subset(df_u, skew == 1 & r2 == 0.001 & s.e == 0.01) +
    coord_cartesian(ylim = c(0,0.03))

ps11 = ggdraw(ps11_base) +
    draw_plot(inlay_ps11_1, x = 0.12, y = 0.75, width = 0.2, height = 0.18) +
    draw_plot(inlay_ps11_2, x = 0.12, y = 0.52, width = 0.2, height = 0.18)

ps12_base <- psX_base %+% subset(df_u, skew == 0.5 & s.e %in% c(0.01,0.1,0.5)) +
    labs(title = expression("Partial Environmental Induction ("*varphi*" = 0.5)"))  # Update the title to varphi = 0.5

inlay_ps12_1 = inlay_ps10_1 %+% subset(df_u, skew == 0.5 & r2 == 0 & s.e == 0.01) +
    coord_cartesian(ylim = c(0,0.02))
inlay_ps12_2 = inlay_ps10_1 %+% subset(df_u, skew == 0.5 & r2 == 0.001 & s.e == 0.01) +
    coord_cartesian(ylim = c(0,0.02))

ps12 = ggdraw(ps12_base) +
    draw_plot(inlay_ps12_1, x = 0.12, y = 0.75, width = 0.2, height = 0.18) +
    draw_plot(inlay_ps12_2, x = 0.12, y = 0.52, width = 0.2, height = 0.18)

ps13_base <- psX_base %+% subset(df_u, skew == -1 & s.e %in% c(0.01,0.1,0.5)) +
    labs(title = expression("Maladaptive Environmental Induction ("*varphi*" = -1)"))  # Update the title to varphi = -1

inlay_ps13_1 = inlay_ps10_1 %+% subset(df_u, skew == -1 & r2 == 0 & s.e == 0.01) +
    coord_cartesian(ylim = c(-0.03,0.03))
inlay_ps13_2 = inlay_ps10_1 %+% subset(df_u, skew == -1 & r2 == 0.001 & s.e == 0.01) +
    coord_cartesian(ylim = c(-0.03,0.03))

ps13 = ggdraw(ps13_base) +
    draw_plot(inlay_ps13_1, x = 0.12, y = 0.77, width = 0.2, height = 0.16) +
    draw_plot(inlay_ps13_2, x = 0.12, y = 0.54, width = 0.2, height = 0.16)

ps14_base = psX_base %+% subset(df_u_pre, skew == 1 & s.e %in% c(0.01,0.1,0.5)) +
    labs(title = expression("Pre-selection Environmental Induction ("*varphi*" = 1)")) 

inlay_ps14_1 = inlay_ps10_1 %+% subset(df_u_pre, skew == 1 & r2 == 0 & s.e == 0.01) +
    coord_cartesian(ylim = c(0,0.12))
inlay_ps14_2 = inlay_ps10_1 %+% subset(df_u_pre, skew == 1 & r2 == 0.001 & s.e == 0.01) +
    coord_cartesian(ylim = c(0,0.12))

ps14 = ggdraw(ps14_base) +
    draw_plot(inlay_ps14_1, x = 0.12, y = 0.75, width = 0.2, height = 0.18) +
    draw_plot(inlay_ps14_2, x = 0.12, y = 0.52, width = 0.2, height = 0.18)

ps15_base = psX_base %+% subset(df_u_2nd, skew == 0 & s.e %in% c(0.01,0.1,0.5)) +
    labs(title = expression("Secondary Contact, Neutral Epimutation ("*varphi*" = 0)")) 

inlay_ps15_1 = inlay_ps10_1 %+% subset(df_u_2nd, skew == 0 & r2 == 0 & s.e == 0.01) +
    coord_cartesian(ylim = c(0,0.03))
inlay_ps15_2 = inlay_ps10_1 %+% subset(df_u_2nd, skew == 0 & r2 == 0.001 & s.e == 0.01) +
    coord_cartesian(ylim = c(0,0.03))

ps15 = ggdraw(ps15_base) +
    draw_plot(inlay_ps15_1, x = 0.12, y = 0.75, width = 0.2, height = 0.18) +
    draw_plot(inlay_ps15_2, x = 0.12, y = 0.52, width = 0.2, height = 0.18)

ps16_base = psX_base %+% subset(df_u_2nd, skew == 1 & s.e %in% c(0.01,0.1,0.5)) +
    labs(title = expression("Secondary Contact, Environmental Induction ("*varphi*" = 1)")) 

inlay_ps16_1 = inlay_ps10_1 %+% subset(df_u_2nd, skew == 1 & r2 == 0 & s.e == 0.01) +
    coord_cartesian(ylim = c(0,0.12))
inlay_ps16_2 = inlay_ps10_1 %+% subset(df_u_2nd, skew == 1 & r2 == 0.001 & s.e == 0.01) +
    coord_cartesian(ylim = c(0,0.12))

ps16 = ggdraw(ps16_base) +
    draw_plot(inlay_ps16_1, x = 0.12, y = 0.75, width = 0.2, height = 0.18) +
    draw_plot(inlay_ps16_2, x = 0.12, y = 0.52, width = 0.2, height = 0.18)

ggsave("figS10.pdf",ps10,width = 240,height = 250, units = "mm")
ggsave("figS11.pdf",ps11,width = 240,height = 250, units = "mm")
ggsave("figS12.pdf",ps12,width = 240,height = 250, units = "mm")
ggsave("figS13.pdf",ps13,width = 240,height = 250, units = "mm")
ggsave("figS14.pdf",ps14,width = 240,height = 250, units = "mm")
ggsave("figS15.pdf",ps15,width = 240,height = 250, units = "mm")
ggsave("figS16.pdf",ps16,width = 240,height = 250, units = "mm")

#### Figure 2 ####

p2a = ggplot(subset(df_skew, s.e == 0 &
                        u %in% c(0.02,0.1,0.5) &
                        skew %in% c(0,0.25,0.5,0.75,1) &
                        r2 == 0.5),
             aes(x = m,
                 y = prop_E_dif,
                 color = skew_fact,
                 group = u_fact:skew_fact,
                 lty = u_fact,
                 size = u_fact)) +
    scale_linetype_manual(values = c("33","solid","12")) +
    scale_size_manual(values = c(1,0.5,1.3), guide = "none") +
    scale_color_manual(values = viridis_pal_p2) +
    geom_line() +
    labs(color = expression(varphi),
         y = expression("Epigenetic divergence ("*E[1] - E[2]*" )"),
         x = "Migration (m)",
         lty = expression(mu),
         title = "Epimutation-migration Balance") +
    theme(legend.position = "none") +
    guides(linetype = guide_legend(override.aes = list(color = "black")))

p2b = ggplot(subset(df_s, m == 0 &
                        u %in% c(0.02,0.1,0.5) &
                        skew %in% c(0,0.25,0.5,0.75,1) &
                        r2 == 0.5),
             aes(x = s.e,
                 y = prop_E_dif,
                 color = skew_fact,
                 group = u_fact:skew_fact,
                 lty = u_fact,
                 size = u_fact)) +
    scale_color_manual(values = viridis_pal_p2) +
    scale_linetype_manual(values = c("12","33","solid")) +
    scale_size_manual(values = c(1.3,1,0.5), guide = "none") +
    geom_line() +
    labs(color = expression(varphi),
         y = expression("Epigenetic divergence ("*E[1] - E[2]*" )"),
         x = "Selection (s)",
         lty = expression(mu),
         title = "Epimutation-selection Balance") +
    guides(linetype = guide_legend(override.aes = list(color = "black")))

p2 = plot_grid(p2a,p2b,
               ncol = 2,
               labels = c("A","B"),
               rel_widths = c(1,1.2),
               align = "h")

ggsave("fig2.pdf",p2,width = 240,height = 120, units = "mm")

#### Figure 3 and Figure S9 ####

Dnorm_cols = grep("Dnorm_",colnames(df_LD))
Dnorm_cols_string = colnames(df_LD)[Dnorm_cols]

# Function to find the year when values reach zero (within error range)
find_zero_gen <- function(df, error = 1e-12) {
    # Extract years from column names
    col_names <- colnames(df)
    gens <- as.numeric(gsub(".*_gen_([0-9]+)$", "\\1", col_names))
    
    # Initialize result vector
    result <- rep(NA, nrow(df))
    
    # For each row in the dataframe
    for (i in 1:nrow(df)) {
        # Check each column in order
        for (j in 1:length(col_names)) {
            val <- df[i, j]
            
            # If value is close to zero (within error range), record the year and break
            if (!is.na(val) && abs(val) <= error) {
                result[i] <- gens[j]
                break
            }
        }
        
        # if(is.na(result[i])){result[i] = 2e6}
    }
    
    return(result)
}

# Add new column with zero year
df_LD$zero_gen <- find_zero_gen(df_LD[,Dnorm_cols_string])

p3 = ggplot(subset(df_LD, s.e == 0.01 &
                        m %in% c(0.001) &
                        skew %in% c(0,1)),
             aes(x = r2,
                 y = log(zero_gen-2e6+1),
                 color = u_fact,
                 group = u_fact:skew_fact,
                 lty = skew_fact,
                 size = skew_fact)) +
    scale_color_manual(values = viridis_pal_p3,
                       breaks = c(0,0.02,0.1,0.2,0.3,0.4,0.5)) +
    geom_line() +
    guides(lty = guide_legend(override.aes = list(linewidth = c(1.3,0.5)))) +
    scale_linetype_manual(values = c("12","solid")) +
    scale_size_manual(values = c(1.3,0.5), guide = "none") +
    labs(x = "Recombination (r)",
         y = "log(Gen. to Linkage Equilibrium)",
         color = expression(mu),
         lty = expression(varphi))

ggsave("fig3.pdf",p3,width = 120,height = 120, units = "mm")

ps9 = ggplot(subset(df_LD, skew %in% c(0,1)),
            aes(x = r2,
                y = log(zero_gen-2e6+1),
                color = u_fact,
                group = u_fact:skew_fact,
                lty = skew_fact,
                size = skew_fact)) +
    scale_color_manual(values = viridis_pal_p3,
                       breaks = c(0,0.02,0.1,0.2,0.3,0.4,0.5)) +
    geom_line() +
    guides(lty = guide_legend(override.aes = list(linewidth = c(1.3,0.5)))) +
    scale_linetype_manual(values = c("12","solid")) +
    scale_size_manual(values = c(1.3,0.5), guide = "none") +
    labs(x = "Recombination (r)",
         y = "log(Gen. to Linkage Equilibrium)",
         color = expression(mu),
         lty = expression(varphi)) +
    facet_wrap(~as.factor(s.e):as.factor(m), labeller = facet_labels_ps9)

ggsave("figS9.pdf",ps9,width = 240,height = 250, units = "mm")

#### Figure 4 ####

df_r$m_fact = factor(df_r$m, levels = c(unique(df_r$m)[5:1]))

p4a <- ggplot(subset(df_r,
                    skew == 0 &
                        m %in% c(0.001,0.1,0.5) &
                        s.e %in% c(0.5)),
             aes(x = r2,
                 y = RI,
                 color = u_fact,
                 group = u_fact,
                 size = u0,
                 lty = u0,
                 order = u0)) +
    scale_color_manual(values = viridis_pal_p3,
                       breaks = c(0,0.02,0.1,0.2,0.3,0.4,0.5)) +
    scale_linetype_manual(values = c("solid","12")) +
    scale_size_manual(values = c(0.5,1.3), guide = "none") +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    geom_line() +
    guides(lty = "none",
           color = guide_legend(
               override.aes = list(
                   linetype = c("12", "solid", "solid", "solid", "solid", "solid", "12"), 
                   linewidth = c(1.3, 0.5, 0.5, 0.5, 0.5, 0.5, 1.3)
               )
           )) +
    labs(x = "Recombination (r)",
         y = "RI", color = expression(mu),
         title = expression("Neutral Epimutation ("*varphi == 0*")")) +
    ylim(c(-0.01,1)) +
    facet_wrap(~m_fact, labeller = facet_labels_p4a)

df_r$m_fact = factor(df_r$m, levels = c(unique(df_r$m)[5:1]))

p4b <- ggplot(subset(df_r,
                    skew == 1 &
                        m %in% c(0.001,0.1,0.5) &
                        s.e %in% c(0.5)),
             aes(x = r2,
                 y = RI,
                 color = u_fact,
                 group = u_fact,
                 lty = u0,
                 size = u0)) +
    scale_color_manual(values = viridis_pal_p3,
                       breaks = c(0,0.02,0.1,0.2,0.3,0.4,0.5)) +
    scale_linetype_manual(values = c("solid","12")) +
    scale_size_manual(values = c(0.5,1.3), guide = "none") +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    geom_line() +
    guides(lty = "none",
           color = guide_legend(
               override.aes = list(
                   linetype = c("12", "solid", "solid", "solid", "solid", "solid", "12"), 
                   linewidth = c(1.3, 0.5, 0.5, 0.5, 0.5, 0.5, 1.3)
               )
           )) +
    labs(x = "Recombination (r)",
         y = "RI", color = expression(mu),
         title = expression("Environmental Induction ("*varphi == 1*")")) +
    ylim(c(-0.01,1)) +
    facet_wrap(~m_fact, labeller = facet_labels_p4b) +
    theme(legend.position = "none")

p4a_leg = get_legend(p4a)
p4a_noleg = p4a + theme(legend.position = "none")

p4_noleg = plot_grid(p4a_noleg,p4b,ncol = 1)
p4 = plot_grid(p4_noleg,p4a_leg, ncol = 2, rel_widths = c(1,0.1))

ggsave("fig4.pdf",p4,width = 240,height = 200, units = "mm")

#### Figure 5 ####

p5_base = psX_base %+% subset(df_u, skew == 1 & r2 %in% c(0.001,0.5) & s.e == 0.05) +
    # scale_color_manual(values = viridis_pal_p3) +
    facet_null() +
    labs(title = "Weak Selection (s = 0.05)") +
    ylim(c(0,1)) +
    facet_wrap(~r2, labeller = facet_labels_p5)

# Create inlay plot for r2 = 0.001 (with y limits 0-0.2)
inlay_001 = psX_base %+% subset(df_u, skew == 1 & r2 == 0.001 & s.e == 0.05) +
    facet_null() +
    theme(title = element_blank(),
          legend.position = "none") +
    coord_cartesian(ylim = c(0,0.10))

# Create inlay plot for r2 = 0.5 (with y limits 0-0.2)
inlay_05 = psX_base %+% subset(df_u, skew == 1 & r2 == 0.5 & s.e == 0.05) +
    facet_null() +
    theme(title = element_blank(),
          legend.position = "none") +
    coord_cartesian(ylim = c(0,0.10))

p5 = ggdraw(p5_base) +
    # Add inlay for r2 = 0.001 facet (left side)
    draw_plot(inlay_001, x = 0.15, y = 0.2, width = 0.33, height = 0.60) +
    # Add inlay for r2 = 0.5 facet (right side)
    draw_plot(inlay_05, x = 0.55, y = 0.2, width = 0.33, height = 0.60)

ggsave("fig5.pdf",p5,width = 240,height = 120, units = "mm")

#### Figure 6 ####

p6a = psX_base %+% subset(df_u_2nd, skew == 1 & r2 == 0.001 & s.e == 0.5) +
    facet_null() +
    labs(title = "Secondary Contact") +
    theme(plot.title = element_text(size = 12)) +
    ylim(c(0,1))

p6b = psX_base %+% subset(df_u, skew == 1 & r2 == 0.001 & s.e == 0.5) +
    facet_null() +
    labs(title = "Primary Contact") +
    theme(plot.title = element_text(size = 12),
          legend.position = "none") +
    ylim(c(0,1))

p6c = ggplot(subset(df_all,s.e == 0.1 & skew == 1 & r2 == 0.001 & m %in% c(0.1)),
              aes(x = allo_num,
                  y = RI,
                  color = u_fact,
                  lty = u0,
                  size = u0,
                  group = u_fact)) +
    scale_color_manual(values = viridis_pal_p3,
                       breaks = c(0,0.02,0.1,0.2,0.3,0.4,0.5)) +
    scale_linetype_manual(values = c("solid","12")) +
    scale_size_manual(values = c(0.5,1.3), guide = "none") +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    geom_line() +
    guides(lty = "none",
           color = guide_legend(
               override.aes = list(
                   linetype = c("12", "solid", "solid", "solid", "solid", "solid", "12"), 
                   linewidth = c(1.3, 0.5, 0.5, 0.5, 0.5, 0.5, 1.3)
               )
           )) +
    xlim(c(0,100)) +
    ylim(c(0,0.35)) +
    # geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    labs(x = "Generations of allopatry",
         color = expression(mu),
         title = expression("m = 0.1; s = 0.1; r = 0.001")) +
    theme(legend.position = "none")

p6d = ggplot(subset(df_all,s.e == 0.5 & skew == 1 & r2 == 0.5 & m %in% c(0.5)),
             aes(x = allo_num,
                 y = RI,
                 color = u_fact,
                 lty = u0,
                 size = u0)) +
    scale_color_manual(values = viridis_pal_p3,
                       breaks = c(0,0.02,0.1,0.2,0.3,0.4,0.5)) +
    scale_linetype_manual(values = c("solid","12")) +
    scale_size_manual(values = c(0.5,1.3), guide = "none") +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    geom_line() +
    guides(lty = "none",
           color = guide_legend(
               override.aes = list(
                   linetype = c("12", "solid", "solid", "solid", "solid", "solid", "12"), 
                   linewidth = c(1.3, 0.5, 0.5, 0.5, 0.5, 0.5, 1.3)
               )
           )) +
    xlim(c(0,25)) +
    ylim(c(0,0.35)) +
    # geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    labs(x = "Generations of allopatry",
         color = expression(mu),
         title = expression("m = 0.5; s = 0.5; r = 0.5")) +
    theme(legend.position = "none")

p6a_leg = get_legend(p6a)
p6a_noleg = p6a + theme(legend.position = "none")

# p6cd = plot_grid(p6c,p6d, ncol = 2, labels = c("c","d"))
# p6bcd = plot_grid(p6b,p6cd, ncol = 1, rel_heights = c(1,0.8), labels = c("b",""))
p6_noleg = plot_grid(p6a_noleg,p6b,p6c,p6d,
               labels = c("AUTO"),
               ncol = 2)
p6 = plot_grid(p6_noleg,p6a_leg, ncol = 2,
               rel_widths = c(1,0.1))

ggsave("fig6.pdf",p6,width = 240,height = 200, units = "mm")

#### Figures S17 to S19 ####

# make new variables which are NA outside the region of interest for plotting
df_all$RI_crop = df_all$RI
df_all$RI_crop[df_all$s.e == 0.5 & df_all$allo_num > 25] = NA
df_all$RI_crop[df_all$s.e == 0.1 & df_all$allo_num > 100] = NA
df_all$allo_crop = df_all$allo_num
df_all$allo_crop[df_all$s.e == 0.5 & df_all$allo_num > 25] = NA
df_all$allo_crop[df_all$s.e == 0.1 & df_all$allo_num > 100] = NA

# Plot with m = 0.01
ps17 = ggplot(subset(df_all, m == 0.01),
            aes(size = u0, x = allo_crop, y = RI_crop, color = u_fact, lty = u0)) +
    scale_color_manual(values = viridis_pal_p3,
                       breaks = c(0,0.02,0.1,0.2,0.3,0.4,0.5)) +
    scale_linetype_manual(values = c("solid","12")) +
    scale_size_manual(values = c(0.5,1.3), guide = "none") +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    geom_line() +
    guides(lty = "none",
           color = guide_legend(
               override.aes = list(
                   linetype = c("12", "solid", "solid", "solid", "solid", "solid", "12"), 
                   linewidth = c(1.3, 0.5, 0.5, 0.5, 0.5, 0.5, 1.3)
               )
           )) +
    # geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    facet_wrap(~as.factor(r2):as.factor(s.e), ncol = 3, scales = "free", labeller = facet_labels) +  # Apply custom labeller
    labs(x = "Generations of allopatry",
         y = "RI",
         color = expression(mu),
         title = expression("RI with a period of allopatry ("*varphi*"= 1, m = 0.01)"))

# Plot with m = 0.1
ps18 = ggplot(subset(df_all, m == 0.1),
            aes(size = u0, x = allo_crop, y = RI_crop, color = u_fact, lty = u0)) +
    scale_color_manual(values = viridis_pal_p3,
                       breaks = c(0,0.02,0.1,0.2,0.3,0.4,0.5)) +
    scale_linetype_manual(values = c("solid","12")) +
    scale_size_manual(values = c(0.5,1.3), guide = "none") +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    geom_line() +
    guides(lty = "none",
           color = guide_legend(
               override.aes = list(
                   linetype = c("12", "solid", "solid", "solid", "solid", "solid", "12"), 
                   linewidth = c(1.3, 0.5, 0.5, 0.5, 0.5, 0.5, 1.3)
               )
           )) +facet_wrap(~as.factor(r2):as.factor(s.e), ncol = 3, scales = "free", labeller = facet_labels) +  # Apply custom labeller
    labs(x = "Generations of allopatry",
         y = "RI",
         color = expression(mu),
         title = expression("RI with a period of allopatry ("*varphi*"= 1, m = 0.1)"))

# Plot with m = 0.5
ps19 = ggplot(subset(df_all, m == 0.5),
            aes(size = u0, x = allo_crop, y = RI_crop, color = u_fact, lty = u0)) +
    scale_color_manual(values = viridis_pal_p3,
                       breaks = c(0,0.02,0.1,0.2,0.3,0.4,0.5)) +
    scale_linetype_manual(values = c("solid","12")) +
    scale_size_manual(values = c(0.5,1.3), guide = "none") +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    geom_line() +
    guides(lty = "none",
           color = guide_legend(
               override.aes = list(
                   linetype = c("12", "solid", "solid", "solid", "solid", "solid", "12"), 
                   linewidth = c(1.3, 0.5, 0.5, 0.5, 0.5, 0.5, 1.3)
               )
           )) +facet_wrap(~as.factor(r2):as.factor(s.e), ncol = 3, scales = "free", labeller = facet_labels) +  # Apply custom labeller
    labs(x = "Generations of allopatry",
         y = "RI",
         color = expression(mu),
         title = expression("RI with a period of allopatry ("*varphi*"= 1, m = 0.5)"))

ggsave("figS17.pdf",ps17,width = 240,height = 250, units = "mm")
ggsave("figS18.pdf",ps18,width = 240,height = 250, units = "mm")
ggsave("figS19.pdf",ps19,width = 240,height = 250, units = "mm")

#### analytical comparison plots ####

#### Analytical predictions from the geometric series approximation ####

compute_geom_RI <- function(m, p, s, u, y, h, r) {
    
    q = 1-p
    i=u*(1+y)
    re = r+(1-r)*i
    
    # Compute each component
    term1 <- p^2
    term2 <- q^2 * (1 - s) * i
    term3 <- q^2 * ((1 - s) * (1 - i) * (1 - s * h) * re) / (1 - (1 - s * h) * (1 - re))
    term4 <- p * q * ((1 - s * h) * re) / (1 - (1 - s * h) * (1 - re)) + p * q * (1 - s * h)
    
    # Compute the final value
    me <- m * (term1 + term2 + term3 + term4)
    RI = 1-me/m
    
    # Return the result
    return(RI)
}

geom_RI = numeric()
for(i in 1:nrow(df_r)) {
    geom_RI[i] = compute_geom_RI(m = df_r$m[i],
                               p = 1-df_r$prop_E_1[i],
                               s = df_r$s.e[i],
                               u = as.numeric(as.character(df_r$u[i])),
                               y = as.numeric(as.character(df_r$skew[i])),
                               h = df_r$h.e[i],
                               r = df_r$r2[i])
}

df_r$geom_RI = geom_RI

#### plotting against simulation results (Figure S2 & S3) ####

dfl = pivot_longer(df_r,cols = c("geom_RI","RI"),names_to = "RI_type",values_to = "RI_all")
dfl$RI_type = factor(dfl$RI_type, levels = c("geom_RI","RI"))

dfl$cat = paste0(dfl$RI_type,dfl$u)
dfl$cat = factor(dfl$cat, levels = unique(dfl$cat)[c(seq(2,14,2),seq(1,13,2))])

ps2 = ggplot(subset(dfl, u %in% c(0,0.02,0.1,0.3,0.5) &
                       skew %in% c(1) &
                        m %in% c(0.000001,0.1,0.5) &
                       RI_type %in% c("RI","geom_RI") & h.e %in% c(0.5)),
            aes(x = r2,
                y = RI_all,
                color = cat,
                group = cat,
                lty = RI_type,
                size = RI_type)) +
    scale_linetype_manual(values = c("12","solid"), labels = c("Analytical","Simulation")) +
    scale_size_manual(values = c(1.3,0.5), labels = c("Analytical","Simulation"), guide = "none") +
    scale_color_manual(values = viridis_pal_ps2, labels = c(0,0.02,0.1,0.3,0.5,0,0.02,0.1,0.3,0.5)) +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    geom_line() +
    guides(lty = guide_legend(override.aes = list(linewidth = c(1.3,0.5)))) +
    labs(color = expression(mu),
         lty = "Method",
         y = "RI",
         x = "Recombination (r)",
         title = expression("Geometric series approximation ("*varphi~"= 1)")) +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.text = element_text(size = 10),
          legend.key.height = unit(3.5,"mm"),
          legend.key.width = unit(5,"mm")) +
    facet_wrap(~as.factor(s.e):as.factor(m), scales = "free_y",ncol = 3,
               labeller = facet_labels_ps2)

ps3 = ggplot(subset(dfl, u %in% c(0,0.02,0.1,0.3,0.5) &
                        skew %in% c(0.5) &
                        m %in% c(0.000001,0.1,0.5) &
                        RI_type %in% c("RI","geom_RI") & h.e %in% c(0.5)),
             aes(x = r2,
                 y = RI_all,
                 color = cat,
                 group = cat,
                 lty = RI_type,
                 size = RI_type)) +
    scale_linetype_manual(values = c("12","solid"), labels = c("Analytical","Simulation")) +
    scale_size_manual(values = c(1.3,0.5), labels = c("Analytical","Simulation"), guide = "none") +
    scale_color_manual(values = viridis_pal_ps2, labels = c(0,0.02,0.1,0.3,0.5,0,0.02,0.1,0.3,0.5)) +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    geom_line() +
    guides(lty = guide_legend(override.aes = list(linewidth = c(1.3,0.5)))) +
    labs(color = expression(mu),
         lty = "Method",
         y = "RI",
         x = "Recombination (r)",
         title = expression("Geometric series approximation ("*varphi~"= 0.5)")) +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.text = element_text(size = 10),
          legend.key.height = unit(3.5,"mm"),
          legend.key.width = unit(5,"mm")) +
    facet_wrap(~as.factor(s.e):as.factor(m), scales = "free_y",ncol = 3,
               labeller = facet_labels_ps2)

ggsave("figS2.pdf",ps2,width = 240,height = 250, units = "mm")
ggsave("figS3.pdf",ps3,width = 240,height = 250, units = "mm")

#### Analytical predictions from no transmission approximation ####

no_trans_RI <- function(m, s, h, y) {
    numerator <-   4-s-2*s*h-2*s*y-s*y^2+2*s*h*y^2
    denominator <- 4-s-2*s*h+2*s*y-s*y^2+2*s*h*y^2-4*m*s*y
    RI_value <- 1 - (numerator / denominator)
    return(RI_value)
}

df_skew$no_trans_RI = no_trans_RI(
    df_skew$m,
    df_skew$s.e,
    df_skew$h.e,
    df_skew$skew
)

#### plotting against simulation results (Figure S4 & S5) ####

dfl2 = pivot_longer(df_skew,cols = c("no_trans_RI","RI"),names_to = "RI_type",values_to = "RI_all")
dfl2$RI_type = factor(dfl2$RI_type, levels = c("no_trans_RI","RI"))

dfl2$cat = paste0(dfl2$RI_type,dfl2$skew)

ps4 = ggplot(subset(dfl2, u == 0.5 &
                          s.e != 0 &
                          skew %in% c(0,0.25,0.5,0.75,1) &
                          RI_type %in% c("RI","no_trans_RI") &
                          h.e == 0.5),
               aes(x = m,
                   y = RI_all,
                   color = cat,
                   group = cat,
                   lty = RI_type,
                   size = RI_type)) +
    scale_linetype_manual(values = c("12","solid"), labels = c("Analytical","Simulation")) +
    scale_size_manual(values = c(1.3,0.5), labels = c("Analytical","Simulation"), guide = "none") +
    scale_color_manual(values = viridis_pal_ps2[c(6:10,1:5)], labels = c(0,0.25,0.5,0.75,1,0,0.25,0.5,0.75,1)) +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    guides(lty = guide_legend(override.aes = list(linewidth = c(1.3,0.5)))) +
    geom_line() +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    labs(color = expression(varphi),
         lty = "Method",
         size = "Method",
         y = "RI",
         x = "Migration (m)",
         title = expression("No Transmission Approximation ("*mu~"= 0.5)")) +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.text = element_text(size = 10),
          legend.key.height = unit(3.5,"mm"),
          legend.key.width = unit(5,"mm")) +
    facet_wrap(~as.factor(s.e):as.factor(r2), scales = "free_y", ncol = 3,
               labeller = facet_labels_ps4)


ps5 = ggplot(subset(dfl2, u == 0.1 &
                          s.e != 0 &
                           skew %in% c(0,0.25,0.5,0.75,1) &
                           RI_type %in% c("RI","no_trans_RI") &
                           h.e == 0.5),
                aes(x = m,
                    y = RI_all,
                    color = cat,
                    group = cat,
                    lty = RI_type,
                    size = RI_type)) +
    scale_linetype_manual(values = c("12","solid"), labels = c("Analytical","Simulation")) +
    scale_color_manual(values = viridis_pal_ps2[c(6:10,1:5)], labels = c(0,0.25,0.5,0.75,1,0,0.25,0.5,0.75,1)) +
    scale_size_manual(values = c(1.3,0.5), labels = c("Analytical","Simulation"), guide = "none") +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    guides(lty = guide_legend(override.aes = list(linewidth = c(1.3,0.5)))) +
    geom_line() +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    labs(color = expression(varphi),
         lty = "Method",
         y = "RI",
         x = "Migration (m)",
         title = expression("No Transmission Approximation ("*mu~"= 0.1)")) +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.text = element_text(size = 10),
          legend.key.height = unit(3.5,"mm"),
          legend.key.width = unit(5,"mm")) +
    facet_wrap(~as.factor(s.e):as.factor(r2), scales = "free_y", ncol = 3,
               labeller = facet_labels_ps4)

ggsave("figS4.pdf",ps4,width = 240,height = 200, units = "mm")
ggsave("figS5.pdf",ps5,width = 240,height = 200, units = "mm")

#### Analytical predictions from asexual approximation ####

# see Kickstarter_revisions_analytical_solving.ipynb for generating the USM/MUS matricies

# Core function to compute equilibrium state vector using Matrix package
compute_equilibrium <- function(M, X0, n = NULL, use_fallbacks = TRUE) {
    # If n is not specified, use the matrix dimension
    if(is.null(n)) n <- nrow(M)
    
    tryCatch({
        # Convert to dense matrix format for eigenvalue calculation
        M_dense <- as(M, "dgeMatrix")
        
        # Compute eigendecomposition using Matrix's methods
        M_eig <- eigen(M_dense)
        
        # Create lambda matrix with better numerical precision
        lambda_values <- M_eig$values
        tol <- max(dim(M_dense)) * max(abs(lambda_values)) * .Machine$double.eps
        
        # For dominant eigenvalue approach, use this block:
        max_value <- max(Re(lambda_values))
        lambda_matrix <- Matrix(0, nrow = nrow(M), ncol = ncol(M))
        diag_positions <- which(sapply(Re(lambda_values), function(x) 
            isTRUE(all.equal(x, max_value))))
        for (i in diag_positions) {
            lambda_matrix[i, i] <- 1
        }
        
        # Alternative approach using power:
        # lambda_values[abs(lambda_values) < tol] <- 0
        # lambda_matrix <- Diagonal(x = lambda_values^n)
        
        # Use pseudoinverse with better tolerance control
        M_vect <- Matrix(M_eig$vectors, sparse = FALSE)
        
        # Calculate pseudoinverse using SVD with more precise tolerance
        svd_result <- svd(M_vect)
        d <- svd_result$d
        tol <- max(dim(M_vect)) * max(d) * .Machine$double.eps
        d_inv <- ifelse(d > tol, 1/d, 0)
        M_vect_inv <- svd_result$v %*% Diagonal(x = d_inv) %*% t(svd_result$u)
        
        # Compute result
        X_result <- as.numeric(M_vect %*% lambda_matrix %*% M_vect_inv %*% X0)
        X_result <- X_result / sum(X_result)
        
        return(X_result)
    }, error = function(e) {
        warning(paste("Error in primary matrix computation method:", e$message))
        
        # If fallbacks are disabled, return NA immediately
        if (!use_fallbacks) {
            warning("Fallback methods disabled, returning NA")
            return(rep(NA, length(X0)))
        }
        
        # Otherwise try alternative approaches
        warning("Attempting fallback methods...")
        
        # Alternative approach using Matrix power operations
        try({
            warning("Trying Matrix power method...")
            # Use direct matrix power
            M_pow <- M
            for (i in 2:20) {  # Using 20 iterations
                M_pow <- M_pow %*% M
            }
            X_result <- as.numeric(M_pow %*% X0)
            X_result <- X_result / sum(X_result)
            warning("Matrix power method succeeded")
            return(X_result)
        }, silent = TRUE)
        
        # Second fallback using iterative approach
        try({
            warning("Trying iterative convergence method...")
            X_curr <- as.numeric(X0)
            for (i in 1:100) {
                X_next <- as.numeric(M %*% X_curr)
                X_next <- X_next / sum(X_next)
                
                # Check for convergence
                if (sum(abs(X_next - X_curr)) < 1e-10) {
                    warning("Iterative method converged successfully")
                    return(X_next)
                }
                X_curr <- X_next
            }
            warning("Iterative method reached maximum iterations")
            return(X_curr)
        }, silent = TRUE)
        
        warning("All methods failed, returning NA")
        return(rep(NA, length(X0)))
    })
}

# Prepare matrices for computation (convert to Matrix format and regularize)
prepare_matrices <- function(MUS, USM, Mc, epsilon = 1e-10) {
    # Convert to Matrix package format
    MUS <- Matrix(MUS, sparse = FALSE)
    USM <- Matrix(USM, sparse = FALSE)
    Mc <- Matrix(Mc, sparse = FALSE)
    
    # Apply regularization for numerical stability
    diag(MUS) <- diag(MUS) + epsilon
    diag(USM) <- diag(USM) + epsilon
    
    # Compute and check condition numbers
    cond_MUS <- rcond(MUS)
    cond_USM <- rcond(USM)
    
    # Print warnings if matrices are too ill-conditioned
    if (cond_MUS < 1e-14) {
        warning(sprintf("MUS matrix is ill-conditioned: reciprocal condition number = %g", cond_MUS))
    }
    if (cond_USM < 1e-14) {
        warning(sprintf("USM matrix is ill-conditioned: reciprocal condition number = %g", cond_USM))
    }
    
    return(list(MUS = MUS, USM = USM, Mc = Mc))
}

# General function to compute RI given matrices and indices
compute_RI <- function(MUS, USM, Mc, X0, A0_indices, Aeq_indices = NULL, use_fallbacks = TRUE) {
    # If Aeq_indices not provided, use same as A0_indices
    if(is.null(Aeq_indices)) Aeq_indices <- A0_indices
    
    # Prepare matrices
    matrices <- prepare_matrices(MUS, USM, Mc)
    
    # Primary state after first matrix operation
    Xprimary <- compute_equilibrium(matrices$USM, X0, use_fallbacks = use_fallbacks)
    
    # Check if computation succeeded
    if(any(is.na(Xprimary))) {
        warning("Failed to compute primary equilibrium")
        return(NA)
    }
    
    # Apply contact matrix
    Xcontact <- as.numeric(matrices$Mc %*% Xprimary)
    
    # Final equilibrium
    Xequil <- compute_equilibrium(matrices$MUS, Xcontact, use_fallbacks = use_fallbacks)
    
    # Check if computation succeeded
    if(any(is.na(Xequil))) {
        warning("Failed to compute final equilibrium")
        return(NA)
    }
    
    # Calculate RI based on indices
    A0 <- sum(Xcontact[A0_indices$full_indices] + 
                  ifelse(is.null(A0_indices$half_indices), 0, 0.5 * Xcontact[A0_indices$half_indices]))
    
    Aeq <- sum(Xequil[Aeq_indices$full_indices] + 
                   ifelse(is.null(Aeq_indices$half_indices), 0, 0.5 * Xequil[Aeq_indices$half_indices]))
    
    RI <- 1 - Aeq / A0
    
    return(RI)
}

# Specific RI calculation functions using the generalized approach
calc_asex_RI <- function(m, mu, phi, s, use_fallbacks = TRUE) {
    # Matrix definitions for haploid model
    MUS <- matrix(c(
        (1 - m) * (- mu * (1 - phi) + 1), mu * (1 - m) * (1 - s) * (phi + 1), m * (1 - s) * (- mu * (phi + 1) + 1), m * mu * (1 - phi), 0, 0, 0, 0,
        mu * (1 - m) * (1 - phi), (1 - m) * (1 - s) * (- mu * (phi + 1) + 1), m * mu * (1 - s) * (phi + 1), m * (- mu * (1 - phi) + 1), 0, 0, 0, 0,
        m * (- mu * (1 - phi) + 1), m * mu * (1 - s) * (phi + 1), (1 - m) * (1 - s) * (- mu * (phi + 1) + 1), mu * (1 - m) * (1 - phi), 0, 0, 0, 0,
        m * mu * (1 - phi), m * (1 - s) * (- mu * (phi + 1) + 1), mu * (1 - m) * (1 - s) * (phi + 1), (1 - m) * (- mu * (1 - phi) + 1), 0, 0, 0, 0,
        0, 0, 0, 0, (1 - m) * (- mu * (1 - phi) + 1), mu * (1 - m) * (1 - s) * (phi + 1), m * (1 - s) * (- mu * (phi + 1) + 1), m * mu * (1 - phi),
        0, 0, 0, 0, mu * (1 - m) * (1 - phi), (1 - m) * (1 - s) * (- mu * (phi + 1) + 1), m * mu * (1 - s) * (phi + 1), m * (- mu * (1 - phi) + 1),
        0, 0, 0, 0, m * (- mu * (1 - phi) + 1), m * mu * (1 - s) * (phi + 1), (1 - m) * (1 - s) * (- mu * (phi + 1) + 1), mu * (1 - m) * (1 - phi),
        0, 0, 0, 0, m * mu * (1 - phi), m * (1 - s) * (- mu * (phi + 1) + 1), mu * (1 - m) * (1 - s) * (phi + 1), (1 - m) * (- mu * (1 - phi) + 1)
    ), nrow=8, byrow=TRUE)
    
    USM <- matrix(c(
        (1 - m) * (-mu * (1 - phi) + 1), mu * (1 - m) * (1 - s) * (phi + 1), m * (-mu * (1 - phi) + 1), m * mu * (1 - s) * (phi + 1), 0, 0, 0, 0,
        mu * (1 - m) * (1 - phi), (1 - m) * (1 - s) * (-mu * (phi + 1) + 1), m * mu * (1 - phi), m * (1 - s) * (-mu * (phi + 1) + 1), 0, 0, 0, 0,
        m * (1 - s) * (-mu * (phi + 1) + 1), m * mu * (1 - phi), (1 - m) * (1 - s) * (-mu * (phi + 1) + 1), mu * (1 - m) * (1 - phi), 0, 0, 0, 0,
        m * mu * (1 - s) * (phi + 1), m * (-mu * (1 - phi) + 1), mu * (1 - m) * (1 - s) * (phi + 1), (1 - m) * (-mu * (1 - phi) + 1), 0, 0, 0, 0,
        0, 0, 0, 0, (1 - m) * (-mu * (1 - phi) + 1), mu * (1 - m) * (1 - s) * (phi + 1), m * (-mu * (1 - phi) + 1), m * mu * (1 - s) * (phi + 1),
        0, 0, 0, 0, mu * (1 - m) * (1 - phi), (1 - m) * (1 - s) * (-mu * (phi + 1) + 1), m * mu * (1 - phi), m * (1 - s) * (-mu * (phi + 1) + 1),
        0, 0, 0, 0, m * (1 - s) * (-mu * (phi + 1) + 1), m * mu * (1 - phi), (1 - m) * (1 - s) * (-mu * (phi + 1) + 1), mu * (1 - m) * (1 - phi),
        0, 0, 0, 0, m * mu * (1 - s) * (phi + 1), m * (-mu * (1 - phi) + 1), mu * (1 - m) * (1 - s) * (phi + 1), (1 - m) * (-mu * (1 - phi) + 1)
    ), nrow = 8, byrow = TRUE)
    
    Mc <- matrix(c(
        1 - m, 0,    m,    0, 0, 0, 0, 0,
        0,    1 - m, 0,    m, 0, 0, 0, 0,
        0,    0,    1 - m, 0, 0, 0, 0, 0,
        0,    0,    0,    1 - m, 0, 0, 0, 0,
        0,    0,    0,    0, 0, 0, 0, 0,
        0,    0,    0,    0, 0, 0, 0, 0,
        m,    0,    0,    0, 0, 0, 0, 0,
        0,    m,    0,    0, 0, 0, 0, 0
    ), nrow = 8, byrow = TRUE)
    
    # Initial state
    X0 = c(0.5, 0, 0, 0.5, 0, 0, 0, 0)
    
    # Define indices for calculating RI
    A0_indices = list(full_indices = 5:8, half_indices = NULL)
    
    # Compute RI
    return(compute_RI(MUS, USM, Mc, X0, A0_indices, use_fallbacks = use_fallbacks))
}

RI_asex_s <- numeric()
RI_asex_sh <- numeric()

for (i in 1:nrow(df_u)) {
    print(i)
    
    # for matching fitness of heterozygotes
    # s = 0.5*df_u$s.e[i]
    # for matching fitness of homozygotes
    # s = df_u$s.e[i]
    # for matching multiplicative fitness
    # s = 1-(1-df_u$s.e[i])^2
    
    RI_asex_s[i] <- calc_asex_RI(m = df_u$m[i],
                               mu = as.numeric(as.character(df_u$u[i])),
                               phi = as.numeric(as.character(df_u$skew[i])),
                               s = df_u$s.e[i])
    RI_asex_sh[i] <- calc_asex_RI(m = df_u$m[i],
                               mu = as.numeric(as.character(df_u$u[i])),
                               phi = as.numeric(as.character(df_u$skew[i])),
                               s = 0.5*df_u$s.e[i])
}

df_u$RI_asex_s <- RI_asex_s
df_u$RI_asex_sh <- RI_asex_sh

#### Plotting against simulation results (Figures S6-S8) ####

dfl3 = pivot_longer(df_u,cols = c("RI","RI_asex_s","RI_asex_sh"),names_to = "RI_type",values_to = "RI_all")
dfl3$RI_type = factor(dfl3$RI_type, levels = c("RI_asex_s","RI_asex_sh","RI"))

dfl3$cat = paste0(dfl3$RI_type,dfl3$u)
dfl3$cat = factor(dfl3$cat, levels = unique(dfl3$cat)[c(seq(1,19,3),seq(2,20,3),seq(3,21,3))])

ps6 = ggplot(subset(dfl3, u %in% c(0,0.02,0.1,0.3,0.5) &
                        r2 %in% c(0) &
                        s.e %in% c(0.01,0.1,0.5) &
                        skew %in% c(0,0.5,1) &
                        RI_type %in% c("RI","RI_asex_sh") &
                        h.e == 0.5),
             aes(x = m,
                 y = RI_all,
                 color = cat,
                 group = cat,
                 lty = RI_type,
                 size = RI_type)) +
    scale_linetype_manual(values = c("12","solid"), labels = c("Analytical","Simulation")) +
    scale_size_manual(values = c(1.3,0.5), labels = c("Analytical","Simulation"), guide = "none") +
    scale_color_manual(values = viridis_pal_ps2, labels = c(0,0.02,0.1,0.3,0.5,0,0.02,0.1,0.3,0.5)) +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    guides(lty = guide_legend(override.aes = list(linewidth = c(1.3,0.5)))) +
    geom_line() +
    labs(color = expression(mu),
         lty = "Method",
         y = "RI",
         x = "Migration (m)",
         title = expression("Asexual haploid approximation matching heterozygote fitness, ("*r~"= 0)")) +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.text = element_text(size = 10),
          legend.key.height = unit(3.5,"mm"),
          legend.key.width = unit(5,"mm")) +
    facet_wrap(~as.factor(s.e):as.factor(skew), scales = "free_y", ncol = 3,
               labeller = facet_labels_ps6)

ps7 = ggplot(subset(dfl3, u %in% c(0,0.02,0.1,0.3,0.5) &
                        r2 %in% c(0) &
                        s.e %in% c(0.01,0.1,0.5) &
                        skew %in% c(0,0.5,1) &
                        RI_type %in% c("RI","RI_asex_s") &
                        h.e == 0.5),
             aes(x = m,
                 y = RI_all,
                 color = cat,
                 group = cat,
                 lty = RI_type,
                 size = RI_type)) +
    scale_linetype_manual(values = c("12","solid"), labels = c("Analytical","Simulation")) +
    scale_size_manual(values = c(1.3,0.5), labels = c("Analytical","Simulation"), guide = "none") +
    scale_color_manual(values = viridis_pal_ps2, labels = c(0,0.02,0.1,0.3,0.5,0,0.02,0.1,0.3,0.5)) +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    guides(lty = guide_legend(override.aes = list(linewidth = c(1.3,0.5)))) +
    geom_line() +
    labs(color = expression(mu),
         lty = "Method",
         y = "RI",
         x = "Migration (m)",
         title = expression("Asexual haploid approximation matching homozygote fitness ("*r~"= 0)")) +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.text = element_text(size = 10),
          legend.key.height = unit(3.5,"mm"),
          legend.key.width = unit(5,"mm")) +
    facet_wrap(~as.factor(s.e):as.factor(skew), scales = "free_y", ncol = 3,
               labeller = facet_labels_ps6)


ps8 = ggplot(subset(dfl3, u %in% c(0,0.02,0.1,0.3,0.5) &
                        r2 %in% c(0.1) &
                        s.e %in% c(0.01,0.1,0.5) &
                        skew %in% c(0,0.5,1) &
                        RI_type %in% c("RI","RI_asex_sh") &
                        h.e == 0.5),
             aes(x = m,
                 y = RI_all,
                 color = cat,
                 group = cat,
                 lty = RI_type,
                 size = RI_type)) +
    scale_linetype_manual(values = c("12","solid"), labels = c("Analytical","Simulation")) +
    scale_size_manual(values = c(1.3,0.5), labels = c("Analytical","Simulation"), guide = "none") +
    scale_color_manual(values = viridis_pal_ps2, labels = c(0,0.02,0.1,0.3,0.5,0,0.02,0.1,0.3,0.5)) +
    geom_hline(yintercept = 0, lty = 3, color = "darkgrey") +
    guides(lty = guide_legend(override.aes = list(linewidth = c(1.3,0.5)))) +
    geom_line() +
    labs(color = expression(mu),
         lty = "Method",
         y = "RI",
         x = "Migration (m)",
         title = expression("Asexual haploid approximation matching heterozygote fitness ("*r~"= 0.1)")) +
    guides(color = guide_legend(override.aes = list(size = 1.5))) +
    theme(legend.text = element_text(size = 10),
          legend.key.height = unit(3.5,"mm"),
          legend.key.width = unit(5,"mm")) +
    facet_wrap(~as.factor(s.e):as.factor(skew), scales = "free_y", ncol = 3,
               labeller = facet_labels_ps6)

## adding inlays

# Get unique values to identify the first row facets
unique_s_e <- c(0.01,0.1,0.5)
# unique_s_e <- c(0.0050,0.0513,0.2930)
unique_skew <- c(0,0.5,1)
first_row_s_e <- unique_s_e[1]  # Get the s.e value for the first row

# Convert main plot to a drawing canvas
ps6_inlay <- ggdraw(ps6)
ps7_inlay <- ggdraw(ps7)
ps8_inlay <- ggdraw(ps8)

# For each facet in the first row, create and add an inlay
for (i in seq_along(unique_skew)) {
    
    sk <- unique_skew[i]
    
    # Create the inlay plot with x-axis limited to 0-0.02
    inlay6 <- ggplot(subset(dfl3, u %in% c(0,0.02,0.1,0.3,0.5) &
                                r2 %in% c(0) &
                                skew == sk &
                                s.e == first_row_s_e &
                                RI_type %in% c("RI","RI_asex_sh") &
                                h.e == 0.5),
                     aes(x = m,
                         y = RI_all,
                         color = cat,
                         group = cat,
                         lty = RI_type,
                         size = RI_type)) +
        scale_linetype_manual(values = c("12","solid")) +
        scale_size_manual(values = c(1.3,0.5)) +
        scale_color_manual(values = viridis_pal_ps2) +
        geom_line() +
        xlim(0, 0.02) +  # Limit y-axis to 0-0.05 as requested
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.text = element_text(size = 6),
              plot.background = element_rect(fill = "white", color = NA, size = 0.5),
              panel.grid = element_blank(),
              plot.margin = margin(0, 0, 0, 0))
    
    # Create the inlay plot with x-axis limited to 0-0.02
    inlay7 <- ggplot(subset(dfl3, u %in% c(0,0.02,0.1,0.3,0.5) &
                                r2 %in% c(0) &
                                skew == sk &
                                s.e == first_row_s_e &
                                RI_type %in% c("RI","RI_asex_s") &
                                h.e == 0.5),
                     aes(x = m,
                         y = RI_all,
                         color = cat,
                         group = cat,
                         lty = RI_type,
                         size = RI_type)) +
        scale_linetype_manual(values = c("12","solid")) +
        scale_size_manual(values = c(1.3,0.5)) +
        scale_color_manual(values = viridis_pal_ps2) +
        geom_line() +
        xlim(0, 0.02) +  # Limit y-axis to 0-0.05 as requested
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.text = element_text(size = 6),
              plot.background = element_rect(fill = "white", color = NA, size = 0.5),
              panel.grid = element_blank(),
              plot.margin = margin(0, 0, 0, 0))
    
    # Create the inlay plot with x-axis limited to 0-0.02
    inlay8 <- ggplot(subset(dfl3, u %in% c(0,0.02,0.1,0.3,0.5) &
                               r2 %in% c(0.1) &
                               skew == sk &
                               s.e == first_row_s_e &
                               RI_type %in% c("RI","RI_asex_sh") &
                               h.e == 0.5),
                    aes(x = m,
                        y = RI_all,
                        color = cat,
                        group = cat,
                        lty = RI_type,
                        size = RI_type)) +
        scale_linetype_manual(values = c("12","solid")) +
        scale_size_manual(values = c(1.3,0.5)) +
        scale_color_manual(values = viridis_pal_ps2) +
        geom_line() +
        xlim(0, 0.02) +  # Limit y-axis to 0-0.05 as requested
        theme(legend.position = "none",
              axis.title = element_blank(),
              axis.text = element_text(size = 6),
              plot.background = element_rect(fill = "white", color = NA, size = 0.5),
              panel.grid = element_blank(),
              plot.margin = margin(0, 0, 0, 0))
    
    # Position for the inlay (adjust as needed)
    # Calculate positions based on the facet position in the first row
    x_pos_vect = c(0.12,0.4,0.68)
    x_pos <- x_pos_vect[i]  # Adjust horizontal placement
    y_pos <- 0.7                             # Near the top of the plot
    width <- 0.18                             # Width of inlay
    height <- 0.18                            # Height of inlay
    
    # Add the inlay to the main plot
    ps6_inlay <- ps6_inlay + 
        draw_plot(inlay6, x = x_pos, y = y_pos, width = width, height = height)
    ps7_inlay <- ps7_inlay + 
        draw_plot(inlay7, x = x_pos, y = y_pos, width = width, height = height)
    ps8_inlay <- ps8_inlay + 
        draw_plot(inlay8, x = x_pos, y = y_pos, width = width, height = height)
    
}

ggsave("figS6.pdf",ps6_inlay,width = 240,height = 200, units = "mm")
ggsave("figS7.pdf",ps7_inlay,width = 240,height = 200, units = "mm")
ggsave("figS8.pdf",ps8_inlay,width = 240,height = 200, units = "mm")

