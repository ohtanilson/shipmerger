rm(list = ls())
# Overview ----

# This file generates the following figures and tables after cleaning raw data.
# 
# - Figures 2, 3, 4, 5, and 9. (Figures 6, 7, and 8 are generated in ship_merger_score_estimation.jl.)
# - Tables 2, 3, 4, and 5.
# 
# All raw data must be located in input folder. 
# All outputs will be located in figuretable folder or output folder. 
# 
# Note that Table 1 does not contain the information about data and

#
# Raw data includes firm names in Japanese character. 
# Some of acquired small firms cannot be translated into English because I could not find how the names were pronounced correctly even in Japanese.
# Thus, id column in output/data_for_maximum_rank_estimation.csv corresponds with firm identity.
# For additional replications, please refer to output/data_for_maximum_rank_estimation.csv like my main estimation part in Julia.

library(tidyr)
library(dplyr)
library(ggplot2)
library(readr)
library(gridExtra)
library(stargazer)
library(entropy)
library(cowplot)
library(reshape2)
library(networkD3) # for sankey diagram
library(dplyr) # for sankey diagram

# Figures ----

## (Figure 1) The trend of the world's shipping gross tonnage (1000 tons). ----

data_totaltons <- 
  readr::read_csv("input/data_totaltons.csv")
d <- 
  data_totaltons %>%
  dplyr::select(-year_j,
                -Number, 
                -shipping_quantity_world,
                -shipping_quantity_japan) %>%
  dplyr::rename(Japan = grosstons_japan,
                UK = grosstons_UK,
                US = grosstons_US,
                Norway = gtosstons_norway,
                Soviet = grosstons_soviet,
                Greek = grosstons_greek,
                Liberia = grosstons_liberia) %>% 
  tidyr::gather(key = "variable", value = "value", -year) %>% 
  # Liberia is the major country regarding the flag of convenience ships, 
  # but it does not reflect economic activity in the shipping industry.
  dplyr::filter(variable != "Liberia")

x <- 
  ggplot(d, aes(x = year, y = value)) + 
  geom_line(aes(color = variable), alpha = 0.6) +
  geom_point(aes(shape = variable, color = variable)) +
  theme_minimal() + 
  theme_classic() +
  ggtitle("World's shipping tonnage") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("year") + 
  ylab("shipping tonnage (1000 tons)") +
  expand_limits(y=0) + 
  geom_vline(xintercept = 1964, linetype = "longdash") + 
  labs(shape = "Country", color = "Country")

figure_name <- "figuretable/shippingtonnage.png"
ggplot2::ggsave(filename = figure_name,
       plot = x,
       device = "png",
       width = 10,
       height = 7)

# Load main data ----

data <- 
  readr::read_csv("input/210627firmdata_processed.csv")
data <- 
  data[,c(-3,-4,-6:-11,-16,-18)] # pick up D/W
total_ton <- 
  data[,4]+data[,5]+data[,6]+data[,7]
colnames(data)[3:7] = 
  c("type",
    "liner",
    "special",
    "tramper",
    "tanker")
carrier_type_tonnage_table <- 
  data %>% 
  dplyr::select(type, 
                liner,
                special, 
                tramper, 
                tanker)
data_hist <- 
  reshape2::melt(carrier_type_tonnage_table,
                 "type") %>% 
  # drop very small firms for illustration
  dplyr::filter(value >= 0.1) 
#colnames(data)=c("type","liner","special","tramper","tanker")
x <- 
  ggplot(data_hist, aes(x = value, fill = variable)) + 
  geom_histogram(position="dodge", alpha=0.8) +
  theme_classic() +
  #ggtitle("histogram of tonnage sizes") +
  #theme(plot.title = element_text(hjust = 0.5)) +
  xlab("tonnage size") + 
  ylab("count") + 
  labs("carrier type") +
  xlim(0,650000) +
  labs(fill="carrier type")
g1 <- x
figure_name <- 
  "figuretable/carrier_type_tonnage.png"
ggsave(filename = figure_name,
       plot = x,
       device = "png",
       width = 10,
       height = 7)

# Type-based histogram ----

total_tonnage_table <-
  data %>% 
  dplyr::select(type,
                liner,
                special,
                tramper,
                tanker) %>% 
  dplyr::filter(is.na(type)!= 1) %>% 
  dplyr::mutate(total = 
                  liner + special + tramper + tanker) %>% 
  dplyr::mutate(firm_type = 
                  recode(type, 
                         "合併会社" = "main",
                         "系列会社" = "affiliate",
                         "専属会社" = "wholly controlled",
                         "グループ外" = "unmatched"
                  ))
x <- 
  ggplot(total_tonnage_table,
         aes(x = total, fill = firm_type)) + 
  geom_histogram(position="dodge", alpha=0.8) +
  theme_classic() +
  #ggtitle("histogram of tonnage sizes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("total tonnage size") +
  ylab("count") + 
  xlim(0,650000) +
  labs(col="firm type") +
  labs(fill = "firm type")
g2 <- x
#ggplot2::ggsave("total_tonnage.png")
figure_name <-
  "figuretable/total_tonnage.png"
ggsave(filename = figure_name,
       plot = x,
       device = "png",
       width = 10,
       height = 7)

## construct total tonnage histogram ----
total_tonnage_table <-
  data %>% 
  dplyr::select(type,
                liner,
                special,
                tramper,
                tanker) %>% 
  dplyr::filter(is.na(type)!= 1) %>% 
  dplyr::mutate(total = 
                  liner + special + tramper + tanker) %>% 
  dplyr::mutate(firm_type = 
                  recode(type, 
                         "合併会社" = "main",
                         "系列会社" = "affiliate",
                         "専属会社" = "wholly controlled",
                         "グループ外" = "unmatched"
                  ))
x <- 
  ggplot(total_tonnage_table,
         aes(x = total, fill = firm_type)) + 
  geom_histogram(position="dodge", alpha=0.8) +
  theme_classic() +
  #ggtitle("histogram of tonnage sizes") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("total tonnage size") +
  ylab("count") + 
  xlim(0,650000) +
  labs(col="firm type") +
  labs(fill = "firm type")
g2 <- x
#ggplot2::ggsave("total_tonnage.png")
figure_name <-
  "figuretable/total_tonnage.png"
ggsave(filename = figure_name,
       plot = x,
       device = "png",
       width = 10,
       height = 7)

# translate bank names
bank_names = 
  c("Fuji","Kogin","Chogin","Yasuda","Daiichi","Mitsui",
    "MitsuiShintaku","Sumitomo","SumitomoShintaku","Mitsubishi",
    "Kangin","MitsubishiShintaku","Kobe","Daiwa","Sanwa","Tokai",
    "YasudaShintaku","others")
colnames(data) = 
  c("id", "name", "type",
    "liner","special","tramper","tanker",
    "group",bank_names)

data_english <- 
  data %>% 
  dplyr::filter(is.na(type)!= 1) %>% 
  dplyr::mutate(special = ifelse(is.na(special) == 1,
                                 0, special)) %>%
  dplyr::mutate(tramper = ifelse(is.na(tramper) == 1,
                                 0, tramper)) %>%
  dplyr::mutate(total = liner + special + tramper + tanker) %>% 
  dplyr::mutate(firm_type = recode(type, 
                                   "合併会社" = "main",
                                   "系列会社" = "affiliate",
                                   "専属会社" = "wholly controlled",
                                   "グループ外" = "unmatched"
  )) %>% 
  dplyr::filter(group != "nodata")
data_english$group <- 
  factor(data_english$group, 
         levels = c("Nippon Yusen",
                    "Mitsui O.S.K. Lines",
                    "Japan Lines", 
                    "Kawasaki Kisen Kaisha",
                    "Yamashita-Shinnihon Kisen", 
                    "Showa Line",
                    "unmatched"))





## (Figure 2) Configuration of tonnage size for each carrier type ----

data_english_extracted <- 
  data_english %>% 
  dplyr::select(liner,
                tramper,
                special,
                tanker,
                firm_type,
                group)
data_english_extracted <- 
  reshape2::melt(data_english_extracted)
temp <- 
  data_english_extracted %>% 
  dplyr::group_by(group,
                  variable) %>% 
  dplyr::summarise(value = sum(value)) %>% 
  dplyr::group_by(group) %>% 
  dplyr::mutate(group_sum = sum(value)) %>% 
  dplyr::ungroup()
x <-
  ggplot(temp,
         aes(x = "", y = value/10000, fill = variable)) + 
  geom_bar(width = 1, stat = "identity", colour = "black") +
  #ggtitle("Histogram of tonnage size for each merger group") +
  theme(plot.title = element_text(hjust = 0.5)) +
  #theme_classic() +
  ylab("total tonnage size (10,000 (D/W) tons)") +
  xlab("") +
  labs(fill = "carrier type") +
  geom_hline(yintercept = 100, linetype = "longdash") + 
  facet_wrap(vars(group), ncol = 4) + 
  theme(legend.position = c(0.9, 0.2)) +
  geom_text(aes(x = "", y = value/10000,
                label = scales::percent(value/group_sum)),
            #position=position_fill(vjust=0.5),
            position=position_stack(0.5),
            size = unit(3.0,"cm"))
figure_name <- 
  "figuretable/carrier_composition_eachgroup.png"
ggsave(filename = figure_name,
       plot = x,
       device = "png",
       width = 10,
       height = 7)

## (Figure 3) Distribution of tonnage size for each firm type. ----

data_english_matched <- 
  data_english %>% 
  dplyr::filter(group != "unmatched")
x <- 
  ggplot(data_english_matched, 
         aes(x = total/10000, fill = firm_type)) + 
  geom_histogram(position="dodge", 
                 alpha=0.8, binwidth = 6, colour = "black") +
  #theme_classic() +
  #ggtitle("Histogram of tonnage size for each merger group") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("total tonnage size (10,000 (D/W) tons)") +
  ylab("count (firm level)") +
  xlim(0,105) +
  ylim(0,7) + 
  labs(fill = "firm type") +
  facet_wrap(vars(group), ncol = 4) +
  theme(legend.position = c(0.9, 0.2))
g3 <- x 
#ggplot2::ggsave("type_dist_eachgroup.png")
figure_name <-
  "figuretable/type_dist_eachgroup.png"
ggsave(filename = figure_name,
       plot = x,
       device = "png",
       width = 10,
       height = 7)

## (Figure 4) Distribution of tonnage size for each carrier type. ----

data_english_extracted <- 
  data_english %>% 
  dplyr::select(liner,
                tramper,
                special,
                tanker,
                firm_type,
                group)
data_english_extracted <- 
  reshape2::melt(data_english_extracted)
x <- 
  ggplot(data_english_extracted, 
         aes(x = value/10000, fill = variable)) + 
  geom_histogram(position="dodge", alpha=0.8,
                 binwidth = 8, colour = "black") +
  #theme_classic() +
  #ggtitle("Histogram of tonnage size for each merger group") +
  theme(plot.title = element_text(hjust = 0.5)) +
  xlab("total tonnage size (10,000 (D/W) tons)") +
  ylab("count (firm-type level)") +
  xlim(0,70) +
  ylim(0,7) + 
  labs(fill = "carrier type") +
  facet_wrap(vars(group), ncol = 4) +
  theme(legend.position = c(0.9, 0.2))
figure_name <- 
  "figuretable/carrier_dist_eachgroup.png"
ggsave(filename = figure_name,
       plot = x,
       device = "png",
       width = 10,
       height = 7)


## (Figure 5) Pie charts ----

g5_liner <- 
  ggplot(data_english,
         aes(x ="", y = liner/10000, fill = group)) + 
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ggtitle("Liner") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("share of tonnage size") +
  xlab("") +
  coord_polar("y", start=0) +
  theme(legend.position = "none")
g5_tramper <- 
  ggplot(data_english,
         aes(x ="", y = tramper/10000, fill = group)) + 
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ggtitle("Tramper") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("share of tonnage size") +
  xlab("") +
  coord_polar("y", start=0) +
  theme(legend.position = "none")
g5_special <-
  ggplot(data_english, 
         aes(x ="", y = special/10000, fill = group)) + 
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ggtitle("Special") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("share of tonnage size") +
  xlab("") +
  coord_polar("y", start=0) +
  theme(legend.position = "none")
g5_tanker <-
  ggplot(data_english, 
         aes(x ="", y = tanker/10000, fill = group)) + 
  geom_bar(width = 1, stat = "identity") +
  theme_classic() +
  ggtitle("Tanker") +
  theme(plot.title = element_text(hjust = 0.5)) +
  ylab("share of tonnage size") + 
  xlab("") +
  coord_polar("y", start=0)
legend <- 
  cowplot::get_legend(g5_tanker)
g5_tanker <-
  g5_tanker +
  theme(legend.position = "none")

x <- 
  gridExtra::grid.arrange(g5_liner, 
                          g5_tramper,
                          legend,
                          g5_special,
                          g5_tanker,
                          ncol = 3,
                          nrow = 2)
#ggplot2::ggsave("carrier_share_eachgroup.png",temp_g)
figure_name <-
  "figuretable/carrier_share_eachgroup.png"
ggsave(filename = figure_name,
       plot = x,
       device = "png",
       width = 10,
       height = 7)

# Descriptive summary ----
## construct share variables ----
data_english <- 
  data_english %>% 
  dplyr::mutate(HHI = (liner/total)^2 + (special/total)^2 +
                  (tramper/total)^2 + (tanker/total)^2) %>% 
  dplyr::mutate(share_liner = liner/total) %>% 
  dplyr::mutate(share_special = special/total) %>%
  dplyr::mutate(share_tramper = tramper/total) %>% 
  dplyr::mutate(share_tanker = tanker/total) %>% 
  dplyr::filter(total > 0)

## (Table 3) Summary statistics for independent variables. ----

# construct summary statistics table
denom = 1e6
temp <- 
  data_english %>%
  # rescale
  dplyr::mutate(tonnage_size = total/denom,
                tonnage_size_liner = liner/denom,
                tonnage_size_special = special/denom,
                tonnage_size_tramper = tramper/denom,
                tonnage_size_tanker = tanker/denom) %>% 
  dplyr::select(tonnage_size, 
                tonnage_size_liner,
                tonnage_size_special,
                tonnage_size_tramper,
                tonnage_size_tanker, 
                share_liner, 
                share_special, 
                share_tramper,
                share_tanker,
                HHI) %>% 
  dplyr::rename(type_HHI = HHI) %>% 
  tidyr::gather("variable", "value") %>%
  dplyr::group_by(variable) %>%
  dplyr::summarise(
    N = sum(!is.na(value)),
    mean = mean(value, na.rm = TRUE),
    sd = sd(value, na.rm = TRUE),
    min = min(value, na.rm = TRUE),
    q25 = quantile(value, probs = .25, na.rm = TRUE),
    q50 = quantile(value, probs = .5, na.rm = TRUE),
    q75 = quantile(value, probs = .75, na.rm = TRUE),
    max = max(value, na.rm = TRUE)) %>%
  dplyr::mutate(
    mean = format(round(mean, 3), 3),
    sd = format(round(sd, 3), 3),
    min = format(round(min, 3), 3),
    q25 = format(round(q25, 3), 3),
    q50 = format(round(q50, 3), 3),
    q75 = format(round(q75, 3), 3),
    max = format(round(max, 3), 3)) %>% 
  dplyr::select(-variable) %>% 
  as.data.frame()
# reorder variables
temp <- rbind(
  temp[5,],
  temp[6:9,],
  temp[1:4,],
  temp[10,])
rownames(temp) <-
  c("total tonnage size",
    "total tonnage size of liner",
    "total tonnage size of special",
    "total tonnage size of tanker",
    "total tonnage size of tramper",
    "share of liner",
    "share of special",
    "share of tanker",
    "share of tramper",
    "HHI based on carrier types")
# make a latex table
table_name <-
  "figuretable/shipping_market_stats_table.tex"
x <- Hmisc::latex(
  file = table_name,
  title = "",
  booktabs = TRUE,
  object = temp,
  col.just = rep("r", dim(temp)[2]),
  center = "none",
  table.env = FALSE,
  cgroupTexCmd = "itshape",
  colnamesTexCmd = "itshape",
  collabel.just = rep("r", dim(temp)[2]),
  rgroupTexCmd = "itshape",
  n.rgroup = c(5,5),
  rgroup = c("Size variables (million (D/W) tons)", 
             "Specialization variables (percentage share)")
)

## (Table 4) Regression ----

### construct data for preliminary regression ----
buyer_id <- 
  data_english$id
target_id <-
  data_english$id
id_new <- 
  expand.grid(buyer_id,
              target_id)
colnames(id_new) <-
  c("buyer_id",
    "target_id")
data_english_for_reg <- 
  id_new %>% 
  dplyr::full_join(data_english,
                   by = c("buyer_id" = "id")) %>%
  dplyr::full_join(data_english,
                   by = c("target_id" = "id")) %>% 
  dplyr::filter(buyer_id != target_id) %>% 
  dplyr::mutate(match = 
                  ifelse(group.x == group.y, 1, 0)) %>% 
  dplyr::mutate(Fuji_xy =
                  ifelse(Fuji.x == Fuji.y, 1, 0)) %>% 
  dplyr::mutate(Kogin_xy =
                  ifelse(Kogin.x == Kogin.y, 1, 0)) %>% 
  dplyr::mutate(Chogin_xy =
                  ifelse(Chogin.x == Chogin.y, 1, 0)) %>% 
  dplyr::mutate(Yasuda_xy =
                  ifelse(Yasuda.x == Yasuda.y, 1, 0)) %>% 
  dplyr::mutate(Daiichi_xy = 
                  ifelse(Daiichi.x == Daiichi.y, 1, 0)) %>% 
  dplyr::mutate(Mitsui_xy = 
                  ifelse(Mitsui.x == Mitsui.y, 1, 0)) %>% 
  dplyr::mutate(MitsuiShintaku_xy =
                  ifelse(MitsuiShintaku.x == MitsuiShintaku.y, 1, 0)) %>% 
  dplyr::mutate(Sumitomo_xy = 
                  ifelse(Sumitomo.x == Sumitomo.y, 1, 0)) %>% 
  dplyr::mutate(SumitomoShintaku_xy =
                  ifelse(SumitomoShintaku.x == SumitomoShintaku.y, 1, 0)) %>% 
  dplyr::mutate(Mitsubishi_xy =
                  ifelse(Mitsubishi.x == Mitsubishi.y, 1, 0)) %>% 
  dplyr::mutate(Kangin = 
                  ifelse(Kangin.x == Kangin.y, 1, 0)) %>% 
  dplyr::mutate(MitsubishiShintaku_xy = 
                  ifelse(MitsubishiShintaku.x == MitsubishiShintaku.y, 1, 0)) %>% 
  dplyr::mutate(Kobe_xy = 
                  ifelse(Kobe.x == Kobe.y, 1, 0)) %>% 
  dplyr::mutate(Daiwa_xy = 
                  ifelse(Daiwa.x == Daiwa.y, 1, 0)) %>% 
  dplyr::mutate(Sanwa_xy = 
                  ifelse(Sanwa.x == Sanwa.y, 1, 0)) %>%
  dplyr::mutate(Tokai_xy = 
                  ifelse(Tokai.x == Tokai.y, 1, 0)) %>% 
  dplyr::mutate(YasudaShintaku_xy = 
                  ifelse(YasudaShintaku.x == YasudaShintaku.y, 1, 0)) %>%
  dplyr::mutate(others_xy = 
                  ifelse(others.x == others.y, 1, 0)) 
data_english_for_reg <- 
  data_english_for_reg %>% 
  replace(is.na(.), 0) %>%
  dplyr::mutate(bank_coverage = 
                  Fuji_xy + Kogin_xy +
                  Chogin_xy + Yasuda_xy + 
                  Daiichi_xy + Mitsui_xy +
                  MitsuiShintaku_xy + Sumitomo_xy +
                  SumitomoShintaku_xy +
                  Mitsubishi_xy + Kangin + 
                  MitsubishiShintaku_xy + Kobe_xy +
                  Daiwa_xy + Sanwa_xy + Tokai_xy +
                  YasudaShintaku_xy + others_xy) %>% 
  # the number of all banks is 19
  dplyr::mutate(bank_coverage_ratio = bank_coverage/19) %>% 
  dplyr::mutate(same_type = 
                  ifelse(firm_type.x == firm_type.y, 1, 0))
### five specifications ----
res1 <- 
  glm(match ~ 
        log(liner.x*liner.y+1) +
        log(tramper.x*tramper.y+1) +
        log(special.x*special.y+1) +
        log(tanker.x*tanker.y+1) +
        log(total.x*total.y+1),
      data = data_english_for_reg, 
      family = "binomial")
res2 <-
  glm(match ~ bank_coverage_ratio +
        log(HHI.x*HHI.y+1) + 
        log(share_liner.x*share_liner.y+1) +
        log(share_special.x*share_special.y+1) +
        log(share_tramper.x*share_tramper.y+1) +
        log(share_tanker.x*share_tanker.y+1),
      data = data_english_for_reg,
      family = "binomial")
res3 <- 
  glm(match ~ log(liner.x*liner.y+1) +
        log(tramper.x*tramper.y+1) +
        log(special.x*special.y+1) + 
        log(tanker.x*tanker.y+1) +
        log(total.x*total.y+1) +
        bank_coverage_ratio + 
        log(HHI.x*HHI.y+1) +
        log(share_liner.x*share_liner.y+1) +
        log(share_special.x*share_special.y+1) + 
        log(share_tramper.x*share_tramper.y+1) +
        log(share_tanker.x*share_tanker.y+1),
      data = data_english_for_reg,
      family = "binomial")
res4 <-
  glm(match ~ log(liner.x*liner.y+1) +
        log(tramper.x*tramper.y+1) +
        log(special.x*special.y+1) +
        log(tanker.x*tanker.y+1) +
        log(total.x*total.y+1) +
        bank_coverage_ratio +
        log(HHI.x*HHI.y+1) + 
        log(share_liner.x*share_liner.y+1) +
        log(share_special.x*share_special.y+1) + 
        log(share_tramper.x*share_tramper.y+1) +
        log(share_tanker.x*share_tanker.y+1) +
        same_type,
      data = data_english_for_reg, 
      family = "binomial")
res_LPM <- 
  glm(match ~ log(liner.x*liner.y+1) +
        log(tramper.x*tramper.y+1) +
        log(special.x*special.y+1) + 
        log(tanker.x*tanker.y+1) +
        log(total.x*total.y+1) +
        bank_coverage_ratio +
        log(HHI.x*HHI.y+1) +
        log(share_liner.x*share_liner.y+1) +
        log(share_special.x*share_special.y+1) + 
        log(share_tramper.x*share_tramper.y+1) +
        log(share_tanker.x*share_tanker.y+1) +
        same_type,
      data = data_english_for_reg)
table_name <-
  "figuretable/regression_matching.tex"
x <-
  stargazer(res1,
            res2,
            res3,
            res4,
            res_LPM,
            header = FALSE,
            float = FALSE,
            type = "latex",
            digits = 3,
            dep.var.labels = c("1(match)"),
            covariate.labels = c("log(liner$_{b}$ *liner$_{t}$+1)",
                                 "log(tramper$_{b}$ *tramper$_{t}$+1)",
                                 "log(special$_{b}$ *special$_{t}$+1)",
                                 "log(tanker$_{b}$ *tanker$_{t}$+1)",
                                 "log(total$_{b}$ *total$_{t}$+1)",
                                 "bank coverage similarity ratio",
                                 "log(HHI$_{b}$ *HHI$_{t}$+1)",
                                 "log(share of liner$_{b}$ *share of liner$_{t}$+1)",
                                 "log(share of special$_{b}$ *share of special$_{t}$+1)",
                                 "log(share of tramper$_{b}$ *share of tramper$_{t}$+1)",
                                 "log(share of tanker$_{b}$ *share of tanker$_{t}$+1)", 
                                 "same type","Intercept"),
            add.lines = list(c("Model",
                               "Logit", 
                               "Logit",
                               "Logit",
                               "Logit",
                               "OLS")),
            omit.stat = c("f", "ll"),
            star.cutoffs = NA,
            omit.table.layout = "n",
            model.names = FALSE,
            out = table_name
  )


## (Table 2) Summary of total tonnage size for each group ----

data_english <-
  data_english[,-c(9:26)]
data_english_for_reg <-
  data_english %>% 
  dplyr::mutate(type = recode(type, 
                              "合併会社" = "(1) main",
                              "系列会社" = "(2) affiliate",
                              "専属会社" = "(3) wholly controlled",
                              "グループ外" = "unmatched"
  )) %>% 
  dplyr::select(-name,
                -firm_type) %>% 
  dplyr::select(id, 
                type,
                group, 
                everything()) %>% 
  dplyr::group_by(group) %>% 
  dplyr::mutate(group_id =
                  cur_group_id()) %>% 
  dplyr::mutate(group_total_tonnage = 
                  sum(total)) %>% 
  dplyr::ungroup() %>% 
  dplyr::arrange(group_id)

data_english_for_reg_table <- 
  data_english_for_reg %>% 
  dplyr::filter(total > 0) %>% 
  dplyr::group_by(group, 
                  type) %>% 
  dplyr::summarise(total_tonnage = sum(total),
                   num_of_firms = dplyr::n())
total_tonnage <-
  data_english_for_reg_table %>% 
  dplyr::group_by(group) %>% 
  dplyr::summarise(total_tonnage_group = sum(total_tonnage))
data_english_for_reg_table <- 
  cbind(data_english_for_reg_table,
        c(total_tonnage$total_tonnage_group[1], "", "",
          total_tonnage$total_tonnage_group[2], "", "",
          total_tonnage$total_tonnage_group[3], "", "",
          total_tonnage$total_tonnage_group[4], "", "",
          total_tonnage$total_tonnage_group[5], "", "",
          total_tonnage$total_tonnage_group[6], "", "",
          total_tonnage$total_tonnage_group[7]))
colnames(data_english_for_reg_table) <-
  c("group",
    "firm type", 
    "total tonnage",
    "number of firms",
    "total tonnage in a group")
data_english_for_reg_table <- 
  data_english_for_reg_table[,-1]
data_english_for_reg_table <-
  rbind(data_english_for_reg_table,
        c("",sum(data_english_for_reg_table[,2]),
          sum(data_english_for_reg_table[,3]),
          sum(data_english_for_reg_table[,2])))
table_name <-
  "figuretable/total_tonnage_size_group_table.tex"
x <- Hmisc::latex(
  file = table_name,
  title = "",
  booktabs = TRUE,
  object = data_english_for_reg_table,
  col.just = c("l", "r", "r", "r"),
  center = "none",
  table.env = FALSE,
  #cgroupTexCmd = "itshape",
  #colnamesTexCmd = "itshape",
  #collabel.just = rep("r", dim(data_english_for_reg_table)[2])#,
  rgroupTexCmd = "itshape",
  n.rgroup = c(3,3,3,3,3,3,1, 1),
  rgroup = c("Nippon Yusen", 
             "Mitsui OSK Line",
             "Japan Line",
             "Kawasaki Kisen Kaisha",
             "Yamashita Shinnihon Kisen",
             "Showa Line",
             "Unmatched",
             "Total")
)


# Export dataset for maximum rank estimator ----

denom = 1e6
data_english_for_reg_for_maximum_rank_estimation <- 
  data_english_for_reg %>% 
  # rescale
  dplyr::mutate(liner = liner/denom,
                special = special/denom,
                tramper = tramper/denom,
                tanker = tanker/denom,
                total = total/denom,
                group_total_tonnage = group_total_tonnage/denom) 
write.csv(x=data_english_for_reg_for_maximum_rank_estimation,
          file="output/data_for_maximum_rank_estimation.csv")


