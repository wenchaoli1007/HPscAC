
library(Seurat)
library(dplyr)
library(readr)
library(HPscAC)
library(glmnet)
library(purrr)

setwd("~/HPscAC/data/")

model_set = readRDS("./all_model.RDS")
feature_set = readRDS("./all_model_inputfeatures.RDS")
input = readRDS("./demo_bcg_cd8t.RDS")

cell_type = "CD8T"
model = model_set[[cell_type]]
marker_gene = feature_set[[cell_type]]

## preprocess the data
preprocessed_df = PreProcess(input, cell_type, model, marker_gene)


### predict the age

predict_res = AgingClockCalculator(preprocessed_df, model, marker_gene)

age_per_donor = Age_Donor(predict_res)

### visualization

# density plot #
ggplot(predict_res, aes(x=Prediction, fill=condition)) +
  geom_density(alpha=0.5)+
  theme_classic() +
  facet_wrap(~group, ncol = 4) +
  theme(text = element_text(size = 16)) +
  geom_vline(aes(xintercept = md, color = condition), size = 0.8) +
  xlab("Predicted Age") +
  xlim(0, 100) +
  scale_fill_manual(values = c("#b0d992", "#f4c28f", "#d2352c")) +
  scale_color_manual(values = c("#b0d992", "#f4c28f", "#d2352c"))

# boxplot #

test_plot = predict_res %>% group_by(donor_id) %>% mutate(md = round(mean(Prediction))) %>% dplyr::select(-Prediction) %>% distinct()

idx = match(test_plot$donor_id, time_data$sample)
test_plot$days = time_data$time[idx]
test_plot$donor = time_data$donor[idx]

test_plot4 = test_plot %>% filter(condition %in% c("control"))
adj_diff = lm(test_plot4$md ~ test_plot4$age)

adj_age = as.numeric(coef(adj_diff)[1]) + test_plot$age * as.numeric(coef(adj_diff)[2])
test_plot$adj_age = adj_age
test_plot$diff = test_plot$md - test_plot$adj_age

test_plot$diff = round(test_plot$diff, 3)
test_plot$condition = factor(test_plot$condition, 
                             levels = c("control", "mild", "severe"))

test_plot$celltype = cell_type
bp_df = test_plot

pers = bp_df %>% filter(condition != "control")

ggplot(bp_df, aes(x = condition, y = diff, fill = condition)) + 
  geom_boxplot(alpha = 0.6) + 
  geom_jitter(color="black", size=0.4, alpha=0.5) +
  theme_classic() +
  theme(axis.title.x = element_text(size = 13),
        axis.title.y = element_text(size = 13),
        strip.text = element_text(size = 13),
        axis.text.y = element_text(size = 13)) +
  facet_wrap(~celltype, ncol = 5) + 
  ylab("TAA") +
  ylim(-25, 82) +
  stat_compare_means(comparisons = my_comparisons, label = "p.signif", size=4)+
  scale_fill_manual(values = c("#b0d992", "#f4c28f", "#d2352c")) +
  theme(axis.text.x = element_text(size = 0), axis.ticks.x = element_blank()) +
  xlab("Condition")



