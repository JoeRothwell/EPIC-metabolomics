# All functions and workflow for alcohol study

# Intake correlation for CS and HCC study and match function
# Get functions
source("Intake_cor_functions.R")
# Approach 1: match HCC and CS features and perform separate correlations on the overlap ----

matched_pos <- CS_HCC_match(mode = "pos")
matched_neg <- CS_HCC_match(mode = "neg")

# Extract vectors of matched features
CSpos <- unique(matched_pos$CS_feat)
CSneg <- unique(matched_neg$CS_feat)
HCCpos <- unique(matched_pos$HCC_feat)
HCCneg <- unique(matched_neg$HCC_feat)

# Run correlations on matched features for both CS and HCC
alcpos_cs <- intakecor_cs(food = "Qe_Alc", incr = F, min.sample = 250, pcutoff = 1, matchvec = CSpos)
alcneg_cs <- intakecor_cs(food = "Qe_Alc", incr = F, pos = F, min.sample = 250, pcutoff = 1, matchvec = CSneg)
alcpos_hcc <- intakecor_hcc(matchvec = HCCpos)
alcneg_hcc <- intakecor_hcc(pos = F, matchvec = HCCneg)

# Join pos and neg data together for each study
pos_cs <- alcpos_cs %>% select(Pcor:CS_feat)
pos_hcc <- alcpos_hcc %>% select(Pcor:HCC_feat)
matched_pos <- matched_pos %>% left_join(pos_cs, by = "CS_feat") %>% left_join(pos_hcc, by = "HCC_feat")

neg_cs <- alcneg_cs %>% select(Pcor:CS_feat)
neg_hcc <- alcneg_hcc %>% select(Pcor:HCC_feat)
matched_neg <- matched_neg %>% left_join(neg_cs, by = "CS_feat") %>% left_join(neg_hcc, by = "HCC_feat")

# Plot correlations (4165 in both modes)
plot(matched_pos$Pcor.x, matched_pos$Pcor.y, xlab = "Correlation CS", ylab = "Correlation HCC")
points(matched_neg$Pcor.x, matched_neg$Pcor.y, xlab = "Correlation CS", ylab = "Correlation HCC")
# See Feature_matching for ggplot plot




# Approach 2: Match significant CS features with all HCC ----

# Get significant CS features and extract feature table number
pos_cs <- intakecor_cs(food = "Qe_Alc", incr = F, min.sample = 250, pcutoff = 0.05, matchvec = NULL)
neg_cs <- intakecor_cs(food = "Qe_Alc", incr = F, pos = F, min.sample = 250, pcutoff = 0.05, matchvec = NULL)
filt1 <- pos_cs$feat
filt2 <- neg_cs$feat

# Match with all HCC features and get vector of matched features
matched_pos_filt <- CS_HCC_match(mode = "pos", filt = filt1)
matched_neg_filt <- CS_HCC_match(mode = "neg", filt = filt2)
f1 <- matched_pos_filt$HCC_feat
f2 <- matched_neg_filt$HCC_feat

# Run correlations on subset
pos_hcc <- intakecor_hcc(matchvec = f1)
neg_hcc <- intakecor_hcc(pos = F, matchvec = f2)

# Join model data for CS and HCC
pos_cs <- pos_cs %>% select(feat, Pcor:pval)
neg_cs <- neg_cs %>% select(feat, Pcor:pval)
pos_hcc <- pos_hcc %>% select(Pcor:HCC_feat)
neg_hcc <- neg_hcc %>% select(Pcor:HCC_feat)

matched_pos_filt1 <- matched_pos_filt %>% left_join(pos_cs, by = c("CS_feat" = "feat")) %>% 
  left_join(pos_hcc, by = "HCC_feat", suffix = c(".cs", ".hcc"))
matched_neg_filt1 <- matched_neg_filt %>% left_join(neg_cs, by = c("CS_feat" = "feat")) %>% 
  left_join(neg_hcc, by = "HCC_feat", suffix = c(".cs", ".hcc"))



