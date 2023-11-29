ancova_post_hoc <- function(file = "../../HKFI-FoodIntake/code/ABCdata.tsv",
                           name = "Figure1"){
  library(multcomp)
  library(emmeans)
  require(tidyverse)
  require(rstatix)
  require(ggpubr)
  require(broom)
  require(patchwork)
  library(interactions)
  
  dataABC = read_tsv(file)


  # ========================== Assumption check ======================
  # TYPE 2 sum of square is appropriate for the model with no signifcant interaction
  
  path = paste0("../../HKFI-FoodIntake/results/",name)
  
  if (!dir.exists(path)){
    dir.create(path, recursive = T)
  }
  
  
  # == check 1
  step1_anova = dataABC %>% anova_test(formula = Bodyweight ~ Group)
  
  if (  step1_anova$p < 0.05 ){
    print ("======== warning: The covariate (body weight) and the treatment (Group) are dependent ======")
  }

  # == check 2
  step2_anova = dataABC %>% anova_test(formula = foodIntake ~ Group * Bodyweight,type = 2)
  
  if (  step2_anova[3,"p"] < 0.05 ){
    print ("======== warning: The covariate (body weight) and the treatment (Group) have interaction. We will use III SS======")
    method = "III"
  }else{
    method = "II"
  }
  
  # == check 3
  
  # Fit the model, the covariate goes first
  model <- lm(foodIntake ~  Bodyweight + Group, data = dataABC)
  # Inspect the model diagnostic metrics
  model.metrics <- augment(model)  # Remove details
  # Assess normality of residuals using shapiro wilk test
  step3_shapiro = shapiro_test(model.metrics$.resid)
  
  if (  step3_shapiro$p.value < 0.05 ){
    print ("======== warning: we can not assume normality of residuals based on The Shapiro Wilk test ======")
  }
  
  # == check 4
  
  step4_levene = model.metrics %>% levene_test(.resid ~ Group)
  
  if (  step3_shapiro$p.value < 0.05 ){
    print ("======== warning: we can not assume homogeneity of the residual variances for all groups based on The Levene’s test  ======")
  }

  # == check 5
  
  # which model is better?
  # Fit the model, the covariate goes first
  lm_model1 <- lm(foodIntake ~  Bodyweight * Group, data = dataABC)
  result_1 = summary(lm_model1)
  
  lm_model2 <-  update( lm_model1, . ~ Bodyweight + Group) 
  result_2 = summary(lm_model2)
  
  if (  result_2$adj.r.squared > result_1$adj.r.squared ){
    print(result_1$adj.r.squared )
    print(result_2$adj.r.squared )
    print ("======== warning: The model without interaction may work better according to the R squared ======")
  }
  
  # == check 6
  
  if (  result_1$coefficients[2,4] > 0.05){

    print ("======== warning: Body weight have no significant effect on food intake =====")
  }
  
  
  print ("======== End of assumption check  ======")
  # ========================== Computation ======================

  if (method == "III"){
    
    # This is slightly more involved than the type II results. 
    # First, it is necessary to set the contrasts option in R. 
    # Because the multi-way ANOVA model is over-parameterised, it is necessary to choose a contrasts setting that sums to zero, 
    # otherwise the ANOVA analysis will give incorrect results with respect to the expected hypothesis. (The default contrasts type does not satisfy this requirement.)
    
    options(contrasts = c("contr.sum","contr.poly"))
    
    # res.aov = rstatix::Anova(mod = ancova_model, type = "III")
    dataABC$Group = factor(dataABC$Group)
    ancova_model <- aov(foodIntake ~ Group * Bodyweight, data = dataABC)
    res.aov = rstatix::Anova(mod = ancova_model, type = "III")
  
    res.aov <- dataABC %>% anova_test( foodIntake ~  Bodyweight + Group + Bodyweight:Group, type = 3)
    
    # ===== if we found interaction effect, 
    # then we will use simple effect analysis to further see the difference in certain level of body weight ===
    covariate_level = quantile(dataABC$Bodyweight, probs = c(0.25,0.5,0.75))
    
    #covariate_level = quantile(dataABC$Bodyweight, probs = seq(.05, .95, by = .1))
    
    emm <- emmeans(ancova_model, ~ Group | Bodyweight, at = list(Bodyweight = covariate_level)) 
    # summary(emm) 

    # manipulate the result
    res.emmeans <- emm %>% tibble::as_tibble() %>% 
      dplyr::rename(se = "SE", conf.low = "lower.CL", conf.high = "upper.CL") %>% 
      mutate(method = "Emmeans test")
     
    # pairwise_comparisons <- pairwise_comparison_after_emmeans(emm)
    # by default, pairs() function use "tukey", we will use "bonferroni" in our project
    # res.compare <- summary(pairs(emm), adjust = "bonferroni") %>% as_tibble()
  
    comparisons <- emmeans::contrast(emm, method = "pairwise", adjust = "bonferroni")
    comparisons <- tidy(comparisons, conf.int = TRUE, conf.level = 0.95)
    res.comparisons <- comparisons %>% tidyr::separate(col = "contrast", 
                                                   into = c("group1", "group2"), sep = " - ") %>% 
      add_significance("adj.p.value")
                                                                                                                           
    to.remove <- c("estimate", "std.error", "conf.low", "conf.high", "null.value")
    to.keep <- setdiff(colnames(res.comparisons), to.remove)
    res.comparisons  <- res.comparisons [, to.keep]
    
    # merge emmean 
    res.anno = res.emmeans %>% dplyr::select("Group","Bodyweight","emmean")
    
    res.comparisons <- res.comparisons %>% left_join(res.anno,by = c("group1" = "Group","Bodyweight")) %>% 
      dplyr::rename(emmean_group1 = emmean ) %>% 
      left_join(res.anno,by = c("group2" = "Group", "Bodyweight")) %>% 
      dplyr::rename(emmean_group2 = emmean )
    
    # === use Johnson-Neyman intervals ===
    # USE INTERACTIONS to find the interval
    # Convert the factor variable 'Group' into numeric dummy variables
    dummy_vars <- model.matrix(~ Group - 1, dataABC)
    # Add the new dummy variables to your data frame
    new_data <- cbind(dataABC, dummy_vars)
    
    points = dataABC %>% mutate(group = Group)
    if (is.factor(points[,"group"])){
      group_list = unique(points[,"group"]) %>% unlist()%>% as.vector()
    }else{
      group_list = unique(points[,"group"]) %>% unlist() %>% as.vector()
    }
    possible_combinations = sets::set_power(group_list) %>% as.vector() 
    compared_group = list()
    for (i in seq(length(possible_combinations))){
      result_c = possible_combinations[[i]] %>% as.character()
      if (length(result_c) == length(group_list)-2){
        sub_compared_group = group_list[!group_list %in% result_c]
        compared_group[[length(compared_group)+1]] =  sub_compared_group 
      }
    }

    # formula 
    
    stat.test = data.frame(group1 = character(), group2 = character(),pval_or_low = character(), psign_or_high = character(), interactions = character())
    for (x in compared_group){
      part_data = new_data %>% filter(Group %in% x)
      part_data$Group = factor(part_data$Group,levels = x)

      library(interactions)
      # Update the model using the new data frame with numeric dummy variables
      part_anova = part_data %>% anova_test(formula = foodIntake ~ Group * Bodyweight,type = 2)
      print (x)
      # ===== check if the interaction also stay siginificant in two groups ====
      if (   part_anova[3,"p"] < 0.05 ){
        print ("======== warning: The covariate (body weight) and the treatment (Group) have interaction======")
        # rename the group name we want to 
        # we will use the first one as the referecen, e.g. c(A,B), A will be the reference
        target_name = paste0("Group",x[[2]])
        part_data <- part_data %>% rename(dum_group = `target_name`)
        
        ancova_model_updated <- lm(foodIntake ~ Bodyweight * dum_group, data = part_data)
        # Perform the region of significance analysis using the new model
        region_of_significant <- interactions::sim_slopes(ancova_model_updated, 
                                                          pred = dum_group, 
                                                          modx = Bodyweight,    
                                                          johnson_neyman = TRUE,
                                                          control.fdr = TRUE,
                                                          jnalpha = 0.05)
        res.region = region_of_significant$jn[[1]]
        
        # JN interval
        res.region$bounds
        
        # JN plot   
        p_report3 = res.region$plot
        
        stat.test[nrow(stat.test) + 1,] = c(x[1],x[2],res.region$bounds[1],res.region$bounds[2],"interaction")  
        
        ggsave(p_report3, filename = paste0("../../HKFI-FoodIntake/results/",name,"/ABC_JN_interval.png"), dpi = 320)
        
        
      }else{

        print ("======== notification: The covariate (body weight) and the treatment (Group) have no interaction======")
        # We will use normal TukeyHSD to compare
        part_ancova_model = aov(foodIntake ~ Bodyweight + Group + Bodyweight:Group,data = part_data)
        
        #posthoc = glht(part_ancova_model, linfct = mcp( Group = "Tukey"))
        #res.comparisons = tidy(summary(posthoc)) %>% tidyr::separate(col = "lhs", 
        #                                                             into = c("group1", "group2"), sep = " - ") %>% 
        #  add_significance("p.value") %>% select(-"rhs")
        
        res.comparisons_part <-  part_data %>% 
          emmeans_test(
            foodIntake ~ Group, covariate = Bodyweight,
            p.adjust.method = "bonferroni"
          )
        
        res.emmeans_part =  get_emmeans(res.comparisons_part)
  
        
        stat.test[nrow(stat.test) + 1,] = c(x[1],x[2],res.comparisons_part$p.adj,res.comparisons_part$p.adj.signif,"no interaction")  
      
        
      }
    }
    
    write_tsv(stat.test,paste0("../../HKFI-FoodIntake/results/",name,"/ABC_pairwise_interaction_analysis_result.tsv"))

   }else if (method == "II"){
    stat.test= NULL
    # ==== identical with those in rstatix
    #dataABC$Group = factor(dataABC$Group)
    #ancova_model <- aov(data = dataABC, foodIntake ~ Bodyweight + Group) 
    #res.aov = rstatix::Anova(mod = ancova_model, type = "II")
    
    res.aov <- dataABC %>% anova_test( foodIntake ~  Bodyweight * Group, type = 2)
    # posthoc = glht(ancova_model, linfct = mcp( Group = "Tukey"))
    # res.comparisons = tidy(summary(posthoc)) 
    
    # === we get the marginal mean by use the median bodyweight 
    # covariate_level = quantile(dataABC$Bodyweight, probs = c(0.5))
    
    # emm <- emmeans(ancova_model, ~ Group | Bodyweight, at = list(Bodyweight = covariate_level)) 

    # res.emmeans <- emm %>% tibble::as_tibble() %>% 
    #  dplyr::rename(se = "SE", conf.low = "lower.CL", conf.high = "upper.CL") %>% 
    #  mutate(method = "Emmeans test")
    
    res.comparisons <- dataABC %>% 
      emmeans_test(
        foodIntake ~ Group, covariate = Bodyweight,
        p.adjust.method = "bonferroni"
      )
    
    res.emmeans =  get_emmeans(res.comparisons)
    
  }

  write_tsv(get_anova_table(res.aov),paste0("../../HKFI-FoodIntake/results/",name,"/ABC_ancova_result.tsv"))
  write_tsv(res.comparisons ,paste0("../../HKFI-FoodIntake/results/",name,"/ABC_pairwise_comparison_result.tsv"))
  write_tsv(res.emmeans ,paste0("../../HKFI-FoodIntake/results/",name,"/ABC_adjusted_mean_result.tsv"))
  
  print ("======== End of computation ======")
  
  # ========================== Final visualization ======================

  # Pairwise comparisons through emmeans
  # simple effect analysis at the level of mean body weight 
  library(emmeans)
  pwc <- dataABC %>% 
    emmeans_test(
      foodIntake ~ Group, covariate = Bodyweight,
      p.adjust.method = "bonferroni"
    )
  
  # get_emmeans(pwc)
  pwc <- pwc %>% add_xy_position(x = "Group", fun = "mean_se")
  
  p_report2 = ggline(get_emmeans(pwc), x = "Group", y = "emmean") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
    stat_pvalue_manual(pwc, hide.ns = TRUE, tip.length = FALSE) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE, row = 2 ),
      caption = get_pwc_label(pwc),
      y = "Least-squares means"
    )
  
  
  ggsave(p_report2, filename = paste0("../../HKFI-FoodIntake/results/",name,"/ABC_Least_squares_means.png"), dpi = 320)
  
  
  
  plot1 = ggscatter(
    dataABC, x = "Bodyweight", y = "foodIntake", 
    add = "reg.line",size = 3, color = "Group")+
    geom_point(aes(color = Group), size = 3) + 
    geom_point(shape = 1, color = "black", size = 3) + 
    stat_regline_equation(
      aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Group)
    )  +
    geom_rug()+ 
    theme_pubr() + 
    labs(y = "Food intake (g)", x = "Body weight (g)") + 
    ggsci::scale_color_lancet() +
    theme(legend.position = "bottom") 
  
  # if we have interaction effect, we will draw the JN interval 
  if (!is.null(stat.test) & length(unique(stat.test$interactions)) == 2 & table(stat.test$interactions)[1] == 1 ){
    upper_bound = stat.test %>% filter(interactions == "interaction") %>% pull("psign_or_high") %>% as.numeric()
    lower_bound = stat.test %>% filter(interactions == "interaction") %>% pull("pval_or_low") %>% as.numeric()
    min_bound  = min(dataABC$Bodyweight)
    max_bound  = max(dataABC$Bodyweight)
    
    y_min_food_intake = max(dataABC$foodIntake) + 0.03* max(dataABC$foodIntake)
    y_max_food_intake = max(dataABC$foodIntake) + 0.02* max(dataABC$foodIntake)
    plot1 <- plot1 +
      annotate(
        xmin =  c(-Inf,upper_bound),     # plus one 
        xmax = c(upper_bound,Inf),
        ymin =  y_min_food_intake,
        ymax =  y_max_food_intake,
        geom = "rect",
        fill = c("red", "blue")
      )
  }
  
  
  dens1 <- ggplot(dataABC, aes(x = Bodyweight, fill = Group)) + 
    geom_density(alpha = 0.4) + 
    theme_void() + 
    theme(legend.position = "none")+ 
    ggsci::scale_fill_lancet()
  
  dens2 <- ggplot(dataABC, aes(x = foodIntake, fill = Group)) + 
    geom_density(alpha = 0.4) + 
    theme_void() + 
    theme(legend.position = "none") + 
    coord_flip()+ 
    ggsci::scale_fill_lancet()
  
  
  p_report1  = dens1 + plot_spacer() + plot1 + dens2 + 
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
  
  ggsave(p_report1 , filename = paste0("../../HKFI-FoodIntake/results/",name,"/ABC_scatter_with_density.png"), dpi = 320)
  
}















pairwise_comparison_after_emmeans <- function (res.emmeans, 
                                               grouping.vars = NULL, 
                                               method = "pairwise", 
                                               p.adjust.method = "bonferroni", conf.level = 0.95) 
{
  comparisons <- emmeans::contrast(res.emmeans, by = grouping.vars, 
                                   method = method, adjust = "none")
  comparisons <- tidy(comparisons, conf.int = TRUE, conf.level = conf.level)
  comparisons <- comparisons %>% tidyr::separate(col = "contrast", 
                                                 into = c("group1", "group2"), sep = "-") %>% dplyr::rename(se = "std.error", 
                                                                                                            p = "p.value") %>% dplyr::select(!!!syms(grouping.vars), 
                                                                                                                                             everything())
  p.adjusted <- emmeans::contrast(res.emmeans, by = grouping.vars, 
                                  method = method, adjust = p.adjust.method) %>% as.data.frame() %>% 
    pull("p.value")
  comparisons <- comparisons %>% mutate(p.adj = p.adjusted) %>% 
    add_significance("p.adj")
  res.emmeans.tbl <- tibble::as_tibble(res.emmeans)
  variables <- intersect(colnames(res.emmeans.tbl), colnames(comparisons))
  for (variable in variables) {
    if (is.factor(res.emmeans.tbl[[variable]])) {
      comparisons[[variable]] <- factor(comparisons[[variable]], 
                                        levels = levels(res.emmeans.tbl[[variable]]))
    }
  }
  comparisons <- base::droplevels(comparisons)
  comparisons %>% dplyr::arrange(!!!syms(grouping.vars))
  return(comparisons)
}

ancova_post_hoc_EE <- function(file = "../../HKFI-FoodIntake/code/ABCdata.tsv",
                            name = "Figure1"){
  library(multcomp)
  library(emmeans)
  require(tidyverse)
  require(rstatix)
  require(ggpubr)
  require(broom)
  require(patchwork)
  library(interactions)
  
  dataABC = read_csv(file)
  
  dataABC <- dataABC %>% mutate(Group = paste0("Genotype ",Group))
  
  # ========================== Assumption check ======================
  # TYPE 2 sum of square is appropriate for the model with no signifcant interaction
  
  path = paste0("../../HKFI-FoodIntake/results/",name)
  
  if (!dir.exists(path)){
    dir.create(path, recursive = T)
  }
  
  
  # == check 1
  step1_anova = dataABC %>% anova_test(formula = Bodyweight ~ Group)
  
  
  if (  step1_anova$p < 0.05 ){
    print ("======== warning: The covariate (body weight) and the treatment (Group) are dependent ======")
  }
  
  # == check 2
  step2_anova = dataABC %>% anova_test(formula = EE ~ Group * Bodyweight,type = 2)

  # homogeneity of regression slopes
  if (  step2_anova[3,"p"] < 0.05 ){
    print ("======== warning: The covariate (body weight) and the treatment (Group) have interaction. We will use III SS======")
    method = "III"
  }else{
    print ("======== warning: The covariate (body weight) and the treatment (Group) have no interaction. We will use II SS======")
    method = "II"
  }
  
  # == check 3
  
  # Fit the model, the covariate goes first
  model <- lm(EE ~  Bodyweight + Group, data = dataABC)
  # Inspect the model diagnostic metrics
  model.metrics <- augment(model)  # Remove details
  # Assess normality of residuals using shapiro wilk test
  step3_shapiro = shapiro_test(model.metrics$.resid)
  
  if (  step3_shapiro$p.value < 0.05 ){
    print ("======== warning: we can not assume normality of residuals based on The Shapiro Wilk test ======")
  }
  
  # == check 4
  
  step4_levene = model.metrics %>% levene_test(.resid ~ Group)
  
  if (  step3_shapiro$p.value < 0.05 ){
    print ("======== warning: we can not assume homogeneity of the residual variances for all groups based on The Levene’s test  ======")
  }
  
  # == check 5
  
  # which model is better?
  # Fit the model, the covariate goes first
  lm_model1 <- lm(EE ~  Bodyweight * Group, data = dataABC)
  result_1 = summary(lm_model1)
  
  lm_model2 <-  update( lm_model1, . ~ Bodyweight + Group) 
  result_2 = summary(lm_model2)
  
  if (  result_2$adj.r.squared > result_1$adj.r.squared ){
    print(result_1$adj.r.squared )
    print(result_2$adj.r.squared )
    print ("======== warning: The model without interaction may work better according to the R squared ======")
  }
  
  # == check 6
  
  if (   step2_anova[2,"p"] > 0.05){
    
    print ("======== warning: Body weight have no significant effect on food intake =====")
  }
  
  
  print ("======== End of assumption check  ======")
  # ========================== Computation ======================
  
  if (method == "III"){
    
    # This is slightly more involved than the type II results. 
    # First, it is necessary to set the contrasts option in R. 
    # Because the multi-way ANOVA model is over-parameterised, it is necessary to choose a contrasts setting that sums to zero, 
    # otherwise the ANOVA analysis will give incorrect results with respect to the expected hypothesis. (The default contrasts type does not satisfy this requirement.)
    
    options(contrasts = c("contr.sum","contr.poly"))
    
    # res.aov = rstatix::Anova(mod = ancova_model, type = "III")
    dataABC$Group = factor(dataABC$Group)
    ancova_model <- aov(EE ~ Group * Bodyweight, data = dataABC)
    res.aov = rstatix::Anova(mod = ancova_model, type = "III")
    
    res.aov <- dataABC %>% anova_test( EE ~  Bodyweight + Group + Bodyweight:Group, type = 3)
    
    # ===== if we found interaction effect, 
    # then we will use simple effect analysis to further see the difference in certain level of body weight ===
    covariate_level = quantile(dataABC$Bodyweight, probs = c(0.25,0.5,0.75))
    
    #covariate_level = quantile(dataABC$Bodyweight, probs = seq(.05, .95, by = .1))
    
    emm <- emmeans(ancova_model, ~ Group | Bodyweight, at = list(Bodyweight = covariate_level)) 
    # summary(emm) 
    
    # manipulate the result
    res.emmeans <- emm %>% tibble::as_tibble() %>% 
      dplyr::rename(se = "SE", conf.low = "lower.CL", conf.high = "upper.CL") %>% 
      mutate(method = "Emmeans test")
    
    # pairwise_comparisons <- pairwise_comparison_after_emmeans(emm)
    # by default, pairs() function use "tukey", we will use "bonferroni" in our project
    # res.compare <- summary(pairs(emm), adjust = "bonferroni") %>% as_tibble()
    
    comparisons <- emmeans::contrast(emm, method = "pairwise", adjust = "bonferroni")
    comparisons <- tidy(comparisons, conf.int = TRUE, conf.level = 0.95)
    res.comparisons <- comparisons %>% tidyr::separate(col = "contrast", 
                                                       into = c("group1", "group2"), sep = " - ") %>% 
      add_significance("adj.p.value")
    
    to.remove <- c("estimate", "std.error", "conf.low", "conf.high", "null.value")
    to.keep <- setdiff(colnames(res.comparisons), to.remove)
    res.comparisons  <- res.comparisons [, to.keep]
    
    # merge emmean 
    res.anno = res.emmeans %>% dplyr::select("Group","Bodyweight","emmean")
    
    res.comparisons <- res.comparisons %>% left_join(res.anno,by = c("group1" = "Group","Bodyweight")) %>% 
      dplyr::rename(emmean_group1 = emmean ) %>% 
      left_join(res.anno,by = c("group2" = "Group", "Bodyweight")) %>% 
      dplyr::rename(emmean_group2 = emmean )
    
    # === use Johnson-Neyman intervals ===
    # USE INTERACTIONS to find the interval
    # Convert the factor variable 'Group' into numeric dummy variables
    dummy_vars <- model.matrix(~ Group - 1, dataABC)
    # Add the new dummy variables to your data frame
    new_data <- cbind(dataABC, dummy_vars)
    
    points = dataABC %>% mutate(group = Group)
    if (is.factor(points[,"group"])){
      group_list = unique(points[,"group"]) %>% unlist()%>% as.vector()
    }else{
      group_list = unique(points[,"group"]) %>% unlist() %>% as.vector()
    }
    possible_combinations = sets::set_power(group_list) %>% as.vector() 
    compared_group = list()
    for (i in seq(length(possible_combinations))){
      result_c = possible_combinations[[i]] %>% as.character()
      if (length(result_c) == length(group_list)-2){
        sub_compared_group = group_list[!group_list %in% result_c]
        compared_group[[length(compared_group)+1]] =  sub_compared_group 
      }
    }
    
    # formula 
    
    stat.test = data.frame(group1 = character(), group2 = character(),pval_or_low = character(), psign_or_high = character(), interactions = character())
    for (x in compared_group){
      part_data = new_data %>% filter(Group %in% x)
      part_data$Group = factor(part_data$Group,levels = x)
      
      library(interactions)
      # Update the model using the new data frame with numeric dummy variables
      part_anova = part_data %>% anova_test(formula = EE ~ Group * Bodyweight,type = 2)
      print (x)
      # ===== check if the interaction also stay siginificant in two groups ====
      if (   part_anova[3,"p"] < 0.05 ){
        print ("======== warning: The covariate (body weight) and the treatment (Group) have interaction======")
        # rename the group name we want to 
        # we will use the first one as the referecen, e.g. c(A,B), A will be the reference
        target_name = paste0("Group",x[[2]])
        part_data <- part_data %>% rename(dum_group = `target_name`)
        
        ancova_model_updated <- lm(EE ~ Bodyweight * dum_group, data = part_data)
        # Perform the region of significance analysis using the new model
        region_of_significant <- interactions::sim_slopes(ancova_model_updated, 
                                                          pred = dum_group, 
                                                          modx = Bodyweight,    
                                                          johnson_neyman = TRUE,
                                                          control.fdr = TRUE,
                                                          jnalpha = 0.05)
        res.region = region_of_significant$jn[[1]]
        
        # JN interval
        res.region$bounds
        
        # JN plot   
        p_report3 = res.region$plot
        
        stat.test[nrow(stat.test) + 1,] = c(x[1],x[2],res.region$bounds[1],res.region$bounds[2],"interaction")  
        
        ggsave(p_report3, filename = paste0("../../HKFI-FoodIntake/results/",name,"/ABC_JN_interval.png"), dpi = 320)
        
        
      }else{
        
        print ("======== notification: The covariate (body weight) and the treatment (Group) have no interaction======")
        # We will use normal TukeyHSD to compare
        part_ancova_model = aov(EE ~ Bodyweight + Group + Bodyweight:Group,data = part_data)
        
        #posthoc = glht(part_ancova_model, linfct = mcp( Group = "Tukey"))
        #res.comparisons = tidy(summary(posthoc)) %>% tidyr::separate(col = "lhs", 
        #                                                             into = c("group1", "group2"), sep = " - ") %>% 
        #  add_significance("p.value") %>% select(-"rhs")
        
        res.comparisons_part <-  part_data %>% 
          emmeans_test(
            EE ~ Group, covariate = Bodyweight,
            p.adjust.method = "bonferroni"
          )
        
        res.emmeans_part =  get_emmeans(res.comparisons_part)
        
        
        stat.test[nrow(stat.test) + 1,] = c(x[1],x[2],res.comparisons_part$p.adj,res.comparisons_part$p.adj.signif,"no interaction")  
        
        
      }
    }
    
    write_tsv(stat.test,paste0("../../HKFI-FoodIntake/results/",name,"/ABC_pairwise_interaction_analysis_result.tsv"))
    
  }else if (method == "II"){
    stat.test= NULL
    # ==== identical with those in rstatix
    #dataABC$Group = factor(dataABC$Group)
    #ancova_model <- aov(data = dataABC, foodIntake ~ Bodyweight + Group) 
    #res.aov = rstatix::Anova(mod = ancova_model, type = "II")
    
    res.aov <- dataABC %>% anova_test( EE ~  Bodyweight * Group, type = 2)
    # posthoc = glht(ancova_model, linfct = mcp( Group = "Tukey"))
    # res.comparisons = tidy(summary(posthoc)) 
    
    # === we get the marginal mean by use the median bodyweight 
    # covariate_level = quantile(dataABC$Bodyweight, probs = c(0.5))
    
    # emm <- emmeans(ancova_model, ~ Group | Bodyweight, at = list(Bodyweight = covariate_level)) 
    
    # res.emmeans <- emm %>% tibble::as_tibble() %>% 
    #  dplyr::rename(se = "SE", conf.low = "lower.CL", conf.high = "upper.CL") %>% 
    #  mutate(method = "Emmeans test")
    
    res.comparisons <- dataABC %>% 
      emmeans_test(
        EE ~ Group, covariate = Bodyweight,
        p.adjust.method = "bonferroni"
      )
    
    res.emmeans =  get_emmeans(res.comparisons)
    
  }
  
  write_tsv(get_anova_table(res.aov),paste0("../../HKFI-FoodIntake/results/",name,"/ABC_ancova_result.tsv"))
  write_tsv(res.comparisons ,paste0("../../HKFI-FoodIntake/results/",name,"/ABC_pairwise_comparison_result.tsv"))
  write_tsv(res.emmeans ,paste0("../../HKFI-FoodIntake/results/",name,"/ABC_adjusted_mean_result.tsv"))
  
  print ("======== End of computation ======")
  
  # ========================== Final visualization ======================
  
  # Pairwise comparisons through emmeans
  # simple effect analysis at the level of mean body weight 
  library(emmeans)
  pwc <- dataABC %>% 
    emmeans_test(
      EE ~ Group, covariate = Bodyweight,
      p.adjust.method = "bonferroni"
    )
  
  # get_emmeans(pwc)
  pwc <- pwc %>% add_xy_position(x = "Group", fun = "mean_se")
  
  p_report2 = ggline(get_emmeans(pwc), x = "Group", y = "emmean") +
    geom_errorbar(aes(ymin = conf.low, ymax = conf.high), width = 0.2) + 
    stat_pvalue_manual(pwc, hide.ns = TRUE, tip.length = FALSE) +
    labs(
      subtitle = get_test_label(res.aov, detailed = TRUE, row = 2 ),
      caption = get_pwc_label(pwc),
      y = "Least-squares means"
    )
  
  
  ggsave(p_report2, filename = paste0("../../HKFI-FoodIntake/results/",name,"/ABC_Least_squares_means.png"), dpi = 320)
  
  
  
  plot1 = ggscatter(
    dataABC, x = "Bodyweight", y = "EE", 
    add = "reg.line",size = 3, color = "Group")+
    geom_point(aes(color = Group), size = 3) + 
    geom_point(shape = 1, color = "black", size = 3) + 
    stat_regline_equation(
      aes(label =  paste(..eq.label.., ..rr.label.., sep = "~~~~"), color = Group)
    )  +
    geom_rug()+ 
    theme_pubr() + 
    labs(y = "Energy expenditure", x = "Body weight (g)") + 
    ggsci::scale_color_lancet() +
    theme(legend.position = "bottom") 
  
  # if we have interaction effect, we will draw the JN interval 
  if (!is.null(stat.test) & length(unique(stat.test$interactions)) == 2 & table(stat.test$interactions)[1] == 1 ){
    upper_bound = stat.test %>% filter(interactions == "interaction") %>% pull("psign_or_high") %>% as.numeric()
    lower_bound = stat.test %>% filter(interactions == "interaction") %>% pull("pval_or_low") %>% as.numeric()
    min_bound  = min(dataABC$Bodyweight)
    max_bound  = max(dataABC$Bodyweight)
    
    y_min_food_intake = max(dataABC$EE) + 0.03* max(dataABC$EE)
    y_max_food_intake = max(dataABC$EE) + 0.02* max(dataABC$EE)
    plot1 <- plot1 +
      annotate(
        xmin =  c(-Inf,upper_bound),     # plus one 
        xmax = c(upper_bound,Inf),
        ymin =  y_min_food_intake,
        ymax =  y_max_food_intake,
        geom = "rect",
        fill = c("red", "blue")
      )
  }
  
  
  dens1 <- ggplot(dataABC, aes(x = Bodyweight, fill = Group)) + 
    geom_density(alpha = 0.4) + 
    theme_void() + 
    theme(legend.position = "none")+ 
    ggsci::scale_fill_lancet()
  
  dens2 <- ggplot(dataABC, aes(x = EE, fill = Group)) + 
    geom_density(alpha = 0.4) + 
    theme_void() + 
    theme(legend.position = "none") + 
    coord_flip()+ 
    ggsci::scale_fill_lancet()
  
  
  p_report1  = dens1 + plot_spacer() + plot1 + dens2 + 
    plot_layout(ncol = 2, nrow = 2, widths = c(4, 1), heights = c(1, 4))
  
  ggsave(p_report1 , filename = paste0("../../HKFI-FoodIntake/results/",name,"/ABC_scatter_with_density.png"), dpi = 320)
  
}



source_1 <- function(){
  

  
  require(reghelper)
  # ancova_model
  model = lm(foodIntake ~ Bodyweight * Group, data = dataABC)
  
  simple_slopes(model)
  
  as.data.frame(state.x77)
  
  
  # === use linear model == 
  points = dataABC %>% mutate(group = Group)
  if (is.factor(points[,"group"])){
    group_list = unique(points[,"group"]) %>% unlist()%>% as.vector()
  }else{
    group_list = unique(points[,"group"]) %>% unlist() %>% as.vector()
  }
  possible_combinations = sets::set_power(group_list) %>% as.vector() 
  compared_group = list()
  for (i in seq(length(possible_combinations))){
    result_c = possible_combinations[[i]] %>% as.character()
    if (length(result_c) == length(group_list)-2){
      sub_compared_group = group_list[!group_list %in% result_c]
      compared_group[[length(compared_group)+1]] =  sub_compared_group 
    }
  }
  
  print(compared_group)  
  
  stat.test = data.frame(group1 = character(), group2 = character(), pval = character(), direction = character())
  for (x in compared_group){
    model = lm(data = dataABC %>% dplyr::filter(Group %in% x), formula = foodIntake ~ Group + Bodyweight)
    result = summary(model)
    stat.test[nrow(stat.test) + 1,] = c(x[1],x[2],result$coefficients[2,4],sign(result$coefficients[2,1]))                                    
  }
  
  
}


