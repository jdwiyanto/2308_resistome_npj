library(tidyverse)
library(phyloseq)
library(readxl)
library(broom)
library(purrr)
library(vegan)
library(Maaslin2)
library(corrplot)
library(ggraph)
library(igraph)
library(tidygraph)
library(lme4)

load('env_rev1.Rdata')

# import metaphlan output
meta = read.delim(file = 'input/meta_merged.tsv', sep = '\t')

# clean up names
namelist = 
  meta[1,-1] %>% 
  as.character() %>% 
  str_extract(., pattern = 'AMD_?[0-9][0-9][0-9]') %>% 
  gsub('_', '', .) %>%
  gsub('AMD', '', .) %>%
  gsub('^0+', '', .) %>%
  paste0('JD-', .)

# import metadata
df = read_xlsx(path = 'input/qs_data_final_deidentified.xlsx',
               sheet = 'Sheet1') %>% filter(id %in% namelist) %>% 
  column_to_rownames('id')

# filter metaphlan output to include only those in metadata
colnames(meta) = c('ID', namelist)
meta = 
  meta %>% 
  slice(-1) %>% 
  filter(grepl('s__', ID)) %>% 
  mutate(across(!ID, .fns = as.numeric)) %>%
  select(ID, rownames(df)) 

# generate taxa table based on metaphlan output
tax = meta$ID %>% str_split_fixed(pattern = '\\|', n=7) %>% data.frame
colnames(tax) = c('kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species')
rownames(tax) = tax$species
tax2 = tax %>% as.matrix

# revise meta rownames as species name 
rownames(meta) = tax$species

# revise meta to remove ID column
meta = meta %>% select(-ID) 

# import args-oap output
arg_gene = read_table(file = 'input/args_soap_out/normalized_cell.gene.txt') %>% column_to_rownames('gene') 
arg_sub = read_table(file = 'input/args_soap_out/normalized_cell.subtype_nospace.txt') %>% column_to_rownames('subtype')
arg_type = read_table(file = 'input/args_soap_out/normalized_cell.type.txt') %>% column_to_rownames('type') 

colnames(arg_gene) =
  colnames(arg_gene) %>% 
  gsub('_', '', .) %>%
  gsub('AMD', '', .) %>%
  gsub('^0+', '', .) %>%
  paste0('JD-', .)

colnames(arg_sub) =
  colnames(arg_sub) %>% 
  gsub('_', '', .) %>%
  gsub('AMD', '', .) %>%
  gsub('^0+', '', .) %>%
  paste0('JD-', .)

colnames(arg_type) =
  colnames(arg_type) %>% 
  gsub('_', '', .) %>%
  gsub('AMD', '', .) %>%
  gsub('^0+', '', .) %>%
  paste0('JD-', .)

# only keep arg samples in metadata
arg_gene = arg_gene %>% select(rownames(df)) 
arg_sub= arg_sub %>% select(rownames(df)) 
arg_type = arg_type %>% select(rownames(df)) 

# create phyloseq objects
# metaphlan abundance
ps = phyloseq(sample_data(df),
              otu_table(meta, taxa_are_rows = T),
              tax_table(tax2))

# resistome gene
arg_gene_tax = tibble(gene = rownames(arg_gene)) %>% as.matrix()
rownames(arg_gene_tax) = rownames(arg_gene)

ps.arg.gene = 
  phyloseq(sample_data(df),
         otu_table(arg_gene, taxa_are_rows = T),
         tax_table(arg_gene_tax))

# resistome type
arg_type_tax = tibble(gene = rownames(arg_type)) %>% as.matrix()
rownames(arg_type_tax) = rownames(arg_type)

ps.arg.type = 
  phyloseq(sample_data(df),
           otu_table(arg_type, taxa_are_rows = T),
           tax_table(arg_type_tax))

# resistome sub - problematic -- substitute single quote with dot to solve issue
arg_sub_newtaxname = gsub(pattern = "'", replacement = "\\.", rownames(arg_sub))
rownames(arg_sub) = arg_sub_newtaxname

arg_sub_tax = tibble(gene = rownames(arg_sub)) %>% separate(col = 'gene', into = c('type', 'sub'), sep = '__') %>% as.matrix()
rownames(arg_sub_tax) = rownames(arg_sub)

arg_sub_otu = arg_sub %>% otu_table(taxa_are_rows = T) 

ps.arg.sub = 
  phyloseq(sample_data(df),
           arg_sub_otu,
           tax_table(arg_sub_tax))

# manuscript section - resistome profile of the segamat cohort
# no of resistance genes identified

arg_sub %>% rowSums %>% mean
arg_sub %>% rowSums %>% sd

# filter off <0.5% abundance and <10% prevalence of gene
ps.arg.sub %>% transform_sample_counts(fun = function(x) { x / sum(x) }) %>% microbiome::core(detection = 0.5/100, prevalence = 10/100) 

arg.filt.names = ps.arg.sub %>% transform_sample_counts(fun = function(x) { x / sum(x) }) %>% microbiome::core(detection = 0.5/100, prevalence = 10/100)  %>% taxa_names()

# how many antibiotic resistance genes types do the filtered prevalence confer
ps.arg.sub %>% prune_taxa(arg.filt.names, .) %>% tax_glom('type')

# figure1a - ARG abundance

fig.arg.sum = 
  ps.arg.sub %>% otu_table() %>% data.frame(check.names = F) %>% rownames_to_column('id') %>% pivot_longer(!id) %>%
  rename(gene=id, id=name) %>%
  inner_join(df.clean %>% select(id,ethnicity), by='id') %>%
  mutate(ethnicity = ifelse(ethnicity=='Native', 'Jakun', ethnicity)) %>%
  ggplot(aes(x=reorder(id,value), y=value, color=ethnicity)) +
  geom_point() +
  ggsci::scale_color_d3(name='') +
  theme_bw() +
  labs(y='ARG copy per cell', x='Subject') +
  theme(axis.text.x = element_blank())

# ARG type abundance description
arg.top20 = 
  ps.arg.sub %>% prune_taxa(arg.filt.names, .) %>% 
  otu_table() %>% 
  data.frame %>%
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>%
  separate(col = 'rowname', into = c('family', 'gene'), sep = '__') %>%
  group_by(gene) %>%
  summarise(mean = mean(value)) %>%
  slice_max(mean, n=20) %>% 
  pull(gene)

# mean and sd abundance of each family 
ps.arg.sub %>% prune_taxa(arg.filt.names, .) %>% tax_glom('type') %>%
  otu_table %>%
  data.frame  %>% 
  rownames_to_column() %>%
  pivot_longer(!rowname) %>%
  separate(col = 'rowname', into = c('family', 'gene'), sep = '__') %>%
  group_by(family) %>%
  summarise(mean = mean(value),
            sd = sd(value)) %>%
  arrange(-mean)

# mean and sd abundance of each gene
ps.arg.sub %>% prune_taxa(arg.filt.names, .) %>% 
  otu_table %>%
  data.frame  %>% 
  rownames_to_column() %>%
  pivot_longer(!rowname) %>%
  separate(col = 'rowname', into = c('family', 'gene'), sep = '__') %>%
  group_by(gene) %>%
  summarise(mean = mean(value),
            sd = sd(value)) %>%
  arrange(-mean)

# mean and sd abundance of ecoli genes
ps.arg.sub %>% prune_taxa(arg.filt.names, .) %>% 
  otu_table %>%
  data.frame  %>% 
  rownames_to_column() %>%
  pivot_longer(!rowname) %>%
  separate(col = 'rowname', into = c('family', 'gene'), sep = '__') %>%
  group_by(gene) %>%
  summarise(mean = mean(value),
            sd = sd(value)) %>%
  arrange(-mean) %>%
  filter(grepl('Escherichia', gene))

# visualise top 20 gene subtype
fig.arg.top20 = 
  ps.arg.sub %>% prune_taxa(arg.filt.names, .) %>% 
  otu_table() %>% 
  data.frame %>%
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>%
  separate(col = 'rowname', into = c('family', 'gene'), sep = '__') %>%
  filter(gene %in% arg.top20) %>%
  mutate(gene = ifelse(grepl('Bifidobacteria', gene), 'ileS', gene)) %>%
  mutate(family = ifelse(family=='macrolide-lincosamide-streptogramin', 'MLS', family)) %>%
  ggplot(aes(x=value, y=reorder(gene,value), fill=family)) + 
  geom_jitter(width = 0.1, size= 0.05, alpha=0.8, color='grey') + 
  geom_boxplot(outlier.shape = NA) +
  ggsci::scale_fill_d3(name='') + 
  scale_x_continuous(limits = c(0,0.6)) +
  theme_bw() + 
  labs(x='Gene copy per cell', y='Antibiotic resistance gene') 

# ARG sub abundance description
arg.filt.sub.names = ps.arg.sub %>% transform_sample_counts(fun = function(x) { x / sum(x) }) %>% microbiome::core(detection = 0.5/100, prevalence = 10/100) %>% taxa_names()

ps.arg.sub %>% prune_taxa(arg.filt.sub.names, .) %>% 
  otu_table() %>% 
  data.frame %>%
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>%
  # mutate(value = log2(value+1/1000000)) %>%
  ggplot(aes(x=value, y=reorder(rowname,value))) + geom_boxplot()
  
# how many  ecoli genes are here
ps.arg.sub %>% prune_taxa(arg.filt.sub.names, .) %>% 
  otu_table() %>% 
  data.frame %>%
  rownames_to_column() %>% 
  pivot_longer(!rowname) %>%
  filter(grepl('Escherichia', rowname)) %>%
  ggplot(aes(x=value, y=reorder(rowname,value))) + geom_boxplot()

# import cleaned metadata
df.clean = read_table(file = 'input/suppl_table3.tsv')
df.clean = df %>% rownames_to_column('id') %>% select(id, household) %>% inner_join(df.clean, by='id') %>% tibble

# import demographic metadata only
df.demo = read_table(file='input/suppl_table1.tsv')
df.demo = df.demo %>% inner_join(df.clean %>% select(id, income), by='id')

# resistome diversity
rich.arg = ps.arg.sub %>% estimate_richness(measures='Shannon') %>% rownames_to_column('id')

# resistome diversity across demo var
jd_bmi_cat_asian = 
  function (df, bmivar) {
    df[[bmivar]] = as.numeric(df[[bmivar]])
    df$bmi_cat <- ifelse(df[[bmivar]] < 18.5, "underweight", 
                         ifelse(between(df[[bmivar]], 18.5, 22.99), "normal", 
                                ifelse(between(df[[bmivar]], 23, 24.99), "overweight",
                                       ifelse(between(df[[bmivar]], 25, 29.99), "obese1", "obese2"))))
    df$bmi_cat = as.factor(df$bmi_cat)
    df$bmi_cat = factor(df$bmi_cat, levels = c("underweight", "normal", "overweight", "obese1", "obese2"))
    return(df$bmi_cat)
  }

jd_age_dec = function (df, agevar) {
  df[[agevar]] <- as.numeric(df[[agevar]])
  df$age_cat <- ifelse(df[[agevar]] < 11, "<11", ifelse(
    between(df[[agevar]], 11, 20), "11-20", ifelse(
      between(df[[agevar]], 21, 30), "21-30", ifelse(
        between(df[[agevar]], 31, 40), "31-40", ifelse(
          between(df[[agevar]], 41, 50), "41-50", ifelse(
            between(df[[agevar]], 51, 60), "51-60", ifelse(
              between(df[[agevar]], 61, 70), "61-70", ifelse(
                between(df[[agevar]], 71, 80), "71-80", ifelse(
                  between(df[[agevar]], 81, 90), "81-90", ">90")))))))))
  
  df$age_cat = as.factor(df$age_cat)
  df$age_cat <- factor(df$age_cat, levels = c("<11", "11-20", "21-30", "31-40", "41-50", "51-60", "61-70", "71-80", "81-90", ">90"))
  return(df$age_cat)
}

df.demo.res = 
  df.demo %>% inner_join(rich.arg, by='id') %>% 
  mutate(income = as.factor(income), bmi = as.numeric(bmi), household = as.factor(household)) %>%
  mutate(bmi_cat = jd_bmi_cat_asian(df = ., bmivar = 'bmi'), age_dec = jd_age_dec(df = ., agevar = 'age'))

# transform age data using boxcox, since better than log2 and square root in normalising age dist although not perfect
boxcox.age = MASS::boxcox(lm(df.demo.res$age ~ 1))
boxcox.age.lambda = boxcox.age$x[which.max(boxcox.age$y)]
((df.demo.res$age ^ boxcox.age.lambda) / boxcox.age.lambda) %>% shapiro.test()

df.demo.res$age %>% sqrt %>% shapiro.test()
df.demo.res$age %>% log2 %>% shapiro.test()

df.demo.res$age_transformed = ((df.demo.res$age ^ boxcox.age.lambda) / boxcox.age.lambda)

# transform bmi using log2
df.demo.res$bmi_transformed = df.demo.res$bmi %>% log2()


# linear model for continuous var
df.demo.res %>% select(Shannon, age_transformed, bmi_transformed) %>%
  pivot_longer(!Shannon) %>%
  mutate(Shannon = log2(Shannon)) %>%
  group_by(name) %>%
  nest(!name) %>%
  mutate(test = map(.x = data,
                    .f = ~lm(Shannon ~ value, data  = .x) %>% tidy)) %>%
  ungroup() %>%
  unnest(test) %>%
  filter(term != '(Intercept)') %>%
  select(-data)

# kruskal.test for categorical var
kw.arg.cat = 
  df.demo.res %>% select(-id,-household, -bmi, -bmi_transformed, -age_transformed) %>%
  mutate(age_quant = gtools::quantcut(age)) %>%
  select(-age) %>% 
  pivot_longer(!Shannon) %>%
  group_by(name) %>%
  nest(!name) %>%
  mutate(test = purrr::map(.x = data,
                           .f = ~kruskal.test(Shannon ~ value, data = .x) %>% tidy)) %>%
  ungroup() %>%
  unnest(test) %>%
  select(-data)
  
  
# posthoc mann whitney for age dec
df.demo.res %>% 
  mutate(age_quant = gtools::quantcut(age)) %>%
  select(-id,-household, -bmi, -age, -bmi_transformed, -age_transformed) %>%
  pivot_longer(!Shannon) %>%
  filter(name == 'age_quant') %>%
  group_by(name) %>%
  nest(!name) %>%
  mutate(test = map(.x = data,
                    .f = ~pairwise.wilcox.test(x = .x$Shannon,
                                               g = .x$value,
                                               p.adjust.method = 'BH') %>% tidy)) %>%
  ungroup() %>%
  unnest(test) %>%
  select(-data)
df.demo.res$Shannon %>% hist

# diversity test for all lifestyle parameter, continuous var
lm.cont.out = 
  df.clean %>% select(id, age, bmi, exercise, ends_with('frequency'), starts_with('house_')) %>% inner_join(rich.arg, by='id') %>%
  pivot_longer(!c(id, Shannon)) %>%
  mutate(value = log2(value+(1/1000000))) %>%
  group_by(name) %>%
  nest(!name) %>%
  mutate(test = map(.x = data,
                    .f = ~lm(Shannon ~ value, data  = .x) %>% tidy)) %>%
  ungroup() %>%
  unnest(test) %>%
  filter(term != '(Intercept)') %>%
  select(-data)

lm.cont.out %>% filter(p.value<0.1)


# diversity test for all lifestyle parameter, categorical var
kw.cat.out =
  df.clean %>% 
  mutate(age_quant = gtools::quantcut(age)) %>%
  select(!ends_with('frequency')) %>% 
  select(!starts_with('house_')) %>% 
  select(-exercise) %>% 
  mutate(income = as.factor(income)) %>%
  select(-batch, -household, -bmi, -age) %>%
  inner_join(rich.arg, by='id') %>%
  pivot_longer(!c(id, Shannon)) %>% 
  drop_na() %>%
  filter(value != 'N/A') %>% 
  mutate(value = str_to_sentence(value)) %>%
  group_by(name) %>%
  nest(!name) %>%
  mutate(test = purrr::map(.x = data,
                           .f = ~kruskal.test(Shannon ~ value, data = .x) %>% tidy)) %>%
  ungroup() %>%
  unnest(test) %>%
  select(-data) 

fig.arg.cat = 
  df.clean %>% 
  mutate(age_quant = gtools::quantcut(age)) %>%
  inner_join(rich.arg, by='id') %>% 
  select(id, Shannon, kw.cat.out %>% filter(p.value < 0.1) %>% arrange(p.value) %>% pull(name)) %>%
  pivot_longer(!c(id,Shannon)) %>%
  drop_na() %>%
  filter(value != 'N/A') %>%
  mutate(value = str_to_sentence(value)) %>%
  mutate(name = ifelse(name=='age_quant', 'Age (Quantiles)', ifelse(name=='agecat', 'Age (Decade)', ifelse(name=='animal_chicken', 'Chicken near house', ifelse(name=='chronic_diabetes', 'Diabetic', ifelse(name=='drinker', 'Alcohol drinking habit', ifelse(name=='esbl', 'Stool ESBL culture', ifelse(name=='med_suppl', 'On active medication', ifelse(name=='mould', 'Stool mould culture', ifelse(name=='toilet_location_outside', 'Toilet outside main house structure', ifelse(name=='toilet_nature_sit', 'Seating toilet type', ifelse(name=='water_processing_no', 'Water not processed', ifelse(name=='yeast', 'Stool yeast culture', 'notdefined'))))))))))))) %>%
  ggplot(aes(x=value, y=Shannon)) +
  geom_jitter(width=0.1, size=0.2) +
  geom_boxplot(alpha=0.1) + 
  facet_wrap(~name, scales='free_x') +
  theme_bw()
  
# beta diversity

# ordination
ps.arg.sub %>% 
  microbiome::transform('clr') %>%
  plot_ordination(physeq = ., ordination = ordinate(physeq = ., method = 'PCoA', distance = 'euclidean'), color='ethnicity' ) +
  theme_bw() +
  ggsci::scale_color_d3(name='') +
  theme(legend.position = 'top') +
  stat_ellipse()

########################################################
### PERMANOVA OF DEMOGRAPHIC AND LIFESTYLE VARIABLES ###
########################################################

ps.arg.sub.cleandf = phyloseq(ps.arg.sub %>% otu_table,
                              ps.arg.sub %>% tax_table,
                              df.clean %>%
                                mutate(across(.fns = ~ifelse(.=='N/A', NA, .)),
                                       gender = str_to_sentence(gender),
                                       household = as.factor(household),
                                       batch = as.factor(batch)) %>% 
                                column_to_rownames('id') %>% sample_data)


df.arg.cont = 
  df.clean %>%
  mutate(across(.fns = ~ifelse(.=='N/A', NA, .)),
         gender = str_to_sentence(gender),
         household = as.factor(household),
         batch = as.factor(batch),
         income = as.factor(income)) %>%
  select(id, household, batch, age, exercise, bmi, starts_with('house_'), ends_with('frequency')) %>%
  mutate(across(.fns = ~ifelse(is.na(.), median(., na.rm = T), .))) %>%
  column_to_rownames('id')

ps.arg.sub.cleandf.adoniscont = 
  phyloseq(sample_data(df.arg.cont),
           otu_table(ps.arg.sub),
           tax_table(ps.arg.sub))

adonis.cont.out = 
  lapply(names(df.arg.cont)[-c(1:2)], function(var) {
    
    otu = ps.arg.sub.cleandf.adoniscont %>% microbiome::transform('clr') %>% otu_table() %>% t
    form = paste0('otu ~ batch + household + ', var) %>% as.formula()
    out = adonis2(formula = form, 
                  data = df.arg.cont,
                  method = 'euclidean')
    out2 = tibble(var = var,
                  df = out$Df[3],
                  ss = out$SumOfSqs[3],
                  r2 = out$R2[3],
                  f = out$F[3],
                  p = out$`Pr(>F)`[3])
    return(out2)
    
  }) %>% bind_rows()

adonis.cont.out %>% filter(p<0.1)

# categorical variable
Mode <- function(x) {
  uniqx <- unique(x)
  uniqx[which.max(tabulate(match(x, uniqx)))]
}

df.arg.cat = 
  df.clean %>%
  mutate(across(.fns = ~ifelse(.=='N/A', NA, .)),
         gender = str_to_sentence(gender),
         household = as.factor(household),
         batch = as.factor(batch),
         income = as.factor(income),
         calbicans = ifelse(is.na(calbicans), 'Negative', calbicans),
         esbl = ifelse(is.na(esbl), 'Negative', esbl),
         yeast = ifelse(is.na(yeast), 'Negative', yeast)) %>%
  select(-age, -exercise, -bmi, -starts_with('house_'), -ends_with('frequency')) %>% 
  mutate(across(.fns = ~ifelse(is.na(.), Mode(.), .)),
         income = as.factor(income)) %>% column_to_rownames('id')

ps.arg.sub.cleandf.adoniscat = 
  phyloseq(sample_data(df.arg.cat),
         otu_table(ps.arg.sub),
         tax_table(ps.arg.sub))

adonis.cat.out = 
  lapply(names(df.arg.cat)[-c(1:2)], function(var) {
  
  otu = ps.arg.sub.cleandf.adoniscat %>% microbiome::transform('clr') %>% otu_table() %>% t
  form = paste0('otu ~ batch + household + ', var) %>% as.formula()
  out = adonis2(formula = form, 
                data = df.arg.cat,
                method = 'euclidean')
  out2 = tibble(var = var,
                df = out$Df[3],
                ss = out$SumOfSqs[3],
                r2 = out$R2[3],
                f = out$F[3],
                p = out$`Pr(>F)`[3])
  return(out2)
  
}) %>% bind_rows()

adonis.cat.out %>% filter(p<0.1) 

# bind cont and cat
adonis.out.all = bind_rows(adonis.cont.out, adonis.cat.out) %>% filter(p<0.05) %>% arrange(-r2)
# exclude several var to remove collinearity
adonis.out.all.filt = adonis.out.all %>% filter(!(var %in% c('Location', 'calbicans', 'toilet_location_inside', 'toilet_type_flush_sewer')))

###################################
### PERMANOVA SUBGROUP ANALYSIS ###
###################################

# Univariate PERMANOVA subgroup analysis - age group 
df.age_quant = df.clean %>% select(id, age) %>% mutate(age_quant = gtools::quantcut(age)) %>% select(-age)

adonis.uni.cont.age =
  lapply(df.age_quant$age_quant %>% levels, function(var) {
  
  tokeep = df.age_quant %>% filter(age_quant == var) %>% pull(id)
  loopdf = df.arg.cont[tokeep,]
  loopps = prune_samples(tokeep, ps.arg.sub.cleandf.adoniscont)
  
  loop.cont.out = 
    lapply(names(df.arg.cont)[-c(1:2)], function(var2) {
      
      otu = loopps %>% microbiome::transform('clr') %>% otu_table() %>% t
      form = paste0('otu ~ batch + household + ', var2) %>% as.formula()
      out = adonis2(formula = form, 
                    data = loopdf,
                    method = 'euclidean')
      out2 = tibble(subgroup = var,
                    var = var2,
                    df = out$Df[3],
                    ss = out$SumOfSqs[3],
                    r2 = out$R2[3],
                    f = out$F[3],
                    p = out$`Pr(>F)`[3])
      return(out2)
      
    }) %>% bind_rows() 
  
  return(loop.cont.out)
  
}) %>% bind_rows

adonis.uni.cont.age %>% filter(p<0.1)

adonis.uni.cat.age =
  lapply(df.age_quant$age_quant %>% levels, function(var) {
    
    tokeep = df.age_quant %>% filter(age_quant == var) %>% pull(id)
    loopdf = df.arg.cat[tokeep,]
    loopps = prune_samples(tokeep, ps.arg.sub.cleandf.adoniscat)
    
    loop.cat.out = 
      lapply(names(df.arg.cat)[-c(1:2)], function(var2) {
        
        if ((loopdf[[var2]] %>% unique %>% length) > 1 ) {
        
        otu = loopps %>% microbiome::transform('clr') %>% otu_table() %>% t
        form = paste0('otu ~ batch + household + ', var2) %>% as.formula()
        out = adonis2(formula = form, 
                      data = loopdf,
                      method = 'euclidean')
        out2 = tibble(subgroup = var,
                      var = var2,
                      df = out$Df[3],
                      ss = out$SumOfSqs[3],
                      r2 = out$R2[3],
                      f = out$F[3],
                      p = out$`Pr(>F)`[3])
        return(out2)
        
        } 
        
      }) %>% bind_rows() 
    
    return(loop.cat.out)
    
  }) %>% bind_rows

bind_rows(adonis.uni.cat.age, adonis.uni.cont.age) %>% filter(p<0.1) %>% 
  ggplot(aes(x=r2, y=var)) + geom_col() + facet_wrap(~subgroup, nrow=1) +
  theme_bw() +
  labs(x='PERMANOVA R2', y='Lifestyle variable')

# univariate PERMANOVA analysis - ethnicity
adonis.uni.cont.eth =
  lapply(df.arg.cat$ethnicity %>% unique, function(var) {
    
    tokeep = df.arg.cat %>% filter(ethnicity == var) %>% rownames()
    loopdf = df.arg.cont[tokeep,]
    loopps = prune_samples(tokeep, ps.arg.sub.cleandf.adoniscont)
    
    loop.cont.out = 
      lapply(names(df.arg.cont)[-c(1:2)], function(var2) {
        
        otu = loopps %>% microbiome::transform('clr') %>% otu_table() %>% t
        form = paste0('otu ~ batch + household + ', var2) %>% as.formula()
        out = adonis2(formula = form, 
                      data = loopdf,
                      method = 'euclidean')
        out2 = tibble(subgroup = var,
                      var = var2,
                      df = out$Df[3],
                      ss = out$SumOfSqs[3],
                      r2 = out$R2[3],
                      f = out$F[3],
                      p = out$`Pr(>F)`[3])
        return(out2)
        
      }) %>% bind_rows() 
    
    return(loop.cont.out)
    
  }) %>% bind_rows

adonis.uni.cont.eth %>% filter(p<0.1)

adonis.uni.cat.eth =
  lapply(df.arg.cat$ethnicity %>% unique, function(var) {
    
    tokeep = df.arg.cat %>% filter(ethnicity == var) %>% rownames()
    loopdf = df.arg.cat[tokeep,]
    loopps = prune_samples(tokeep, ps.arg.sub.cleandf.adoniscat)
    
    loop.cat.out = 
      lapply(names(df.arg.cat)[-c(1:2)], function(var2) {
        
        if ((loopdf[[var2]] %>% unique %>% length) > 1 ) {
          
          otu = loopps %>% microbiome::transform('clr') %>% otu_table() %>% t
          form = paste0('otu ~ batch + household + ', var2) %>% as.formula()
          out = adonis2(formula = form, 
                        data = loopdf,
                        method = 'euclidean')
          out2 = tibble(subgroup = var,
                        var = var2,
                        df = out$Df[3],
                        ss = out$SumOfSqs[3],
                        r2 = out$R2[3],
                        f = out$F[3],
                        p = out$`Pr(>F)`[3])
          return(out2)
          
        } 
        
      }) %>% bind_rows() 
    
    return(loop.cat.out)
    
  }) %>% bind_rows

bind_rows(adonis.uni.cat.eth, adonis.uni.cont.eth) %>% filter(p<0.1) %>% 
  ggplot(aes(x=r2, y=var)) + geom_col() + facet_wrap(~subgroup, nrow=1) +
  theme_bw() +
  labs(x='PERMANOVA R2', y='Lifestyle variable')

# subgroup permanova analysis - income
adonis.uni.cont.inc =
  lapply(df.arg.cat$income %>% levels, function(var) {
    
    tokeep = df.arg.cat %>% filter(income == var) %>% rownames()
    loopdf = df.arg.cont[tokeep,]
    loopps = prune_samples(tokeep, ps.arg.sub.cleandf.adoniscont)
    
    loop.cont.out = 
      lapply(names(df.arg.cont)[-c(1:2)], function(var2) {
        
        otu = loopps %>% microbiome::transform('clr') %>% otu_table() %>% t
        form = paste0('otu ~ batch + household + ', var2) %>% as.formula()
        out = adonis2(formula = form, 
                      data = loopdf,
                      method = 'euclidean')
        out2 = tibble(subgroup = var,
                      var = var2,
                      df = out$Df[3],
                      ss = out$SumOfSqs[3],
                      r2 = out$R2[3],
                      f = out$F[3],
                      p = out$`Pr(>F)`[3])
        return(out2)
        
      }) %>% bind_rows() 
    
    return(loop.cont.out)
    
  }) %>% bind_rows

adonis.uni.cont.inc %>% filter(p<0.1)

adonis.uni.cat.inc =
  lapply(df.arg.cat$income %>% levels, function(var) {
    
    tokeep = df.arg.cat %>% filter(income == var) %>% rownames()
    loopdf = df.arg.cat[tokeep,]
    loopps = prune_samples(tokeep, ps.arg.sub.cleandf.adoniscat)
    
    loop.cat.out = 
      lapply(names(df.arg.cat)[-c(1:2)], function(var2) {
        
        if ((loopdf[[var2]] %>% unique %>% length) > 1 ) {
          
          otu = loopps %>% microbiome::transform('clr') %>% otu_table() %>% t
          form = paste0('otu ~ batch + household + ', var2) %>% as.formula()
          out = adonis2(formula = form, 
                        data = loopdf,
                        method = 'euclidean')
          out2 = tibble(subgroup = var,
                        var = var2,
                        df = out$Df[3],
                        ss = out$SumOfSqs[3],
                        r2 = out$R2[3],
                        f = out$F[3],
                        p = out$`Pr(>F)`[3])
          return(out2)
          
        } 
        
      }) %>% bind_rows() 
    
    return(loop.cat.out)
    
  }) %>% bind_rows

adonis.uni.cat.inc %>% filter(p<0.1)

# subgroup permanova analysis - bmi
adonis.uni.cont.bmi =
  lapply(df.arg.cat$bmicat %>% unique, function(var) {
    
    tokeep = df.arg.cat %>% filter(bmicat == var) %>% rownames()
    loopdf = df.arg.cont[tokeep,]
    loopps = prune_samples(tokeep, ps.arg.sub.cleandf.adoniscont)
    
    loop.cont.out = 
      lapply(names(df.arg.cont)[-c(1:2)], function(var2) {
        
        otu = loopps %>% microbiome::transform('clr') %>% otu_table() %>% t
        form = paste0('otu ~ batch + household + ', var2) %>% as.formula()
        out = adonis2(formula = form, 
                      data = loopdf,
                      method = 'euclidean')
        out2 = tibble(subgroup = var,
                      var = var2,
                      df = out$Df[3],
                      ss = out$SumOfSqs[3],
                      r2 = out$R2[3],
                      f = out$F[3],
                      p = out$`Pr(>F)`[3])
        return(out2)
        
      }) %>% bind_rows() 
    
    return(loop.cont.out)
    
  }) %>% bind_rows

adonis.uni.cont.bmi %>% filter(p<0.1)

adonis.uni.cat.bmi =
  lapply(df.arg.cat$bmicat %>% unique, function(var) {
    
    tokeep = df.arg.cat %>% filter(bmicat == var) %>% rownames()
    loopdf = df.arg.cat[tokeep,]
    loopps = prune_samples(tokeep, ps.arg.sub.cleandf.adoniscat)
    
    loop.cat.out = 
      lapply(names(df.arg.cat)[-c(1:2)], function(var2) {
        
        if ((loopdf[[var2]] %>% unique %>% length) > 1 ) {
          
          otu = loopps %>% microbiome::transform('clr') %>% otu_table() %>% t
          form = paste0('otu ~ batch + household + ', var2) %>% as.formula()
          out = adonis2(formula = form, 
                        data = loopdf,
                        method = 'euclidean')
          out2 = tibble(subgroup = var,
                        var = var2,
                        df = out$Df[3],
                        ss = out$SumOfSqs[3],
                        r2 = out$R2[3],
                        f = out$F[3],
                        p = out$`Pr(>F)`[3])
          return(out2)
          
        } 
        
      }) %>% bind_rows() 
    
    return(loop.cat.out)
    
  }) %>% bind_rows

# subgroup permanova analysis - sex

adonis.uni.cont.sex =
  lapply(df.arg.cat$gender %>% unique, function(var) {
    
    tokeep = df.arg.cat %>% filter(gender == var) %>% rownames()
    loopdf = df.arg.cont[tokeep,]
    loopps = prune_samples(tokeep, ps.arg.sub.cleandf.adoniscont)
    
    loop.cont.out = 
      lapply(names(df.arg.cont)[-c(1:2)], function(var2) {
        
        otu = loopps %>% microbiome::transform('clr') %>% otu_table() %>% t
        form = paste0('otu ~ batch + household + ', var2) %>% as.formula()
        out = adonis2(formula = form, 
                      data = loopdf,
                      method = 'euclidean')
        out2 = tibble(subgroup = var,
                      var = var2,
                      df = out$Df[3],
                      ss = out$SumOfSqs[3],
                      r2 = out$R2[3],
                      f = out$F[3],
                      p = out$`Pr(>F)`[3])
        return(out2)
        
      }) %>% bind_rows() 
    
    return(loop.cont.out)
    
  }) %>% bind_rows

adonis.uni.cat.sex =
  lapply(df.arg.cat$gender %>% unique, function(var) {
    
    tokeep = df.arg.cat %>% filter(gender == var) %>% rownames()
    loopdf = df.arg.cat[tokeep,]
    loopps = prune_samples(tokeep, ps.arg.sub.cleandf.adoniscat)
    
    loop.cat.out = 
      lapply(names(df.arg.cat)[-c(1:2)], function(var2) {
        
        if ((loopdf[[var2]] %>% unique %>% length) > 1 ) {
          
          otu = loopps %>% microbiome::transform('clr') %>% otu_table() %>% t
          form = paste0('otu ~ batch + household + ', var2) %>% as.formula()
          out = adonis2(formula = form, 
                        data = loopdf,
                        method = 'euclidean')
          out2 = tibble(subgroup = var,
                        var = var2,
                        df = out$Df[3],
                        ss = out$SumOfSqs[3],
                        r2 = out$R2[3],
                        f = out$F[3],
                        p = out$`Pr(>F)`[3])
          return(out2)
          
        } 
        
      }) %>% bind_rows() 
    
    return(loop.cat.out)
    
  }) %>% bind_rows


# merge to images for each age, ethnicity, and income subgroup analyses
adonis.sub.age = bind_rows(adonis.uni.cat.age, adonis.uni.cont.age) %>% filter(p<0.05) %>% mutate(group = 'Age')
adonis.sub.inc = bind_rows(adonis.uni.cat.inc, adonis.uni.cont.inc) %>% filter(p<0.05) %>% mutate(group = 'Income')
adonis.sub.eth = bind_rows(adonis.uni.cat.eth, adonis.uni.cont.eth) %>% filter(p<0.05) %>% mutate(group = 'Ethnicity')
adonis.sub.bmi = bind_rows(adonis.uni.cat.bmi, adonis.uni.cont.bmi) %>% filter(p<0.05) %>% mutate(group = 'BMI') %>% mutate(subgroup = ifelse(grepl('obese', subgroup), 'obese', subgroup))
adonis.sub.sex = bind_rows(adonis.uni.cat.sex, adonis.uni.cont.sex) %>% filter(p<0.05) %>% mutate(group = 'Sex')

fig.adonis.subgroup = 
  bind_rows(adonis.sub.age, adonis.sub.inc, adonis.sub.eth, adonis.sub.bmi, adonis.sub.sex) %>%
  mutate(subgroup = ifelse(subgroup=='1', '<MYR400', ifelse(subgroup=='2', 'MYR401-700', ifelse(subgroup=='3', 'MYR701-1000', ifelse(subgroup=='4', 'MYR1001-5000', ifelse(
                            subgroup=='5', '>MYR5000', subgroup))))),
         subgroup = ifelse(subgroup=='Native', 'Jakun', subgroup)) %>%
  mutate(category = ifelse(grepl('water|toilet|utensil|handwashing', var), 'Hygiene', ifelse(
                            grepl('yeast|calbicans|mould|vre', var), 'Colonisation', ifelse(
                              grepl('med_suppl|diabetes|chronic|bpcat', var), 'Health', ifelse(
                                grepl('meat|fruit|fermented|coffee|probiotic', var), 'Diet', ifelse(
                                  grepl('ethnicity|bmi|smoker|Location|income|education|agecat', var), 'Demographic', 'Environment')))))) %>% 
  mutate(var = gsub(pattern = '_', replacement = ' ', var),
         var = str_to_sentence(var)) %>% 
  ggplot(aes(x=r2, y=var, fill=category)) +
  geom_col() +
  ggsci::scale_fill_d3(name='') + 
  facet_wrap(group~subgroup, scales='free_y', ncol=3) +
  theme_bw() +
  labs(x='PERMANOVA R<sup>2</sup>', y='') +
  theme(axis.title.x = ggtext::element_markdown(),
        legend.position = 'top')

##################################################################
### univariate permanova accounted for age, sex, and ethnicity ###
##################################################################

df.arg.cont.multi = 
  df.clean %>%
  mutate(across(.fns = ~ifelse(.=='N/A', NA, .)),
         gender = str_to_sentence(gender),
         household = as.factor(household),
         batch = as.factor(batch),
         income = as.factor(income)) %>%
  select(id, household, batch, age, gender, exercise, bmi, ethnicity, income, starts_with('house_'), ends_with('frequency')) %>%
  mutate(across(.fns = ~ifelse(is.na(.), median(., na.rm = T), .))) %>%
  column_to_rownames('id')

adonis.cont.multi.out = 
  lapply(names(df.arg.cont)[-c(1:2)], function(var) {
    
    otu = ps.arg.sub.cleandf.adoniscont %>% microbiome::transform('clr') %>% otu_table() %>% t
    form = paste0('otu ~ batch + household + age + gender + bmi + ethnicity + income + ', var) %>% as.formula()
    out = adonis2(formula = form, 
                  data = df.arg.cont.multi,
                  method = 'euclidean')
    out2 = tibble(var = var,
                  df = out$Df[8],
                  ss = out$SumOfSqs[8],
                  r2 = out$R2[8],
                  f = out$F[8],
                  p = out$`Pr(>F)`[8])
    return(out2)
    
  }) %>% bind_rows()

df.arg.cat.multi = 
  df.clean %>%
  mutate(across(.fns = ~ifelse(.=='N/A', NA, .)),
         gender = str_to_sentence(gender),
         household = as.factor(household),
         batch = as.factor(batch),
         income = as.factor(income),
         calbicans = ifelse(is.na(calbicans), 'Negative', calbicans),
         esbl = ifelse(is.na(esbl), 'Negative', esbl),
         yeast = ifelse(is.na(yeast), 'Negative', yeast)) %>%
  select(-age, -bmi, -exercise, -starts_with('house_'), -ends_with('frequency')) %>% 
  mutate(across(.fns = ~ifelse(is.na(.), Mode(.), .)),
         income = as.factor(income)) %>% column_to_rownames('id') %>% 
  mutate(age = df.arg.cont.multi$age,
         bmi = df.arg.cont.multi$bmi)

adonis.cat.multi.out = 
  lapply(names(df.arg.cat)[-c(1:2)], function(var) {
    
    otu = ps.arg.sub.cleandf.adoniscat %>% microbiome::transform('clr') %>% otu_table() %>% t
    form = paste0('otu ~ batch + household + age + gender + bmi + ethnicity + income + ', var) %>% as.formula()
    out = adonis2(formula = form, 
                  data = df.arg.cat.multi,
                  method = 'euclidean')
    out2 = tibble(var = var,
                  df = out$Df[8],
                  ss = out$SumOfSqs[8],
                  r2 = out$R2[8],
                  f = out$F[8],
                  p = out$`Pr(>F)`[8])
    return(out2)
    
  }) %>% bind_rows()

adonis.multi.out = bind_rows(adonis.cont.multi.out, adonis.cat.multi.out) %>% filter(p<0.05) %>% arrange(-r2)

################################################
### MAASLIN2 DIFFERENTIAL ABUNDANCE ANALYSIS ###
################################################

# differential abundance analysis of resistome genes
df.maaslin2 =
  cbind(df.arg.cat, df.arg.cont) %>% data.frame %>% mutate(ethnicity = ifelse(ethnicity == 'Native', 'Jakun', ethnicity)) %>%
  select(batch, household, age, bmi, ethnicity, income, adonis.multi.out$var) %>%
  mutate(yeast = ifelse(yeast=='N/A', 'Negative', yeast),
         calbicans = ifelse(calbicans=='N/A', 'Negative', calbicans)) %>%
  select(-calbicans)

ml2.res.out3 = 
  Maaslin2(input_data = ps.arg.sub %>% otu_table %>% data.frame(check.names = F),
           input_metadata = df.maaslin2,
           transform = 'LOG',
           normalization = 'NONE',
           analysis_method = 'LM',
           cores = 8,
           plot_heatmap = FALSE,
           plot_scatter = FALSE,
           fixed_effects = c('age', 'bmi', 'gender', 'ethnicity', 'income', adonis.multi.out$var[-4]),
           random_effects = c('batch', 'household'),
           reference = c('ethnicity,Jakun', 'income,1'),
           output = 'output/rev1/ml2_resistome5')

fig.arg.ml2 =
  ml2.res.out3$results %>% 
  filter(qval<0.1) %>% 
  mutate(metadata = ifelse(metadata=='animal_chicken', 'Chicken near house', ifelse(metadata=='bmi', 'BMI', ifelse(metadata=='chronic_diabetes', 'Diabetic', ifelse(metadata=='ethnicity', 'Ethnicity', ifelse(metadata=='income', 'Income', ifelse(metadata=='toilet_location_outside', 'Toilet outside main house structure', ifelse(metadata=='yeast', 'Stool yeast culture', 'Age')))))))) %>%
    mutate(value = ifelse(value=='bmi', '', ifelse(value=='age', '', ifelse(value=='2', 'MYR401-700', ifelse(value=='3', 'MYR701-1000', ifelse(value=='4', 'MYR1001-5000', value)))))) %>% 
  mutate(metadata_value = paste0(metadata, ': ', value)) %>%
  mutate(feature = gsub('macrolide\\.lincosamide\\.streptogramin', 'MLS', feature)) %>%
  ggplot(aes(x=coef, xmin = coef-stderr, xmax = coef+stderr, y=reorder(feature,coef), color = metadata_value)) +
  geom_pointrange(position = position_dodge(width = 1)) +
  geom_vline(xintercept = 0, linetype='dashed') + 
  theme_bw() +
  ggsci::scale_color_d3(name='', palette = 'category20') +
  labs(x='Maaslin2 Coefficient', y='Antibiotic resistance gene')


################################################
### SPECIES MAASLIN2 DIFFERENTIAL ABUNDANCE ANALYSIS ###
################################################

# differential abundance analysis of resistome genes
df.maaslin2.sp =
  cbind(df.arg.cat, df.arg.cont) %>% data.frame %>% mutate(ethnicity = ifelse(ethnicity == 'Native', 'Jakun', ethnicity)) %>%
  select(batch, household, age, bmi, ethnicity, income, adonis.multi.sp.out$var) %>%
  select(-house_06) # house_06 excluded as collinear with total house member

ml2.sp.out = 
  Maaslin2(input_data = ps %>% microbiome::core(detection=0.5/100, prevalence = 10/100) %>% otu_table %>% data.frame(check.names = F),
           input_metadata = df.maaslin2.sp,
           transform = 'LOG',
           normalization = 'NONE',
           analysis_method = 'LM',
           cores = 8,
           plot_heatmap = FALSE,
           plot_scatter = FALSE,
           fixed_effects = c('age', 'bmi', 'gender', 'ethnicity', 'income', adonis.multi.sp.out$var[-5]),
           random_effects = c('batch', 'household'),
           reference = c('ethnicity,Jakun', 'income,1', 'relationship,head'),
           output = 'output/rev1/ml2_sp')

fig.sp.ml2 =
  ml2.sp.out$results %>% 
  filter(qval<0.1) %>% 
  mutate(metadata = ifelse(metadata=='bmi', 'BMI', ifelse(metadata=='chronic_diabetes', 'Diabetic', ifelse(metadata=='ethnicity', 'Ethnicity', ifelse(metadata=='income', 'Income', ifelse(metadata=='relationship', 'Relationship to head of household', ifelse(metadata=='age', 'Age', ifelse(metadata=='surgery', 'Surgical history in the past year', ifelse(metadata=='chronic_others', 'Chronic disease', ifelse(metadata=='house_total', 'Number of household members', ifelse(metadata=='toilet_type_flush_sewer', 'Having flush sewer toilet', 'Undefined'))))))))))) %>%
  mutate(value = ifelse(value=='bmi', '', ifelse(value=='age', '', ifelse(value=='2', 'MYR401-700', ifelse(value=='3', 'MYR701-1000', ifelse(value=='4', 'MYR1001-5000', ifelse(value=='5', 'MYR5000', ifelse(value=='house_total', '', value)))))))) %>% 
  mutate(metadata_value = paste0(metadata, ': ', value)) %>%
  mutate(feature = gsub('s__', '', feature),
         feature = gsub('_', ' ', feature)) %>%
  mutate(genus = str_split_fixed(feature, pattern = ' ', n=2) %>% data.frame %>% pull(X1),
         sp = str_split_fixed(feature, pattern = ' ', n=4) %>% data.frame %>% pull(X2),
         sp2 = str_split_fixed(feature, pattern = ' ', n=4) %>% data.frame %>% pull(X3),
         sp3 = str_split_fixed(feature, pattern = ' ', n=4) %>% data.frame %>% pull(X4),
         genus_itl = paste0('*', genus, '*'),
         sp_itl = paste0('*', sp, '*'),
         fullname = paste0(genus_itl, ' ', sp_itl, ' ', sp2, sp3)) %>%
  ggplot(aes(x=coef, xmin = coef-stderr, xmax = coef+stderr, y=reorder(fullname,coef), color = metadata_value)) +
  geom_pointrange(position = position_dodge(width = 1)) +
  geom_vline(xintercept = 0, linetype='dashed') + 
  theme_bw() +
  ggsci::scale_color_d3(name='', palette = 'category20') +
  labs(x='Maaslin2 Coefficient', y='') +
    theme(axis.text.y = ggtext::element_markdown())

########################
### NETWORK ANALYSIS ###
########################

# get species maaslin2 output from former maaslin2 run
net_sp = ml2.sp.out$results %>% filter(qval<0.1) %>% pull(feature) %>% unique

# resistome -- run maaslin2 again but now revise taxa name to sequential list to ease extraction for network analysis
ps.rev = ps.arg.sub
taxa_names(ps.rev) = paste0('arg_', 1:length(taxa_names(ps.arg.sub)))

ml2.res.out3.rev = 
  Maaslin2(input_data = ps.rev %>% otu_table %>% data.frame(check.names = F),
           input_metadata = df.maaslin2,
           transform = 'LOG',
           normalization = 'NONE',
           analysis_method = 'LM',
           cores = 8,
           plot_heatmap = FALSE,
           plot_scatter = FALSE,
           fixed_effects = c('age', 'bmi', 'gender', 'ethnicity', 'income', adonis.multi.out$var[-4]),
           random_effects = c('batch', 'household'),
           reference = c('ethnicity,Jakun', 'income,1'),
           output = 'output/rev1/ml2_resistome5_rev')

arg_sig = ml2.res.out3.rev$results %>% filter(qval<0.1) %>% pull(feature) %>% unique
net_arg = tax_table(ps.rev) %>% data.frame %>% rownames_to_column('id') %>% filter(id %in% arg_sig) %>% mutate(name = paste(type,sub,  sep='__')) %>% pull(name)

# compile metadata var associated with resistome and sp alpha and beta div

# resistome beta
net_meta1 = adonis.multi.out %>% pull(var)
# sp beta 
net_meta2 = adonis.multi.sp.out %>% pull(var)
# resistome alpha
net_meta3 = lmer.arg.rich.out %>% arrange(p.value) %>% filter(p.value<0.1) %>% pull(var)
# sp alpha
net_meta4 = lmer.sp.rich.out %>% arrange(p.value) %>% filter(p.value<0.1) %>% pull(var)

net_meta = c(net_meta1, net_meta2, net_meta3, net_meta4) %>% unique

# metadata
netdf_meta = df.arg.cont.multi %>% rownames_to_column('id') %>% inner_join(rich.arg, by='id') %>% 
  inner_join(df.arg.cat.multi %>% rownames_to_column('id') %>% select(-household, -batch, -age, -gender, -bmi, -ethnicity, -income), by='id') %>%
  select(id, net_meta) %>%
  select(-relationship, -house_1317, -house_06) %>%
  rename(frequency_eat_out = dapao_frequency) %>% # remove categorical var since cannot be correlated
  mutate(across(!id, .fns = tolower)) %>% 
  mutate(across(!id, .fns = ~ifelse(.=='no', 0, .)),
         across(!id, .fns = ~ifelse(.=='none', 0, .)),
         across(!id, .fns = ~ifelse(.=='negative', 0, .)),
         across(!id, .fns = ~ifelse(.=='yes', 1, .)),
         across(!id, .fns = ~ifelse(.=='positive', 1, .)),
         across(!id, .fns = ~ifelse(.=='ex', 0, .))) %>%
  mutate(across(!id, .fns = as.numeric))

# species
netdf_sp = 
  prune_species(net_sp, ps) %>% otu_table %>% data.frame(check.names = F) %>% t %>% data.frame %>%
  rownames_to_column('id') %>%
  mutate(across(!id, .fns = function(x) { (x - min(x)) / (max(x) - min(x)) }))

# resistome
netdf_arg = 
  prune_taxa(net_arg, ps.arg.sub) %>% otu_table %>% data.frame(check.names = F) %>% t %>% data.frame %>%
  rownames_to_column('id') %>%
  mutate(across(!id, .fns = function(x) { (x - min(x)) / (max(x) - min(x)) }))

netdf_compiled = 
  netdf_meta %>% inner_join(netdf_sp, by='id') %>% inner_join(netdf_arg, by='id') %>%
  column_to_rownames('id')

# correlation
cor_out = netdf_compiled %>% cor(method = 'spearman') %>% reshape2::melt() %>% mutate(linker = paste0(Var1, '_', Var2)) %>% rename(cor=value) %>% select(-Var1, -Var2)
cor_p_out = netdf_compiled %>% cor.mtest(conf.level = 0.95) 
cor_p_out2 = cor_p_out$p %>% reshape2::melt() %>% mutate(linker = paste0(Var1, '_', Var2)) %>% rename(pval = value) 

unique_pair = names(netdf_compiled) %>% combn(2) %>% t %>% data.frame %>% dplyr::mutate(pair = paste(X1, X2, sep='_'))

cor_out_filt = 
  cor_out %>% inner_join(cor_p_out2, by='linker') %>% 
  filter(linker %in% unique_pair$pair) %>% 
  mutate(cat1 = ifelse(Var1 %in% names(netdf_meta), 'Metadata', ifelse(Var1 %in% names(netdf_sp), 'Species', ifelse(Var1 %in% names(netdf_arg), 'ARG', 'Others'))),
         cat2 = ifelse(Var2 %in% names(netdf_meta), 'Metadata', ifelse(Var2 %in% names(netdf_sp), 'Species', ifelse(Var2 %in% names(netdf_arg), 'ARG', 'Others')))) %>%
  filter(cat1 != cat2) %>%
  mutate(cat12 = paste(cat1, cat2, sep='_')) %>%
  filter(grepl('ARG', cat12)) %>%
  filter(pval<0.1, abs(cor)>0.25 ) %>%
  select(-linker, -pval, -cat1, -cat2, -cat12)

# remove variables with no networks
network_toremove =
  cor_out_filt[,c(2,3,1)] %>%
  as_tbl_graph() %>%
  dplyr::mutate(class = ifelse(name %in% names(netdf_meta), 'Metadata', ifelse(name %in% names(netdf_sp), 'Species', ifelse(name %in% names(netdf_arg), 'ARG', 'Others')))) %>%
  activate(edges) %>%
  dplyr::select(from, to) %>% as_data_frame()

tokeep = c(network_toremove$from, network_toremove$to) %>% unique

# generate graph

set.seed(3212)

# extrafont::font_import()
# extrafont::loadfonts()
fig.network =
  cor_out_filt[,c(2,3,1)] %>%
  as_tbl_graph() %>%
  dplyr::mutate(Class = ifelse(name %in% names(netdf_meta), 'Metadata', ifelse(name %in% names(netdf_sp), 'Species', ifelse(name %in% names(netdf_arg), 'ARG', 'Others')))) %>%
  dplyr::filter(name %in% tokeep) %>%
  activate(edges) %>%
  ggraph(layout='dh') +
  geom_node_point(aes(color=Class), size=4) +
  geom_edge_link(aes(edge_alpha=cor,
                     color=cor),
                 edge_width = 0.8) +
  guides(edge_width='none', edge_alpha='none') +
  geom_node_text(aes(label=name), size=2.5, repel=T) +
  scale_edge_color_gradientn(name = 'Correlation', colors = mixOmics::color.jet(9), limits = c(-1, 1)) +
  ggsci::scale_color_d3() +
  theme_graph() +
  scale_edge_width(range = c(0.2, 1.8))  +
  labs(caption='Significant correlation R>0.25')

fig.network 
ggsave('output/rev1/fig_network_rev.png', plot = fig.network, height=8, width=12)


# how many ARG linked with ecoli
cor_out_filt[,c(2,3,1)] %>% filter(grepl('s__Escherichia_coli', Var1) | grepl('s__Escherichia_coli', Var2))

# how many ARGS linked to yeast
cor_out_filt[,c(2,3,1)] %>% filter(grepl('yeast', Var1) | grepl('yeast', Var2))

# how many antibiotic genes per person
ps.arg.sub %>% otu_table %>% data.frame %>% t %>% data.frame %>% rownames_to_column() %>% tibble %>% pivot_longer(!rowname) %>%
  group_by(rowname) %>%
  summarise(sum = sum(value)) %>%
  pull(sum) %>% mean

ps.arg.sub %>% otu_table %>% data.frame %>% t %>% data.frame %>% rownames_to_column() %>% tibble %>% pivot_longer(!rowname) %>%
  group_by(rowname) %>%
  summarise(sum = sum(value)) %>%
  pull(sum) %>% sd

# visualise microbiome community

support.df = ps %>% tax_table() %>% data.frame %>% rownames_to_column('name')

fig.sp.top20 = 
  ps %>% microbiome::core(detection=0.5/100, prevalence = 10/100) %>% 
  otu_table() %>% 
  data.frame(check.names = F) %>% t %>%
  data.frame %>%
  rownames_to_column('id') %>%
  pivot_longer(!id) %>%
  group_by(name) %>%
  summarise(mean = mean(value),
            sd = sd(value)) %>%
  arrange(-mean) %>%
  slice_max(mean, n=20) %>%
  inner_join(support.df, by='name') %>%
  mutate(name = gsub('s__', '', name),
         name = gsub('_', ' ', name),
         name = paste0('*', name, '*'),
         phylum = gsub('p__', '', phylum)) %>%
  ggplot(aes(x=mean, xmin = mean-sd, xmax = mean+sd, y=reorder(name,mean), color=phylum)) +
  geom_pointrange() +
  ggsci::scale_color_d3(name='', palette = 'category20') + 
  theme_bw() +
  labs(x='Relative abundance (%)', y='') +
  theme(axis.text.y = ggtext::element_markdown(),
        legend.position = 'top')

##########################################################
### Microbial composition - Shannon diversity analysis ###
##########################################################

rich.sp = ps %>% estimate_richness(measures='Shannon') %>% rownames_to_column('id')

df.lmer.sp = 
  df.arg.cont.multi %>% rownames_to_column('id') %>% inner_join(rich.sp, by='id') %>% 
  inner_join(df.arg.cat.multi %>% rownames_to_column('id') %>% select(-household, -batch, -age, -gender, -bmi, -ethnicity, -income), by='id')

lmer.sp.rich.out = 
  lapply(names(df.lmer.sp)[-c(1:3)], function(x) {
  
    if(x == 'age') {
    form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
    form0 = paste0('Shannon ~ 1 + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
  } else if (x=='bmi') {
    form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
    form0 = paste0('Shannon ~ age + gender + ethnicity + 1 + income + (1|household) + (1|batch)') %>% as.formula
  } else if (x=='gender') {
    form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
    form0 = paste0('Shannon ~ age + 1 + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
  } else if (x=='ethnicity') {
    form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
    form0 = paste0('Shannon ~ age + gender + 1 + bmi + income + (1|household) + (1|batch)') %>% as.formula
  } else if (x=='income') {
    form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
    form0 = paste0('Shannon ~ age + gender + ethnicity + bmi + 1 + (1|household) + (1|batch)') %>% as.formula
  } else {
    form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + ', x, ' + (1|household) + (1|batch)') %>% as.formula
    form0 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + 1 + (1|household) + (1|batch)') %>% as.formula
  }
  
  model1 = lmer(form1, data = df.lmer.sp, REML=F)
  model0 = lmer(form0, data = df.lmer.sp, REML=F)
  
  out = anova(model1, model0) %>% tidy %>% mutate(var=x)
  return(out)
  
}) %>% bind_rows() %>% filter(term=='model1')

lmer.sp.rich.out %>% filter(p.value<0.1)

##########################################################
### Resistome  - Shannon diversity analysis ###
##########################################################

rich.arg

df.lmer.arg = 
  df.arg.cont.multi %>% rownames_to_column('id') %>% inner_join(rich.arg, by='id') %>% 
  inner_join(df.arg.cat.multi %>% rownames_to_column('id') %>% select(-household, -batch, -age, -gender, -bmi, -ethnicity, -income), by='id')

lmer.arg.rich.out = 
  lapply(names(df.lmer.arg)[-c(1:3)], function(x) {
    
    if(x == 'age') {
      form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
      form0 = paste0('Shannon ~ 1 + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
    } else if (x=='bmi') {
      form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
      form0 = paste0('Shannon ~ age + gender + ethnicity + 1 + income + (1|household) + (1|batch)') %>% as.formula
    } else if (x=='gender') {
      form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
      form0 = paste0('Shannon ~ age + 1 + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
    } else if (x=='ethnicity') {
      form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
      form0 = paste0('Shannon ~ age + gender + 1 + bmi + income + (1|household) + (1|batch)') %>% as.formula
    } else if (x=='income') {
      form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + (1|household) + (1|batch)') %>% as.formula
      form0 = paste0('Shannon ~ age + gender + ethnicity + bmi + 1 + (1|household) + (1|batch)') %>% as.formula
    } else {
      form1 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + ', x, ' + (1|household) + (1|batch)') %>% as.formula
      form0 = paste0('Shannon ~ age + gender + ethnicity + bmi + income + 1 + (1|household) + (1|batch)') %>% as.formula
    }
    
    model1 = lmer(form1, data = df.lmer.arg, REML=F)
    model0 = lmer(form0, data = df.lmer.arg, REML=F)
    
    out = anova(model1, model0) %>% tidy %>% mutate(var=x)
    return(out)
    
  }) %>% bind_rows() %>% filter(term=='model1')

lmer.arg.rich.out %>% filter(p.value<0.1)

###################################################################################################
### univariate permanova for species profile accounted for age, sex, bmi, ethnicity, and income ###
###################################################################################################

adonis.cont.multi.sp.out = 
  lapply(names(df.arg.cont)[-c(1:2)], function(var) {
    
    otu = ps %>% microbiome::core(detection = 0.5/100, prevalence = 10/100) %>% microbiome::transform('clr') %>% otu_table() %>% t
    form = paste0('otu ~ batch + household + age + gender + bmi + ethnicity + income + ', var) %>% as.formula()
    out = adonis2(formula = form, 
                  data = df.arg.cont.multi,
                  method = 'euclidean')
    out2 = tibble(var = var,
                  df = out$Df[8],
                  ss = out$SumOfSqs[8],
                  r2 = out$R2[8],
                  f = out$F[8],
                  p = out$`Pr(>F)`[8])
    return(out2)
    
  }) %>% bind_rows()

adonis.cat.multi.sp.out = 
  lapply(names(df.arg.cat)[-c(1:2)], function(var) {
    
    otu = ps %>% microbiome::core(detection = 0.5/100, prevalence = 10/100) %>% microbiome::transform('clr') %>% otu_table() %>% t
    form = paste0('otu ~ batch + household + age + gender + bmi + ethnicity + income + ', var) %>% as.formula()
    out = adonis2(formula = form, 
                  data = df.arg.cat.multi,
                  method = 'euclidean')
    out2 = tibble(var = var,
                  df = out$Df[8],
                  ss = out$SumOfSqs[8],
                  r2 = out$R2[8],
                  f = out$F[8],
                  p = out$`Pr(>F)`[8])
    return(out2)
    
  }) %>% bind_rows()

adonis.multi.sp.out = bind_rows(adonis.cont.multi.sp.out, adonis.cat.multi.sp.out) %>% filter(p<0.05) %>% arrange(-r2)

#########################################
### Figure preparation for manuscript ###
#########################################


# figure 1a - Abundance of antibiotic resistance genes from the studied cohort (N=200)
fig.arg.sum

# figure 1b - Abundance of top 20 antibiotic resistance genes detected in the Segamat cohort 
fig.arg.top20

# figure 1c - Factors associated with resistome abundance profile based on PERMANOVA adjusted for batch, household, age, sex, ethnicity, income, and BMI (p<0.05)
fig.arg.adonis = 
  adonis.multi.out %>%
  mutate(
    cat = ifelse(var %in% c('chronic_diabetes', 'med_suppl'), 'Health', 
      ifelse(var %in% c('animal_chicken', 'house_1317', 'toilet_location_outside'), 'Environment', 
        ifelse(var %in% c('yeast', 'calbicans'), 'Stool culture', 
          ifelse(var %in% c('meat_beef'), 'Diet', 'Undefined'          )        )      )    )  ) %>%
  mutate(var2 = ifelse(var=='chronic_diabetes', 'Diabetic', ifelse(var=='animal_chicken', 'Chicken near house', ifelse(var=='yeast', 'Stool yeast culture', ifelse(var=='med_suppl', 'On active medication', ifelse(var=='calbicans', 'Stool *C. albicans* culture', ifelse(var=='toilet_location_outside', 'Toilet outside main house structure', ifelse(var=='house_1317', 'No. of 13-17 y.o in house', ifelse(var=='meat_beef', 'Beef consumption habit', 'Undefined'))))))))) %>% 
  ggplot(aes(x=r2,y=reorder(var2,r2),fill=cat)) + 
  geom_col() +
  labs(x='PERMANOVA R<sup>2</sup>', y='') +
  ggsci::scale_fill_d3(name='') +
  theme_bw() +
  theme(axis.title.x = ggtext::element_markdown(),
        axis.text.y = ggtext::element_markdown()) 

fig1_comp1 = ggpubr::ggarrange(fig.arg.sum, fig.arg.adonis, labels = c('a', 'c'), nrow=2)
fig1_comp2 = ggpubr::ggarrange(fig1_comp1, fig.arg.top20, labels = c('', 'b'))

ggsave(file = 'output/rev1/figure1.pdf', plot = fig1_comp2, height = 8, width = 12)
ggsave(file = 'output/rev2/figure1.png', plot = fig1_comp2, height = 8, width = 12)
ggsave(file = 'output/rev2/figure1.svg', plot = fig1_comp2, height = 8, width = 12)
ggsave(file = 'output/rev3/figure1.pdf', plot = fig1_comp2, height = 8, width = 12)

# figure 2 - Differentially abundant antibiotic resistance genes against lifestyle and demographic factors of the Segamat cohort, analysed using Maaslin2 (q<0.1). Ethnicity and household income variables were compared against Jakun and those earning <MYR400, respectively
fig.arg.ml2
ggsave(file = 'output/rev1/figure2.pdf', plot = fig.arg.ml2, height = 8, width = 12)
ggsave(file = 'output/rev1/figure2.png', plot = fig.arg.ml2, height = 8, width = 12)
ggsave(file = 'output/rev2/figure2.svg', plot = fig.arg.ml2, height = 8, width = 12)
ggsave(file = 'output/rev3/figure2.pdf', plot = fig.arg.ml2, height = 8, width = 12)


# figure 3 - Differentially abundant species against lifestyle and demographic factors of the Segamat cohort, analysed using Maaslin2 (q<0.1). Ethnicity, and household income variables were compared against Jakun and those earning <MYR400, respectively
fig.sp.ml2
ggsave(file = 'output/rev1/figure3.pdf', plot = fig.sp.ml2, height = 8, width = 12)
ggsave(file = 'output/rev1/figure3.png', plot = fig.sp.ml2, height = 8, width = 12)
ggsave(file = 'output/rev2/figure3.svg', plot = fig.sp.ml2, height = 8, width = 12)
ggsave(file = 'output/rev3/figure3.pdf', plot = fig.sp.ml2, height = 8, width = 12)

# figure 4 - Network of metadata, species, and resistome features exhibiting significant spearman correlation (p<0.05), filtered for those with absolute correlation > 0.25
extrafont::font_import()
extrafont::loadfonts()
fig.network
ggsave(file = 'output/rev3/figure4_2.pdf', plot = fig.network, height = 8, width = 12)
ggsave(file = 'output/rev1/figure4.png', plot = fig.network, height = 8, width = 12)
ggsave(file = 'output/rev2/figure4.svg', plot = fig.network, height = 8, width = 12)
ggsave(file = 'output/rev3/figure4.eps', plot = fig.network, height = 8, width = 12)

pdf(file = 'output/rev3/figure4_rev.pdf', height = 8, width= 12)
fig.network
dev.off()


# supplementary figure 1 - Household income distribution of the studied cohort
household_uniq = df.clean %>% pull(household) %>% duplicated

fig.income =
  df.clean %>%
    filter(!household_uniq) %>%
  select(ethnicity, income) %>% 
    mutate(ethnicity = ifelse(ethnicity=='Native', 'Jakun', ethnicity)) %>%
  table %>% prop.table(margin=1) %>% data.frame %>%
  ggplot(aes(x=ethnicity, y=Freq, fill=income)) + geom_col() +
  theme_bw() +
  scale_fill_discrete(name='Income Level', labels = c('<400', '401-700', '701-1000', '1001-5000', '>5000')) +
  labs(x='', y='Proportion')


# supplementary figure 2 - Occupation distribution of the studied cohort classified based on sex
fig_occupation_strat =
  df.clean %>% select(ethnicity, gender, occupation) %>% na.omit() %>%
    mutate(ethnicity = ifelse(ethnicity=='Native', 'Jakun', ethnicity)) %>%
  mutate(ethnicity = str_to_sentence(ethnicity),
         occupation = str_to_sentence(occupation),
         gender = str_to_sentence(gender)) %>% count(occupation, ethnicity, gender) %>% 
  ggplot(aes(y=occupation, x=n, fill=ethnicity)) + geom_col(position= position_dodge(preserve = 'single')) +
  labs(x='', y='') +
  scale_fill_manual(name='Ethnicity', values = c('steelblue', 'orange', 'turquoise', 'indianred2')) +
  theme_bw() +
  facet_wrap(~gender) +
  geom_vline(xintercept = 0, linetype = 'dashed')

# supplementary figure 3 - Lifestyle and demographic factors significantly associated with the Shannon diversity index of the resistome profile (LRT p<0.1)
lmer.arg.rich.sig = lmer.arg.rich.out %>% filter(p.value<0.1) %>% pull(var)

fig.lmer.arg.rich.loop =
  lapply(lmer.arg.rich.sig, function(x) {
    
    if (x=='house_1317') {
      x = 'house_member_13_to_17_years_old_count'
    } else if (x=='dapao_frequency') {
      x = 'frequency_eat_out'
    } else if (x=='yeast') {
      x = 'stool_yeast_culture'
    } else if (x=='med_suppl') {
      x = 'on_active_medication'
    } else if (x=='chronic_diabetes') {
      x = 'diabetic'
    } else if (x=='animal_chicken') {
      x = 'chicken_near_house'
    } else if (x=='water_drinking_pipe') {
      x = 'piped_drinking_water_source'
    } else if (x=='toilet_location_outside') {
      x = 'toilet_outside_main_house'
    } else {
      x = x
    }
    
    df.loop = 
      df.lmer.arg %>% 
      rename(house_member_13_to_17_years_old_count=house_1317,
             frequency_eat_out=dapao_frequency,
             stool_yeast_culture=yeast,
             on_active_medication=med_suppl,
             diabetic=chronic_diabetes,
             chicken_near_house=animal_chicken,
             piped_drinking_water_source=water_drinking_pipe,
             toilet_outside_main_house=toilet_location_outside)
    
    if(is.numeric(df.loop[[x]])) {
      
      out = 
        df.loop %>%
        ggplot(aes_string(x=x, y='log2(Shannon)')) + 
        geom_point() +
        geom_smooth(method='lm') +
        theme_bw()
      
    } else {
      
      out =
        df.loop %>%
        ggplot(aes_string(x=x, y='Shannon')) +
        geom_jitter(width=0.1) +
        geom_boxplot(alpha=0.1) +
        theme_bw() 
    }
    
    return(out)
    
  })

fig.lmer.arg.rich = ggpubr::ggarrange(plotlist = fig.lmer.arg.rich.loop)

# supplementary figure 4 - Subgroup analysis of lifestyle factors and their association with the gut resistome profile based on adjusted PERMANOVA, stratified based on age, sex, ethnicity, income, and BMI
fig.adonis.subgroup_small = 
  fig.adonis.subgroup +
  theme(axis.text.y = element_text(size=5),
        strip.text.x.top = element_text(size=6))
  
# supplementary figure 5 - Top 20 most abundant microbial species in the stool culture of the Segamat cohort
fig.sp.top20

# supplementary figure 6 - Lifestyle and demographic factors significantly associated with the Shannon diversity index of the gut microbiota (LRT p<0.1)
lmer.sp.rich.sig = lmer.sp.rich.out %>% filter(p.value<0.1) %>% pull(var)

fig.lmer.sp.rich.loop =
  lapply(lmer.sp.rich.sig, function(x) {
  
  if(is.numeric(df.lmer.sp[[x]])) {
    
    out = 
      df.lmer.sp %>%
      ggplot(aes_string(x=paste0('log2(', x, ')'), y='Shannon')) + 
      geom_point() +
      geom_smooth(method='lm') +
      theme_bw()
    
  } else {
    
    out =
      df.lmer.sp %>%
      ggplot(aes_string(x=x, y='Shannon')) +
      geom_jitter(width=0.1) +
      geom_boxplot(alpha=0.1) +
      theme_bw() 
  }
  
  return(out)
    
})

fig.lmer.sp.rich = ggpubr::ggarrange(plotlist = fig.lmer.sp.rich.loop, nrow=1)


# supplementary figure 7 - Escherichia coli abundance and the resistome Shannon diversity
ecoli.df = ps %>% psmelt %>% filter(OTU == 's__Escherichia_coli') %>% select(Sample, Abundance) %>% rename(id=Sample)
ecoli.df2 = rich.arg %>% rename(Shannon_resistome=Shannon) %>% left_join(rich.sp) %>% rename(Shannon_sp = Shannon) %>% inner_join(ecoli.df)

fig.arg.ecoli = 
  ecoli.df2 %>% filter(Abundance!=0) %>% ggplot(aes(x=log2(Abundance+0.000001), y=Shannon_resistome)) + geom_point() + 
  ggpubr::stat_regline_equation(label.x = -10, label.y = 5) +
  ggpubr::stat_cor(label.x = -10, label.y = 4.85) +
  theme_bw() +
  geom_smooth(method = 'lm') +
  labs(x='Log<sub>2</sub> *E. coli* relative abundance (%)', y='Resistome Shannon index')+
  theme(axis.title.x = ggtext::element_markdown())

# compile all supplementary figures as pdf
pdf(file = 'output/rev3/supplementary_figures_compiled.pdf', paper = 'a4')

ggpubr::annotate_figure(fig.income, 
                        bottom = ggpubr::text_grob('Supplementary Figure 1\nHousehold income distribution of the studied cohort.', 
                                                   size=10))

ggpubr::annotate_figure(fig_occupation_strat, 
                        bottom = ggpubr::text_grob('Supplementary Figure 2\nOccupation distribution of the studied cohort classified based on sex.', size=10))

ggpubr::annotate_figure(fig.lmer.arg.rich, bottom = ggpubr::text_grob(
  'Supplementary Figure 3\nLifestyle and demographic factors significantly associated with the Shannon diversity index of the resistome\nprofile (LRT p<0.1). The boxplot’s lower and upper boundaries marked the first (25th percentile) and the third\n(75th percentile) quartile of the visualised data, while the middle hinge marked the median. The lower and upper\nwhiskers marked the lowest and highest values no smaller/larger than 1.5×interquartile range of the visualised\ndataset, respectively.', size=10))

ggpubr::annotate_figure(fig.adonis.subgroup_small, bottom = ggpubr::text_grob(
  'Supplementary Figure 4\nSubgroup analysis of lifestyle factors and their association with the gut resistome profile based on adjusted\nPERMANOVA, stratified based on age, sex, ethnicity, income, and BMI.', size=10))

ggpubr::annotate_figure(fig.sp.top20, bottom = ggpubr::text_grob(
  'Supplementary Figure 5\nTop 20 most abundant microbial species in the stool culture of the Segamat cohort. The point refers to the mean\nof the visualised dataset, while the error bar refers to 1× standard deviation of the visualised dataset.', size=10))

ggpubr::annotate_figure(fig.lmer.sp.rich, bottom = ggpubr::text_grob(
  'Supplementary Figure 6\nLifestyle and demographic factors significantly associated with the Shannon diversity index of the gut microbiota.\nThe boxplot’s lower and upper boundaries marked the first (25th percentile) and the third (75th percentile) quartile\nof the visualised data, while the middle hinge marked the median. The lower and upper whiskers marked the\nlowest and highest values no smaller/larger than 1.5×interquartile range of the visualised dataset, respectively.', size=10))

ggpubr::annotate_figure(fig.arg.ecoli, bottom = ggpubr::text_grob(
  'Supplementary Figure 7\nEscherichia coli abundance and the resistome Shannon diversity.', size=10))

dev.off()

# supplementary table 1 - Demographic factors of the studied cohort
df.demo

# supplementary table 2 - List of core antibiotic resistance genes detected from the Segamat cohort obtained based on the output of ARGs-OAP pipeline. 
table.arg.filt = 
  prune_taxa(arg.filt.names, ps.arg.sub) %>% psmelt %>% select(OTU, Sample, Abundance, type, sub) %>% group_by(OTU) %>%
  summarise(mean_abundance = mean(Abundance),
            sd_abundance = sd(Abundance)) %>%
  separate(OTU, into = c('class', 'gene'), sep='__') %>%
  arrange(-mean_abundance)

# supplementary table 3 - Lifestyle and demographic factors tested against resistome Shannon diversity (Linear mixed model, LRT p<0.1)
lmer.arg.rich.out %>% arrange(p.value) 

# supplementary table 4 - PERMANOVA analysis between antibiotic resistance genes and demographic and lifestyle parameters based on the Euclidean distance with 999 permutations, adjusted for household, batch, age, sex, ethnicity, income, and BMI). Resistance gene profile was centred log-ratio transformed before analysis
adonis.multi.out

# supplementary table 5 - Antibiotic resistance genes differentially abundant across lifestyle parameters associated with the resistome profile, analysed using Maaslin2 adjusted for batch and household as random effects and age, sex, ethnicity, BMI, and income as fixed effects. 
ml2.res.out3$results %>% filter(qval<0.1)

# supplementary table 6 - Lifestyle and demographic factors tested against the gut microbiome Shannon diversity (linear mixed model)
lmer.sp.rich.out %>% arrange(p.value) %>% filter(p.value<0.1)

# supplementary table 7- Demographic, batch and household-adjusted PERMANOVA identifying lifestyle variables associated with the gut microbiota profile of the Segamat community
adonis.multi.sp.out

# compile supplementary tables
library(openxlsx)

suppl_tab_list = list(
  suppl_tab1 = df.demo,
  suppl_tab2 = table.arg.filt,
  suppl_tab3 = lmer.arg.rich.out %>% arrange(p.value),
  suppl_tab4 = adonis.multi.out,
  suppl_tab5 = ml2.res.out3$results,
  suppl_tab6 = lmer.sp.rich.out %>% arrange(p.value),
  suppl_tab7 = adonis.multi.sp.out
)

openxlsx::write.xlsx(suppl_tab_list, file = 'output/rev1/supplementary_table_compiled.xlsx')

# direct conversion to pdf, ugly output
library(gridExtra)

pdf(file = 'output/rev1/supplementarty_table1.pdf', height=nrow(df.demo)/3, width = 15)
grid.table(df.demo)
dev.off()
pdf(file = 'output/rev1/supplementarty_table2.pdf', height=nrow(table.arg.filt)/3, width = 15)
grid.table(table.arg.filt)
dev.off()
pdf(file = 'output/rev1/supplementarty_table3.pdf', height=nrow(lmer.arg.rich.out)/3, width = 15)
grid.table(lmer.arg.rich.out %>% arrange(p.value))
dev.off()
pdf(file = 'output/rev1/supplementarty_table4.pdf', height=nrow(adonis.multi.out)/3, width = 10)
grid.table(adonis.multi.out)
dev.off()
pdf(file = 'output/rev1/supplementarty_table5.pdf', height=nrow(ml2.res.out3$results %>% filter(qval<0.1))/3, width = 20)
grid.table(ml2.res.out3$results %>% filter(qval<0.1))
dev.off()
pdf(file = 'output/rev1/supplementarty_table6.pdf', height=nrow(lmer.sp.rich.out)/3, width = 10)
grid.table(lmer.sp.rich.out %>% arrange(p.value))
dev.off()
pdf(file = 'output/rev1/supplementarty_table7.pdf', height=nrow(adonis.multi.sp.out)/3, width = 10)
grid.table(adonis.multi.sp.out)
dev.off()


# save xls sheet to pdf manually, then combine them here

qpdf::pdf_combine(input = c('output/rev1/supplementary_table1.pdf',
                           'output/rev1/supplementary_table2.pdf',
                           'output/rev1/supplementary_table3.pdf',
                           'output/rev1/supplementary_table4.pdf',
                           'output/rev1/supplementary_table5.pdf',
                           'output/rev1/supplementary_table6.pdf',
                           'output/rev1/supplementary_table7.pdf'),
                 output = 'output/rev1/supplementary_table_compiled.pdf')


save.image('env_rev1.Rdata')



