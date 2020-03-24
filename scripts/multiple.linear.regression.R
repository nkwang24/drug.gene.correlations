### multiple.linear.regression.R ##################################################################
#
# Nicholas Wang
# 2/10/20
# 
# Multiple lasso regression 
#

## Load packages ----------------------------------------------------------------------------------
library(tidyverse);
library(readxl);
library(glmnet);
library(ggplot2);

## Read data --------------------------------------------------------------------------------------
# Import GDSC2 drug data
df.drug <- read_excel('./data/GDSC2_fitted_dose_response_15Oct19.xlsx') %>% 
  select(DRUG_NAME, COSMIC_ID, LN_IC50) %>% 
  distinct(DRUG_NAME, COSMIC_ID, .keep_all = TRUE) %>% 
  pivot_wider(names_from = DRUG_NAME, values_from = LN_IC50) %>% 
  select(names(.)[1], sort(names(.)));

# Import Broad CRISPR data
df.gene <- read_csv('./data/Achilles_gene_effect.csv') %>% 
  drop_na();

# Cell line info
df.sample <- read_csv('./data/sample_info.csv', col_types = cols(additional_info = col_character())) %>% 
  select(DepMap_ID, stripped_cell_line_name, COSMIC_ID);

## Build model ------------------------------------------------------------------------------------
# Merge single drug with gene data
drug <- 'Venetoclax';
df.drug.1 <- merge(df.sample, select(df.drug, COSMIC_ID, drug), by = 'COSMIC_ID', all.y = TRUE) %>% 
  merge(df.gene, by.x = 'DepMap_ID', by.y = 'X1') %>% 
  filter(!is.na(!!as.name(drug))) %>% 
  drop_na() %>% 
  as_tibble();

lambdas <- 10^seq(3, -2, by = -0.1);
cvfit <- cv.glmnet(
  as.matrix(df.drug.1[-(1:4)]), 
  pull(df.drug.1, drug), 
  alpha = 1, 
  lambda = lambdas
  );
plot(cvfit);

fit.1se <- glmnet(
  as.matrix(df.drug.1[-(1:4)]), 
  pull(df.drug.1, drug), 
  alpha = 1, 
  lambda = cvfit$lambda.1se
  );
fit.min <- glmnet(
  as.matrix(df.drug.1[-(1:4)]), 
  pull(df.drug.1, drug), 
  alpha = 1, 
  lambda = cvfit$lambda.min
  );

betas <- coef(fit.1se) %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>% 
  filter(s0 != 0, gene != '(Intercept)') %>% 
  arrange(desc(s0));
cors <- cor(df.drug.1[drug], df.drug.1[betas$gene]) %>%
  t();
out <- cbind(betas, cors) %>% 
  rename(beta = s0, cor = drug);
rownames(out) <- NULL;
out

all.cors <- cor(df.drug.1[drug], df.drug.1[-(1:4)]) %>% 
  t() %>% 
  as.matrix() %>% 
  as.data.frame() %>% 
  rownames_to_column('gene') %>% 
  arrange(desc(get(drug)))
head(all.cors, 10)

# Scatter plot of drug-gene pair
plot(df.drug.1[[drug]], df.drug.1[['BRAF (673)']])

ggplot(out, aes(seq_len(length(beta)))) +
  geom_line(aes(y = beta, color = 'beta')) +
  geom_line(aes(y = cor, color = 'cor'))



