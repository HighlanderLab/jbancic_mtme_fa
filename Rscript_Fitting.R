##################################################################
# 
# R script used to fit linear mixed models in Bancic et al. 2022
# 
# Script authors: D.J. Tolhurst and J. Bancic
# 
##################################################################

# Assume that data.df is a multi-trait MET data-frame with the following columns
# Pheno - numeric, n x s observations scaled to unit variance for each trait (stacked trait information)
# Trait - factor, s = 3 levels (YLD, PHT, DTF)
# Env - factor, p = 12 levels
# TraitEnv - factor, s x p = 36 levels, given by the paste of the Trait and Env columns
# ColBlock - factor, cb levels
# RowBlock - factor, rb levels
# Column - factor, c levels
# Row - factor, r levels
# Gkeep - factor, v = 266 levels, and corresponds to those genotypes with marker data
# Gdrop - factor, vd = 25 levels, and corresponds to those genotypes without marker data 
# (Gdrop is not required if all genotypes have marker data)
# Note that data.df is ordered as environments within traits

# Also assume that:
# cb is an object which contains the names of environments in which column block effects are fitted
# rb is an object which contains the names of environments in which row block effects are fitted

# Lastly, assume that Gg is a v x v genomic relationship with column and row names ordered according to Gkeep

#########################################################################################################
# 
# 1 Non-separable models
# 
#########################################################################################################

#--------------------------------------------------------------------------
# diag*
#--------------------------------------------------------------------------

# This includes non-separable diagonal models for the additive and non-additive GET effects,
# as well as a non-separable residual model

# -- Run the model
diagStar.asr <- asreml(Pheno ~ Trait*Env,
                      random = ~ diag(TraitEnv):vm(Gkeep, Gg) +
                                 diag(TraitEnv):ide(Gkeep) +
                                 at(TraitEnv, cb):ColBlock +
                                 at(TraitEnv, rb):RowBlock,
                      residual = ~ dsum(~ar1(Column):ar1(Row) | TraitEnv),
                      sparse = ~ TraitEnv:Gdrop,
                      data = data.df, na.action = na.method(x="include"),
                      maxit = 13, workspace = 6e8)
diagStar.asr <- update(diagStar.asr)


#--------------------------------------------------------------------------
# diag
#--------------------------------------------------------------------------

# This includes non-separable diagonal models for the additive and non-additive GET effects,
# as well as a separable residual model

# -- Obtain start values
diag.sv <- asreml(Pheno  ~ 1,
                  random = ~ diag(TraitEnv):vm(Gkeep, Gg) +
                             diag(TraitEnv):ide(Gkeep) +
                             at(TraitEnv, cb):ColBlock +
                             at(TraitEnv, rb):RowBlock +
                             # Residual terms (Eq. 5)
                             us(Trait):at(Env):ar1v(Column):ar1(Row) +  
                             at(Trait, c(1,3)):Environment:units +   
                             Trait:diag(Environment):units,
                    family = asr_gaussian(dispersion = 1e-5),
                    data = data.df, na.action = na.method(x="include"),
                    start.values = T)
diag.sv <- diag.sv$vparameters.table
rownames(diag.sv) <- diag.sv$Component

# -- Constrain appropriate terms
# Note that there are only (p+s-1) degrees of freedom to estimate  
# at(Trait, c(1,3)):Environment:units  and  Trait:diag(Environment):units
# since they are linear combinations of one another. 
# So we cannot fit the random error term for one trait or one environment.
# We choose DTF since preliminary analysis revealed no random error for this trait

# Constraints due to separable residual structure
diag.sv$Value[grep("YLD:YLD", diag.sv$Component)] <- 1
diag.sv$Constraint[grep("YLD:YLD", diag.sv$Component)] <- "F"

# The process below produces separate spatial models for each environment, 
# but then the residual variances and covariances within and between traits
# are constrained to be equal across environments. Any additional random error
# specific to traits and environments is then captured by the terms including "units"
p <- nlevels(data.df$Env)
s <- nlevels(data.df$Trait)
vcc <- cbind(matrix(1:6, nrow = 6*p),1)
rownames(vcc) <- diag.sv$Component[grep("Column:Row.*Trait_", diag.sv$Component)]
# vcc is a matrix with two columns, the first column defines the
# grouping structure such that components with the same number will be constrained equal.
# Note that there are 6 (unique) combinations of 3 traits.
# The second column is used for scaling, and set to 1 in this case.
# Lastly, note that the rownames of vcc must exactly match the names in 
# diag.sv$Component that are to be constrained

# -- Run the model
diag.asr <- asreml(Pheno ~ Trait*Env,
                   random = ~ diag(TraitEnv):vm(Gkeep, Gg) +
                              diag(TraitEnv):ide(Gkeep) +
                              at(TraitEnv, cb):ColBlock +
                              at(TraitEnv, rb):RowBlock +
                              # Residual terms (Eq. 5)
                              us(Trait):at(Env):ar1v(Column):ar1(Row) + 
                              at(Trait, c(1,3)):Environment:units + 
                              Trait:diag(Environment):units,
                    sparse = ~TraitEnv:Gdrop,
                    family = asr_gaussian(dispersion = 1e-5),
                    G.param = diag.sv, R.param = diag.sv, 
                    vcc = vcc, 
                    data = data.df, na.action = na.method(x="include"),
                    maxit = 13, workspace = 6e8)
diag.asr <- update(diag.asr)


#--------------------------------------------------------------------------
# Non-separable factor analytic linear mixed model (NFA-LMM)
#--------------------------------------------------------------------------

# In this model, the covariances between GET effects are based on a non-separable factor analytic model 
# between traits and environments (Eq. 18).

# Note that the NFAk model is fitted in terms of its two components, 
# that is a reduced rank (rr) + diagonal (diag) model

# -- Run the model
NFA1.asr <- asreml(Pheno ~ Trait*Env,
                   random = ~ rr(TraitEnv,1):vm(Gkeep, Gg) + diag(TraitEnv):vm(Gkeep, Gg) +
                              diag(TraitEnv):ide(Gkeep) +
                              at(TraitEnv, cb):ColBlock +
                              at(TraitEnv, rb):RowBlock +
                              # Residual terms (Eq. 5)
                              us(Trait):at(Env):ar1v(Column):ar1(Row) + 
                              at(Trait, c(1,3)):Environment:units + 
                              Trait:diag(Environment):units,
                   sparse = ~TraitEnv:Gdrop,
                   family = asr_gaussian(dispersion=1e-5),
                   G.param = diag.sv, R.param = diag.sv, # use residual constraints from above
                   vcc = vcc,                            # use vcc matrix from above
                   data = data.df, na.action = na.method(x="include"),
                   maxit = 13, workspace = 6e8)
NFA1.asr <- update(NFA1.asr)


#########################################################################################################
# 
# 2 Partially separable models
# 
#########################################################################################################

#--------------------------------------------------------------------------
# diag:diag
#--------------------------------------------------------------------------

# -- Obtain start values
diagdiag.sv <- asreml(Pheno ~ 1, 
                      random = ~ diag(Trait):diag(Env):vm(Gkeep, Gg) +
                                 diag(TraitEnv):ide(Gkeep) +
                                 at(TraitEnv, cb):ColBlock +
                                 at(TraitEnv, rb):RowBlock +
                                 # Residual terms (Eq. 5)
                                 us(Trait):at(Env):ar1v(Column):ar1(Row) +
                                 at(Trait, c(1,3)):Environment:units + 
                                 Trait:diag(Environment):units,
                     family = asr_gaussian(dispersion = 1e-5),
                     data = data.df, na.action = na.method(x="include"),
                     start.values = T)
diagdiag.sv  <- diagdiag.sv$vparameters.table
rownames(diagdiag.sv) <- diagdiag.sv$Component

# -- Constrain appropriate terms
# Constraints due to separable structure for additive GET effects
diagdiag.sv$Value[grep("Trait.*Env.*Gkeep.*YLD$", diagdiag.sv$Component)] <- 1
diagdiag.sv$Constraint[grep("Trait.*Env.*Gkeep.*YLD$", diagdiag.sv$Component)] <- "F"
# Constraints due to separable structure for residual
diagdiag.sv$Value[grep("YLD:YLD", diagdiag.sv$Component)] <- 1
diagdiag.sv$Constraint[grep("YLD:YLD", diagdiag.sv$Component)] <- "F"

# -- Run the model
diagdiag.asr <- asreml(Pheno ~ Trait*Env,
                       random = ~ diag(Trait):diag(Env):vm(Gkeep, Gg) +
                                  diag(TraitEnv):ide(Gkeep) +
                                  at(TraitEnv, cb):ColBlock +
                                  at(TraitEnv, rb):RowBlock +
                                  # Residual terms (Eq. 5)
                                  us(Trait):at(Env):ar1v(Column):ar1(Row) +
                                  at(Trait, c(1,3)):Environment:units + 
                                  Trait:diag(Environment):units,
                      sparse = ~ TraitEnv:Gdrop,
                      family = asr_gaussian(dispersion = 1e-5),
                      G.param = diagdiag.sv, R.param = diagdiag.sv, 
                      vcc = vcc,   # use vcc matrix from above
                      data = data.df, na.action = na.method(x="include"),
                      maxit = 13, workspace = 6e8)
diagdiag.asr <- update(diagdiag.asr)


#--------------------------------------------------------------------------
# us:diag
#--------------------------------------------------------------------------

# -- Obtain start values
usdiag.sv <- asreml(Pheno ~ 1, 
                    random = ~ us(Trait):diag(Env):vm(Gkeep, Gg) +
                               diag(TraitEnv):ide(Gkeep) +
                               at(TraitEnv, cb):ColBlock +
                               at(TraitEnv, rb):RowBlock +
                               # Residual terms (Eq. 5)
                               us(Trait):at(Env):ar1v(Column):ar1(Row) +
                               at(Trait, c(1,3)):Environment:units + 
                               Trait:diag(Environment):units,
                    family = asr_gaussian(dispersion = 1e-5),
                    data = data.df, na.action = na.method(x="include"),
                    start.values = T)
usdiag.sv  <- usdiag.sv$vparameters.table
rownames(usdiag.sv) <- usdiag.sv$Component

# -- Constrain appropriate terms
# Constraints due to separable structure for additive GET effects and residual
usdiag.sv$Value[grep("YLD:YLD", usdiag.sv$Component)] <- 1
usdiag.sv$Constraint[grep("YLD:YLD", usdiag.sv$Component)] <- "F"

# -- Run the model
usdiag.asr <- asreml(Pheno ~ Trait*Env,
                     random = ~ us(Trait):diag(Env):vm(Gkeep, Gg) +
                                diag(TraitEnv):ide(Gkeep) +
                                at(TraitEnv, cb):ColBlock +
                                at(TraitEnv, rb):RowBlock +
                                # Residual terms (Eq. 5)
                                us(Trait):at(Env):ar1v(Column):ar1(Row) +
                                at(Trait, c(1,3)):Environment:units + 
                                Trait:diag(Environment):units,
                     sparse = ~ TraitEnv:Gdrop,
                     family = asr_gaussian(dispersion = 1e-5),
                     G.param = usdiag.sv, R.param = usdiag.sv, 
                     vcc = vcc,   # use vcc matrix from above
                     data = data.df, na.action = na.method(x="include"),
                     maxit = 13, workspace = 6e8)
usdiag.asr <- update(usdiag.asr)


#--------------------------------------------------------------------------
# Unstructured factor analytic linear mixed model (UFA-LMM)
#--------------------------------------------------------------------------

# In this model, the covariances between GET effects are based on an unstructured factor analytic model 
# between traits and environments (Eq. 19).

# -- Obtain start values
UFA1.sv <- asreml(Pheno ~ 1, 
                  random = ~ us(Trait):fa(Env, 1):vm(Gkeep, Gg) + diag(Trait):id(Env):vm(Gkeep, Gg) +
                             diag(TraitEnv):ide(Gkeep) +
                             at(TraitEnv, cb):ColBlock +
                             at(TraitEnv, rb):RowBlock +
                             # Residual terms (Eq. 5)
                             us(Trait):at(Env):ar1v(Column):ar1(Row) +
                             at(Trait, c(1,3)):Environment:units + 
                             Trait:diag(Environment):units,
                  family = asr_gaussian(dispersion = 1e-5),
                  data = data.df, na.action = na.method(x="include"),
                  start.values = T)
UFA1.sv  <- UFA1.sv$vparameters.table
rownames(UFA1.sv) <- UFA1.sv$Component

# -- Constrain appropriate terms (see Smith et al., 2019)
# Constraints due to separable structure for additive GET effects and residual
UFA1.sv$Value[grep("YLD:YLD", UFA1.sv$Component)] <- 1
UFA1.sv$Constraint[grep("YLD:YLD", UFA1.sv$Component)] <- "F"
# constrain the specific variances to zero (small number) 
UFA1.sv$Value[grep("Trait.*Env.*Gkeep.*var$", UFA1.sv$Component)] <- 1e-5 
UFA1.sv$Constraint[grep("Trait.*Env.*Gkeep.*var$", UFA1.sv$Component)] <- "F"

# -- Run the model
UFA1.asr <- asreml(Pheno ~ Trait*Env,
                   random = ~ us(Trait):fa(Env,1):vm(Gkeep, Gg) + diag(Trait):idv(Env):vm(Gkeep, Gg) + 
                              diag(TraitEnv):ide(Gkeep) +
                              at(TraitEnv, cb):ColBlock +
                              at(TraitEnv, rb):RowBlock +
                              # Residual terms (Eq. 5)
                              us(Trait):at(Env):ar1v(Column):ar1(Row) +
                              at(Trait, c(1,3)):Environment:units + 
                              Trait:diag(Environment):units,
                   sparse = ~ TraitEnv:Gdrop,
                   family = asr_gaussian(dispersion = 1e-5),
                   G.param = UFA1.sv, R.param = UFA1.sv, 
                   vcc = vcc,   # use vcc matrix from above
                   data = data.df, na.action = na.method(x="include"),
                   maxit = 13, workspace = 6e8)
UFA1.asr <- update(UFA1.asr)


#--------------------------------------------------------------------------
# Partially separable factor analytic linear mixed model (SFA-LMM)
#--------------------------------------------------------------------------

# In this model, the covariances between GET effects are based on a partially separable factor analytic 
# model between traits and environments (Eq. 20).

# -- Obtain start values
SFA1.sv <- asreml(Pheno ~ 1, 
                  random = ~ fa(Trait, 1):fa(Env, 1):vm(Gkeep, Gg) +
                             diag(Trait):idv(Env):vm(Gkeep, Gg) +
                             at(Env, c(1:9,11:12)):idv(Trait):vm(Gkeep, Gg) + 
                             diag(TraitEnv):ide(Gkeep) +
                             at(TraitEnv, cb):ColBlock +
                             at(TraitEnv, rb):RowBlock +
                             # Residual
                             us(Trait):at(Env):ar1v(Column):ar1(Row) +
                             at(Trait, c(1,3)):Environment:units + 
                             Trait:diag(Environment):units,
                  family = asr_gaussian(dispersion = 1e-5),
                  data = data.df, na.action = na.method(x="include"),
                  start.values = T)
SFA1.sv  <- SFA1.sv$vparameters.table
rownames(SFA1.sv) <- SFA1.sv$Component

# -- Constrain appropriate terms (see Appendix B)
# Note that there are only (p+s-1) degrees of freedom to estimate  
# diag(Trait):idv(Env):vm(Gkeep, Gg)  and  at(Env, c(1:9,11:12)):idv(Trait):vm(Gkeep, Gg)
# since they are linear combinations of one another. 
# So we cannot fit the residual GET effects for one trait or one environment.
# We choose Environment 10 since preliminary analysis revealed no residual genetic error for this Environment

# Constraints due to separable structure for additive GET effects
UFA1.sv$Value[grep("YLD.*fa1$", UFA1.sv$Component)] <- 1
UFA1.sv$Constraint[grep("YLD.*fa1$", UFA1.sv$Component)] <- "F"
# Constraints due to separable structure for residual
UFA1.sv$Value[grep("YLD:YLD", UFA1.sv$Component)] <- 1
UFA1.sv$Constraint[grep("YLD:YLD", UFA1.sv$Component)] <- "F"
# constrain the trait specific variances to zero
SFA1.sv$Value[grep("Trait.*Env.*Gkeep.*YLD.*var$|Trait.*Env.*Gkeep.*DTF.*var$|Trait.*Env.*Gkeep.*PHT.*var$", SFA1.sv$Component)] <- 0
SFA1.sv$Constraint[grep("Trait.*Env.*Gkeep.*YLD.*var$|Trait.*Env.*Gkeep.*DTF.*var$|Trait.*Env.*Gkeep.*PHT.*var$", SFA1.sv$Component)] <- "F"
# constrain the environmental specific variances to zero (small number) 
SFA1.sv$Value[grep("Trait.*Env.*Gkeep.*MB.*var$|Trait.*Env.*Gkeep.*MV.*var$", SFA1.sv$Component)] <- 1e-5 
SFA1.sv$Constraint[grep("Trait.*Env.*Gkeep.*var$", SFA1.sv$Component)] <- "F"

# -- Run the model
SFA1.asr <- asreml(Pheno ~ Trait*Env,
                   random = ~ fa(Trait, 1):fa(Env, 1):vm(Gkeep, Gg) +
                              diag(Trait):idv(Env):vm(Gkeep, Gg) +
                              at(Env, c(1:9,11:12)):idv(Trait):vm(Gkeep, Gg) + 
                              diag(TraitEnv):ide(Gkeep) +
                              at(TraitEnv, cb):ColBlock +
                              at(TraitEnv, rb):RowBlock +
                              # Residual
                              us(Trait):at(Env):ar1v(Column):ar1(Row) +
                              at(Trait, c(1,3)):Environment:units + 
                              Trait:diag(Environment):units,
                   sparse = ~ TraitEnv:Gdrop,
                   family = asr_gaussian(dispersion = 1e-5),
                   G.param = SFA1.sv, R.param = SFA1.sv, 
                   vcc= vcc,   # use vcc matrix from above
                   data = data.df, na.action = na.method(x="include"),
                   maxit = 13, workspace = 6e8)
SFA1.asr <- update(SFA1.asr)

