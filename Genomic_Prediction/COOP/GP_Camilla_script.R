      #############################################
      #  FA models with generalized main effects  #
      #          Cross-Validation: CV1            #
      #            Additive effects               #   
      #                                           #
      #          Author: Jenifer Camila           #
      #                  2024                     #
      #############################################
      
      #------------------------------------------------------------------------------------------------------------------------------------------------------------------#
      # Removing variables:
      rm(list=ls())

      #------------------------------------------------------------------------------------------------------------------------------------------------------------------#
      # Assessment of predictive accuracy in CV1:
      
      # Ordering:
      data = data[order(data$city),]
      
      # Sample per group:
      k.folds = 5
      rep = 5   
      k = 5
      
      samp = data[!duplicated(data$code),]
      samp = samp %>% mutate(group = ifelse(hg == "C x Dent" | is.na(hg), 1, ifelse(hg == "Dent x Flint", 1, 2)))
      list.folder = cv_1(samp, k.folds, rep)
      
      # Loop:
      acuracy = list()
      results = list()
      
      for (i in 1:rep){
        
           for(j in 1:k.folds){
             
           training_data = data
           
           # Including NA to the response variable in the validation set:
           training_data[training_data$code %in% list.folder[[i]][[j]], "grain_yield"] = NA
           
           # GWS model:
          model_fa = asreml(fixed = grain_yield ~ city + city:trial + city:trial:rep, 
                            random = ~ rr(city, 5):vm(code, h_inv) + 
                            diag(city):vm(code, h_inv) +
                            diag(city):rep:block,
                            residual = ~ dsum(~ id(units) | city),
                            data = training_data, 
                            na.action = na.method(x = "include", y = "include"), 
                            workspace = 6e8,
                            maxit = 100)
           
           # GEBV ordered by environment:
           n_loc = as.character(levels(training_data[["city"]]))
           n_gen = as.factor(levels(data[["code"]]))
           ngeno = length(n_gen); nsite = length(n_loc)
           
           blups = getting_eblups(training_data, model_fa, k, main_effect = FALSE)
           first_part = as.matrix(blups$coef.rr)
           second_part = blups$delta
           
           GEBV = as.matrix(first_part + second_part)
           GEBV = ge_vector(training_data, GEBV)
           
           # Predicted GEBV:
           results[[j]] = GEBV[GEBV[, "code"] %in% list.folder[[i]][[j]], ]
           
           # Results for each fold:
           print("fold")
           print(paste(j,"-------"))  
           print("rep")
           print(paste(i,"-------")) 
           
           }
        
           # GEBV and observed phenotypic data:
           pred_geno = ldply(results)
           pred_geno = pred_geno[order(pred_geno$city, as.numeric(pred_geno$code), decreasing = F),]
           colnames(pred_geno)[4] = "GEBV"
           
           geno_pheno = cbind(pred_geno, pheno) 
           
           # Predictive accuracy per environment:
           ac_env = do.call(rbind, lapply(split(geno_pheno, geno_pheno$city),  
                    function(x) data.frame(env = x$city[1], 
                    ac_env = cor(x$GEBV, x$grain_yield, use = "complete.obs"))))
        
           # Predictive accuracy:
           acuracy[[i]] = data.frame(C.Mourão_1 = ac_env["C. Mourão", "ac_env"], C.Mourão_2 = ac_env["C. Mourão (2)", "ac_env"], 
                                     Dourados = ac_env["Dourados", "ac_env"],
                                     Goiania = ac_env["Goiania", "ac_env"], Londrina_1 = ac_env["Londrina (1)", "ac_env"], 
                                     Londrina_2 = ac_env["Londrina (2)", "ac_env"], N.S_Dores_1 = ac_env["N. S. das Dores", "ac_env"], 
                                     N.S_Dores_2 = ac_env["N. S. das Dores ($)", "ac_env"], Planaltina = ac_env["Planaltina", "ac_env"], 
                                     S.R_Mangabeiras = ac_env["S. Rai. das Mangabeiras", "ac_env"], Sete.Lagoas = ac_env["Sete Lagoas", "ac_env"],
                                     Sinop_1 = ac_env["Sinop (1)", "ac_env"], 
                                     Sinop_2 = ac_env["Sinop (2)", "ac_env"], Vilhena = ac_env["Vilhena", "ac_env"])
           
      }
      
      # Results of cross validation:
      final_acuracy = ldply(acuracy)
      final_acuracy[rep+1, ] = colMeans(final_acuracy)
      final_acuracy[rep+2 ,] = apply(final_acuracy[-(rep+1) ,], 2, sd)  
      final_acuracy = round(final_acuracy, digits = 3)
      final_acuracy = data.frame(cbind(Rep = c(rep(1:rep), "Mean", "SD")), final_acuracy) 
      saveRDS(final_acuracy, file = "6_final_results/cv1_fa_5.rds")
      

      
      
  
      
      
      
      