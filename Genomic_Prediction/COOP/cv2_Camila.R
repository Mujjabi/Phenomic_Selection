        ###########################
        #        Functions        #
        #     Jenifer Camila      #
        #          2024           #
        ###########################
        
        cv_1 = function(samp, k, rep){
           
             samp = split(samp$code, samp$group)
             samp = lapply(samp, droplevels)
             
             ng = sapply(samp, length)
           
             fl = lapply(ng, function(x){
           
             rep(floor(x/k), k) + c(rep(1, x%%k), rep(0, k - x%%k))
           
             })
          
             lapply(1:rep, function(rep){
           
             folds = mapply(function(x,y){
             
             split(x,sample(rep(1:k, y)))
            
             }, x = samp, y = fl, SIMPLIFY = F)
            
             Reduce(function(x,y){mapply(c, x,y)}, folds)
           
             })
        
        }
        
        
        cv_2 = function(samp, k, rep){
           
             samp = split(samp$code_city, samp$group)
             samp = lapply(samp, droplevels)
             
             ng = sapply(samp, length)
           
             fl = lapply(ng, function(x){
           
             rep(floor(x/k), k) + c(rep(1, x%%k), rep(0, k - x%%k))
           
             })
          
             lapply(1:rep, function(rep){
           
             folds = mapply(function(x,y){
             
             split(x, sample(rep(1:k, y)))
            
             }, x = samp, y = fl, SIMPLIFY = F)
            
             Reduce(function(x,y){mapply(c, x,y)}, folds)
           
             })
        
        }
        
        # Examples:
        
        # CV1
        #samp = data.frame(code = 1:200, group = rep(c(1,2), c(117, 83)))
        #folds1 = cv_1(samp, 5, 5)
        #str(folds1)
        
        # CV2
        #samp2 = data.frame(code = paste(rep(1:200, each = 4), rep(LETTERS[1:4], 200), sep = "."), group = rep(c(1,2), 4*c(117, 83)))
        #folds2 = cv_2(samp2, 5, 5)
        #str(folds2)
        
       
        
        
