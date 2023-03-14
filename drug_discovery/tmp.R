
bkup <- av

tmp <- data.frame(mean_a = mean(av$`avExpr_B-cell`),
                  mean_b = mean(av$`avExpr_Solid tumor`),
                  sd_a = sd(av$`avExpr_B-cell`),
                  sd_b = sd(av$`avExpr_Solid tumor`),
                  n_a = length(av$`avExpr_B-cell`),
                  n_b = length(av$`avExpr_Solid tumor`))
sp <- sqrt( ((tmp$n_a-1)*tmp$sd_a^2 + (tmp$n_b-1)*tmp$sd_b^2) / (tmp$n_a+tmp$n_b-2))
tval <- (tmp$mean_a-tmp$mean_b) / sp*(sqrt(1/tmp$n_a+1/tmp$n_b))
pt(tval, (tmp$n_a+tmp$n_b-2), lower.tail = TRUE)
