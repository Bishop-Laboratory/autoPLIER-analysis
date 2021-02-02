# Converts PLIER output to language-agnostic format
load("analysis/PLIER_initial_testing/plierResult.rda")

dir.create("analysis/PLIER_initial_testing/plierResult", 
           showWarnings = FALSE)

write.csv(plierResult$residual, 
          file = "analysis/PLIER_initial_testing/plierResult/residual.csv")
write.csv(plierResult$B, 
          file = "analysis/PLIER_initial_testing/plierResult/B.csv")
write.csv(plierResult$Z, 
          file = "analysis/PLIER_initial_testing/plierResult/Z.csv")
write.csv(plierResult$U, 
          file = "analysis/PLIER_initial_testing/plierResult/U.csv")
write.csv(plierResult$C, 
          file = "analysis/PLIER_initial_testing/plierResult/C.csv")
write.csv(plierResult$summary, 
          file = "analysis/PLIER_initial_testing/plierResult/summary.csv")
write.csv(plierResult$Uauc, 
          file = "analysis/PLIER_initial_testing/plierResult/Uauc.csv")
write.csv(plierResult$Up, 
          file = "analysis/PLIER_initial_testing/plierResult/Up.csv")

