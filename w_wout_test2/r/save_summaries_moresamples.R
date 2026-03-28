fitw <- readRDS("w_wout_test1/data/fitw2000.rds")
fitwout <- readRDS("w_wout_test1/data/fitwout2000.rds")
smryw <- fitw$summary()
smrywout <- fitwout$summary()

saveRDS(smryw, "w_wout_test1/data/smryw2000.rds")
saveRDS(smrywout, "w_wout_test1/data/smrywout2000.rds")
