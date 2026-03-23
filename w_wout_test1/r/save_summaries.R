fitw <- readRDS("w_wout_test1/data/fitw.rds")
fitwout <- readRDS("w_wout_test1/data/fitwout.rds")
smryw <- fitw$summary()
smrywout <- fitwout$summary()

saveRDS(smryw, "w_wout_test1/data/smryw.rds")
saveRDS(smrywout, "w_wout_test1/data/smrywout.rds")
