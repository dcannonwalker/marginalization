mkdir w_wout_test1/data
rig run -f w_wout_test1/r/kparams.r
rig run -f w_wout_test1/r/simulation.r
rig run -f w_wout_test1/r/run_models.r
rig run -f w_wout_test1/r/run_models_moresamples.r
rig run -f w_wout_test1/r/save_summaries.r
rig run -f w_wout_test1/r/save_summaries_moresamples.r

