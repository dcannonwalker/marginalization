# mkdir w_wout_test2/data
# rig run -f w_wout_test2/r/kparams.r
rig run -f w_wout_test2/r/simulation.r
rig run -f w_wout_test2/r/run_models.r
rig run -f w_wout_test2/r/run_models_moresamples.r
rig run -f w_wout_test2/r/save_summaries.r
rig run -f w_wout_test2/r/save_summaries_moresamples.r

