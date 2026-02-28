%% quickBatch.m — Run all 18 scenarios, report summary
setupTrackbench();

results = trackbench.batch.runAllScenarios("default");
