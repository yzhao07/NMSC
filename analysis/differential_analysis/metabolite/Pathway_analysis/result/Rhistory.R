# PID of current job: 1747139
mSet<-InitDataObjects("conc", "pathqea", FALSE)
mSet<-Read.TextData(mSet, "Replacing_with_your_file_path", "rowu", "disc");
mSet<-SanityCheckData(mSet)
mSet<-ReplaceMin(mSet);
mSet<-CrossReferencing(mSet, "name");
mSet<-CreateMappingResultTable(mSet)
mSet<-PerformDetailMatch(mSet, "Caffeoylmalic acid.1");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "Caffeoylmalic acid.1", "Caffeoylmalic acid");
mSet<-PerformDetailMatch(mSet, "D-Lombricine.1");
mSet<-GetCandidateList(mSet);
mSet<-SetCandidate(mSet, "D-Lombricine.1", "D-Lombricine");
mSet<-PreparePrenormData(mSet)
mSet<-Normalization(mSet, "NULL", "NULL", "NULL", ratio=FALSE, ratioNum=20)
mSet<-PlotNormSummary(mSet, "norm_0_", "png", 72, width=NA)
mSet<-PlotSampleNormSummary(mSet, "snorm_0_", "png", 72, width=NA)
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateQeaScore(mSet, "rbc", "gt")
mSet<-PlotPathSummary(mSet, F, "path_view_0_", "png", 72, width=NA, NA, NA )
mSet<-SaveTransformedData(mSet)
mSet<-SetKEGG.PathLib(mSet, "hsa", "current")
mSet<-SetMetabolomeFilter(mSet, F);
mSet<-CalculateQeaScore(mSet, "rbc", "gt")
mSet<-PlotPathSummary(mSet, F, "path_view_1_", "png", 72, width=NA, NA, NA )
mSet<-PlotKEGGPath(mSet, "Citrate cycle (TCA cycle)",576, 480, "png", NULL)
mSet<-RerenderMetPAGraph(mSet, "zoom1699469399915.png",576.0, 480.0, 100.0)
mSet<-PlotKEGGPath(mSet, "Citrate cycle (TCA cycle)",576, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "Arginine biosynthesis",576, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "Citrate cycle (TCA cycle)",576, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "Nitrogen metabolism",576, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "Arginine biosynthesis",576, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "Alanine, aspartate and glutamate metabolism",576, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "Propanoate metabolism",576, 480, "png", NULL)
mSet<-PlotKEGGPath(mSet, "Propanoate metabolism",576, 480, "png", NULL)
mSet<-SaveTransformedData(mSet)
