#' Analysis
#' analysis of various experiment
#'
#' @slot project character. project name.
#' @slot author character. author name.
#' @slot analysis_name character. analysis name.
#' @slot experiment BulkExperiment. experiment need to be analyzed.
#' @slot output list. output, results put in analysis
#' @slot log character. any comment.
#'
#' @return Analysis
#' @export
#'
setClass(
  "Analysis",
  slots=c(
    project="character",
    author="character",
    analysis_name="character",
    experiment="BulkExperiment",
    output="list",
    log="character"
  ),
  prototype=list(
    project="myproject",
    author="",
    analysis_name="myanalysis",
    log=""
  )
)

#' runSurvivalAnalysis
#' run survival analysis, including univariate/multivariate survival analysis, as well as EIS
#'
#' @param obj Analysis. analysis object to be run.
#' @param survival_column character. column for survival days.
#' @param status_column character. column for survival status.
#' @param covariate_columns character. column for covariate passed to univariate,multivariate analysis.
#' @param extra_covariates character. column for extra covariate passed to univariate,multivariate analysis.
#' @param km_covariates character. column for additional km plot.
#' @param palette character. colors for KM curve.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runSurvivalAnalysis",function(obj,survival_column,status_column,covariate_columns,extra_covariates,palette,km_covariates=NULL,save_analysis=F,output_dir,log="run survival analysis,including univariant/multivariate analysis.") standardGeneric("runSurvivalAnalysis"))
setMethod("runSurvivalAnalysis","Analysis",function(obj,survival_column,status_column,covariate_columns,extra_covariates,palette,km_covariates=NULL,save_analysis=F,output_dir,log="run survival analysis,including univariant/multivariate analysis"){
  stopifnot("for survival analysis, patient info should be provided"=!utiltools::is.empty.data.frame(obj@experiment@patient_info))
  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  }
  #browser()
  obj@output[["univariate"]]<-list()
  obj@output[["multivariate"]]<-list()
  patient_info_<-merge(obj@experiment@sample_info,obj@experiment@patient_info,by="Patient_ID",all.x=T)
  stopifnot("Some variates are not in patient info!"=all(c(survival_column,status_column,covariate_columns) %in% colnames(patient_info_)))
  if(!is.null(km_covariates)){
    for(km_covariate in km_covariates){
      km_covariate_survival_data<-data.frame(Survival=patient_info_[[survival_column]],Status=patient_info_[[status_column]],Covariate=patient_info_[[km_covariate]])
      km_covariate_survival<-km_covariate_survival_data %>% survivalAnalysis::analyse_survival(dplyr::vars(Survival,Status),by=Covariate)
      if(missing(palette)){
        palette<-RColorBrewer::brewer.pal(length(unique(patient_info_[["EISG"]])),"Paired")[1:length(unique(km_covariate_survival_data[["Covariate"]]))]
      }

      km_covariate_km_plot<-survivalAnalysis::kaplan_meier_plot(km_covariate_survival,
                                                                break.time.by="breakByMonthYear",
                                                                xlab="Survival (month)",
                                                                legend.title=km_covariate,
                                                                hazard.ratio=TRUE,
                                                                risk.table=TRUE,
                                                                table.layout="clean",
                                                                ggtheme=ggplot2::theme_bw(10),
                                                                palette=palette,
                                                                legend=c(0.8,0.8))
      if(save_analysis){
        dir.create(file.path(output_dir,paste(km_covariate,"_survival",sep="")),recursive = T)
        svg(filename = file.path(output_dir,paste(km_covariate,"_survival",sep=""),paste(km_covariate,"_km_plot.svg",sep="")),width = 12,height=10)
        print(km_covariate_km_plot)
        dev.off()
      }
      obj@output[[paste(km_covariate,"_survival",sep="")]]<-km_covariate_km_plot
    }

  }


  survival_data<-data.frame(Survival=patient_info_[[survival_column]],Status=patient_info_[[status_column]],patient_info_[,covariate_columns])
  if(!missing(extra_covariates)){
    message("Please make sure extra_covariates are consistent with patient info!")
    survival_data<-cbind(survival_data,extra_variates)
  }
  covariates_<-setdiff(colnames(survival_data),c("Survival","Status"))
  #univariate
  univariate_results<-lapply(covariates_, function(covariate){
    res<-survivalAnalysis::analyse_multivariate(survival_data,dplyr::vars(Survival, Status),covariates = list(covariate))
    return(res)})
  univariate_forestplot<-survivalAnalysis::forest_plot(univariate_results,
                                                       endpoint_labeller = c("time"=time),
                                                       labels_displayed = c("endpoint", "factor", "n"),
                                                       orderer = ~order(factor.name),
                                                       ggtheme = ggplot2::theme_bw(base_size = 10),
                                                       relative_widths = c(1, 1.5, 1))
  univariate_results<-do.call(rbind,purrr::map(univariate_results,function(l) l[["summaryAsFrame"]]))
  obj@output$univariate[["univariate_results"]]<-univariate_results
  obj@output$univariate[["univariate_forestplot"]]<-univariate_forestplot

  #multivariate
  multivariate_results<-survivalAnalysis::analyse_multivariate(survival_data,dplyr::vars(Survival, Status),covariates = as.list(covariates_))
  multivariate_forestplot<-survivalAnalysis::forest_plot(multivariate_results,
                                                         endpoint_labeller = c("time"=time),
                                                         labels_displayed = c("endpoint", "factor", "n"),
                                                         orderer = ~order(factor.name),
                                                         ggtheme = ggplot2::theme_bw(base_size = 10),
                                                         relative_widths = c(1, 1.5, 1))
  multivariate_results<-multivariate_results[["summaryAsFrame"]]
  obj@output$multivariate[["multivariate_results"]]<-multivariate_results
  obj@output$multivariate[["multivariate_forestplot"]]<-multivariate_forestplot

  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  if(save_analysis){
    dir.create(file.path(output_dir,"univariate"),recursive = T)
    write.csv(obj@output$univariate[["univariate_results"]],file.path(output_dir,"univariate","univariate_results.csv"),row.names = F,quote = F)
    svg(filename = file.path(output_dir,"univariate","univariate_forestplot.svg"),width = 12,height=10)
    print(obj@output$univariate[["univariate_forestplot"]])
    dev.off()
    dir.create(file.path(output_dir,"multivariate"),recursive = T)
    write.csv(obj@output$multivariate[["multivariate_results"]],file.path(output_dir,"multivariate","multivariate_results.csv"),row.names = F,quote=F)
    svg(filename = file.path(output_dir,"multivariate","multivariate_forestplot.svg"),width = 12,height=10)
    print(obj@output$multivariate[["multivariate_forestplot"]])
    dev.off()
    cat(obj@log,file = file.path(output_dir,"log.txt"))
  }
  return(obj)
}
)


#' runVariantStatisticsAnalysis
#' run analysis about varaint statistics (mutation and copynumber)
#'
#' @param obj  Analysis. analysis object to be run.
#' @param cosmic_tier12_cancer_genes data.frame. cosmic tier1 and tier2 mutations, must have column "Gene Symbol".
#' @param genetic_cancer_genes  data.frame. cosmic tier1 and tier2 mutations, must have column "gene".
#' @param oncokb_cancer_genes  data.frame. cosmic tier1 and tier2 mutations, must have column "Hugo Symbol".
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runVariantStatisticsAnalysis",function(obj,cosmic_tier12_cancer_genes,genetic_cancer_genes,oncokb_cancer_genes,save_analysis=F,output_dir,log="run variant statistics analysis,including snp, indel and cnv.") standardGeneric("runVariantStatisticsAnalysis"))
setMethod("runVariantStatisticsAnalysis","Analysis",function(obj,cosmic_tier12_cancer_genes,genetic_cancer_genes,oncokb_cancer_genes,save_analysis=F,output_dir,log="run variant statistics analysis,including snp, indel and cnv."){
  stopifnot("snpindel should be provided!"=!utiltools::is.empty.data.frame(obj@experiment@snp_indel_assay@assay_data))
  stopifnot("Gistic2 directory should be provided"=dir.exists(obj@experiment@cnv_assay@Gistic2_directory))
  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  }
  DNA_samples<-obj@experiment@sample_info$Sample_ID[obj@experiment@sample_info$Assay=="DNA"]
  filtered_mutations<-obj@experiment@snp_indel_assay@assay_data
  filtered_cd_mutations<-filtered_mutations[!filtered_mutations$Variant_Classification %in% c("3'UTR","5'UTR"),]

  Gistic2_directory<-obj@experiment@cnv_assay@Gistic2_directory
  gistic2_segments_number_file<-file.path(Gistic2_directory,"sample_seg_counts.txt")
  gistic2_copynumber_file<-file.path(Gistic2_directory,"all_thresholded.by_genes.txt")

  segments_number<-read.table(gistic2_segments_number_file,header=T,sep="\t",stringsAsFactors = F,check.names = F)
  segments_number<-segments_number[match(DNA_samples,segments_number$sample),]
  copynumber<-read.table(gistic2_copynumber_file,header=T,sep="\t",stringsAsFactors = F,check.names = F)
  copynumber<-copynumber[,c("Gene Symbol","Locus ID","Cytoband",DNA_samples)]

  #table of samples about segments_number,copynumber_changed_gene_number, AMP_gene_number,DEL_gene_number,total_mutation_number,and various mutation number
  total_mutation_number=filtered_mutations %>% dplyr::group_by(Tumor_Sample_Barcode) %>% dplyr::summarise(total_mutation_number=dplyr::n())
  total_mutation_number<-total_mutation_number[match(DNA_samples,total_mutation_number$Tumor_Sample_Barcode),]
  specific_mutation_number<-filtered_mutations %>% dplyr::group_by(Tumor_Sample_Barcode,Variant_Classification) %>% dplyr::summarise(mutation_number=dplyr::n())
  specific_mutation_number<-tidyr::pivot_wider(specific_mutation_number,id_cols="Tumor_Sample_Barcode",names_from="Variant_Classification",values_from="mutation_number")
  specific_mutation_number<-specific_mutation_number[match(DNA_samples,specific_mutation_number$Tumor_Sample_Barcode),]
  assertthat::are_equal(DNA_samples,total_mutation_number$Tumor_Sample_Barcode)
  assertthat::are_equal(DNA_samples,specific_mutation_number$Tumor_Sample_Barcode)
  assertthat::are_equal(DNA_samples,colnames(copynumber)[-c(1,2,3)])
  assertthat::are_equal(DNA_samples,segments_number$sample)
  variant_statistic_table<-data.frame(Sample_ID=DNA_samples,
                                      Purity=obj@experiment@purity_assay@assay_data$cellularity[match(DNA_samples,obj@experiment@purity_assay@assay_data$Sample_ID)],
                                      Segment_Number=segments_number$segment_count,
                                      copynumberchanged_gene_number=colSums(copynumber[,-c(1,2,3)]!=0),
                                      AMP_gene_number=colSums(copynumber[,-c(1,2,3)]>0),
                                      DEL_gene_number=colSums(copynumber[,-c(1,2,3)]<0),
                                      specific_mutation_number[,-1],
                                      Total_mutation_number=total_mutation_number$total_mutation_number,check.names = F)
  variant_statistic_table$cd_mutation_number<-rowSums(specific_mutation_number[,setdiff(colnames(specific_mutation_number),c("Tumor_Sample_Barcode","3'UTR","5'UTR"))],na.rm=T)

  nonSilent_cosmic_mutation_number<-filtered_cd_mutations %>% dplyr::filter(Hugo_Symbol %in% cosmic_tier12_cancer_genes$`Gene Symbol`, Variant_Classification != "Silent") %>% dplyr::group_by(Tumor_Sample_Barcode) %>% dplyr::summarise(nonSilent_mutation_number=dplyr::n())
  variant_statistic_table$nonSilent_cd_cosmic_mutation_number<-nonSilent_cosmic_mutation_number$nonSilent_mutation_number[match(DNA_samples,nonSilent_cosmic_mutation_number$Tumor_Sample_Barcode)]
  nonSilent_genetic_mutation_number<-filtered_cd_mutations %>% dplyr::filter(Hugo_Symbol %in% genetic_cancer_genes$gene,Variant_Classification != "Silent") %>% dplyr::group_by(Tumor_Sample_Barcode) %>% dplyr::summarise(nonSilent_mutation_number=dplyr::n())
  variant_statistic_table$nonSilent_cd_genetic_mutation_number<-nonSilent_genetic_mutation_number$nonSilent_mutation_number[match(DNA_samples,nonSilent_genetic_mutation_number$Tumor_Sample_Barcode)]
  nonSilent_oncokb_mutation_number<-filtered_cd_mutations %>% dplyr::filter(Hugo_Symbol %in% oncokb_cancer_genes$`Hugo Symbol`,Variant_Classification != "Silent") %>% dplyr::group_by(Tumor_Sample_Barcode) %>% dplyr::summarise(nonSilent_mutation_number=dplyr::n())
  variant_statistic_table$nonSilent_cd_oncokb_mutation_number<-nonSilent_oncokb_mutation_number$nonSilent_mutation_number[match(DNA_samples,nonSilent_oncokb_mutation_number$Tumor_Sample_Barcode)]
  variant_statistic_table$nonSilent_cd_mutation_number<-variant_statistic_table$cd_mutation_number-variant_statistic_table$Silent
  obj@output$variant_statistics<-variant_statistic_table
  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  if(save_analysis){
    write.csv(obj@output$variant_statistics,file.path(output_dir,"variant_statistics.csv"),row.names = F,quote = F)
    cat(obj@log,file=file.path(output_dir,"log.txt"))
  }
  return(obj)
}
)


#' runCINTMBSurvivalAnalysis
#' run survival analysis about chromosome instability and tumor mutation burdon
#'
#' @param obj Analysis. analysis object to be run.
#' @param variantstatisticsanalysis Analysis. variantstatisticsanalysis object.
#' @param CIN_column character. column for chromosome instability.
#' @param TMB_column character. column for tumor mutation burdon.
#' @param cin_cutoff numeric. chromosome instability cutoff.
#' @param tmb_cutoff chromosome instability cutoff.
#' @param survival_column character. column for survival days.
#' @param status_column character. column for survival status.
#' @param palette character. colors for KM curve.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runCINTMBSurvivalAnalysis",function(obj,variantstatisticsanalysis,CIN_column,TMB_column,cin_cutoff,tmb_cutoff,survival_column,status_column,palette=c("High"="red","Low"="blue"),save_analysis=F,output_dir,log="run cin and tmb survival analysis.") standardGeneric("runCINTMBSurvivalAnalysis"))
setMethod("runCINTMBSurvivalAnalysis","Analysis",function(obj,variantstatisticsanalysis,CIN_column,TMB_column,cin_cutoff,tmb_cutoff,survival_column,status_column,palette=c("High"="red","Low"="blue"),save_analysis=F,output_dir,log="run cin and tmb survival analysis."){
  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  }

  message("TMB and CIN data come from variant statistics table, therefore please  run VariantStatistics analysis firstly!")
  cintmb_df<-data.frame(Sample_ID=variantstatisticsanalysis@output$variant_statistics[["Sample_ID"]],
                        CIN=variantstatisticsanalysis@output$variant_statistics[[CIN_column]],
                        TMB=variantstatisticsanalysis@output$variant_statistics[[TMB_column]])
  if(missing(cin_cutoff)){
    cintmb_df$CIN_Group=ifelse(variantstatisticsanalysis@output$variant_statistics[[CIN_column]]>median(variantstatisticsanalysis@output$variant_statistics[[CIN_column]],na.rm=T),"High","Low")
  } else {
    cintmb_df$CIN_Group=ifelse(variantstatisticsanalysis@output$variant_statistics[[CIN_column]]>cin_cutoff,"High","Low")
  }

  if(missing(tmb_cutoff)){
    cintmb_df$TMB_Group=ifelse(variantstatisticsanalysis@output$variant_statistics[[TMB_column]]>median(variantstatisticsanalysis@output$variant_statistics[[TMB_column]],na.rm=T),"High","Low")
  } else {
    cintmb_df$TMB_Group=ifelse(variantstatisticsanalysis@output$variant_statistics[[TMB_column]]>tmb_cutoff,"High","Low")
  }

  cintmb_df<-merge(cintmb_df,obj@experiment@sample_info[,c("Sample_ID","Patient_ID")],by="Sample_ID",all.x=T)
  patient_info<-data.frame(Patient_ID=obj@experiment@patient_info[["Patient_ID"]],
                           Survival=obj@experiment@patient_info[[survival_column]],
                           Status=obj@experiment@patient_info[[status_column]])
  cintmb_df<-merge(cintmb_df,patient_info,by='Patient_ID',all.x=T)
  #CIN survival analysis
  CIN_survival<-cintmb_df %>% survivalAnalysis::analyse_survival(dplyr::vars(Survival,Status),by=CIN_Group)
  CIN_km_plot<-survivalAnalysis::kaplan_meier_plot(CIN_survival,
                                 break.time.by="breakByMonthYear",
                                 xlab="Survival (month)",
                                 legend.title="CIN_Group",
                                 hazard.ratio=TRUE,
                                 risk.table=TRUE,
                                 table.layout="clean",
                                 ggtheme=ggplot2::theme_bw(10),
                                 palette=palette,
                                 legend=c(0.8,0.8))

  TMB_survival<-cintmb_df %>% survivalAnalysis::analyse_survival(dplyr::vars(Survival,Status),by=TMB_Group)
  TMB_km_plot<-survivalAnalysis::kaplan_meier_plot(TMB_survival,
                                 break.time.by="breakByMonthYear",
                                 xlab="Survival (month)",
                                 legend.title="CIN_Group",
                                 hazard.ratio=TRUE,
                                 risk.table=TRUE,
                                 table.layout="clean",
                                 ggtheme=ggplot2::theme_bw(10),
                                 palette=palette,
                                 legend=c(0.8,0.8))

  obj@output$cintmb_df=cintmb_df
  obj@output$CIN_survival=CIN_survival
  obj@output$CIN_km_plot=CIN_km_plot
  obj@output$TMB_survival=TMB_survival
  obj@output$TMB_km_plot=TMB_km_plot
  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  if(save_analysis) {
    write.csv(obj@output$cintmb_df,file.path(output_dir,"cintmb.csv"),row.names = F,quote = F)
    svg(filename = file.path(output_dir,"cin_km.svg"),width = 8,height=8)
    print(obj@output$CIN_km_plot)
    dev.off()
    svg(filename = file.path(output_dir,"tmb_km.svg"),width = 8,height=8)
    print(obj@output$TMB_km_plot)
    dev.off()
    cat(obj@log,file=file.path(output_dir,"log.txt"))
  }
  return(obj)
}
)



#' runCNVAnalysis
#' run analysis about copy number
#'
#' @param obj Analysis. analysis object to be run.
#' @param scale character. parameter passed to stableCluster, see manual of stableCluster.
#' @param method character. parameter passed to stableCluster,see manual of stableCluster.
#' @param subFeatureSize numeric. parameter passed to stableCluster,see manual of stableCluster.
#' @param subSampleSize numeric. parameter passed to stableCluster,see manual of stableCluster.
#' @param noise numeric. parameter passed to stableCluster,see manual of stableCluster.
#' @param cutFUN function.parameter passed to stableCluster,see manual of stableCluster.
#' @param nTimes numeric.parameter passed to stableCluster,see manual of stableCluster.
#' @param clusters numeric. parameter passed to stableCluster,see manual of stableCluster.
#' @param verbose logic. parameter passed to stableCluster,see manual of stableCluster.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param survival_column character. column for survival days.
#' @param status_column character. column for survival status.
#' @param palette character. colors for KM curve.
#' @param cosmic_tier12_cancer_genes data.frame. cosmic tier1 and tier2 mutations, must have column "Gene Symbol".
#' @param threshold numeric.
#' @param genome character. genome name, such as "hg19", "hg38".
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runCNVAnalysis",function(obj,
                                     scale=c("row","column","none"),method=c("bootstrap","perturb","combine"),subFeatureSize=0.8,subSampleSize=1,noise=1,cutFUN,nTimes=100,clusters=2,verbose=F,top_anno,top_anno_sample_order,sample_order=NULL,
                                     survival_column="OS",status_column="Status",palette=c("blue","green","red"),
                                     cosmic_tier12_cancer_genes, threshold=2,genome="hg19",
                                     save_analysis=F,output_dir,
                                     log="run cnv analysis.") standardGeneric("runCNVAnalysis"))
setMethod("runCNVAnalysis","Analysis",function(obj,
                                               scale=c("row","column","none"),method=c("bootstrap","perturb","combine"),subFeatureSize=0.8,subSampleSize=1,noise=1,cutFUN,nTimes=100,clusters=2,verbose=F,top_anno,top_anno_sample_order,sample_order=NULL,
                                               survival_column="OS",status_column="Status",palette=c("blue","green","red"),
                                               cosmic_tier12_cancer_genes, threshold=2,genome="hg19",
                                               save_analysis=F,output_dir,
                                               log="run cnv analysis."){
  scale=match.arg(scale)
  method=match.arg(method)
  stopifnot("Gistic2 directory should be provided"=dir.exists(obj@experiment@cnv_assay@Gistic2_directory))
  stopifnot("Mutsig2CV directory should be provided"=dir.exists(obj@experiment@snp_indel_assay@Mutsig2CV_directory))

  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  }

  Gistic2_directory<-obj@experiment@cnv_assay@Gistic2_directory
  Mutsig2CV_directory<-obj@experiment@snp_indel_assay@Mutsig2CV_directory
  #gistic2_segments_number_file<-file.path(Gistic2_directory,"sample_seg_counts.txt")
  #gistic2_copynumber_file<-file.path(Gistic2_directory,"all_thresholded.by_genes.txt")
  gistic2_all_lesions_file=file.path(Gistic2_directory,"all_lesions.conf_99.txt")
  gistic2_scores_file<-file.path(Gistic2_directory,"scores.gistic")
  gistic2_arm_file=file.path(Gistic2_directory,"broad_values_by_arm.txt")
  gistic2_arm_significance_file<-file.path(Gistic2_directory,"broad_significance_results.txt")
  mutsig2cv_sig_genes_file<-file.path(Mutsig2CV_directory,"sig_genes.txt")

  #DNA_samples<-obj@experiment@sample_info$Sample_ID[obj@experiment@sample_info$Assay=="DNA"]
  DNA_samples<-top_anno_sample_order

  segments<-obj@experiment@cnv_assay@assay_data
  #segments_number<-read.table(gistic2_segments_number_file,header=T,sep="\t",stringsAsFactors = F,check.names = F)
  #segments_number<-segments_number[match(DNA_samples,segments_number$sample),]
  #copynumber<-read.table(gistic2_copynumber_file,header=T,sep="\t",stringsAsFactors = F,check.names = F)
  #copynumber<-copynumber[,c("Gene Symbol","Locus ID","Cytoband",DNA_samples)]
  all_lesions<-read.table(gistic2_all_lesions_file,header=T,sep="\t",stringsAsFactors = F,check.names = F)
  lesions<-all_lesions[(nrow(all_lesions)/2+1):nrow(all_lesions),-c(1:9,ncol(all_lesions))]
  rownames(lesions)<-paste(all_lesions[(nrow(all_lesions)/2+1):nrow(all_lesions),"Descriptor"], gsub(".*\\:","",gsub("\\(.*","",all_lesions[(nrow(all_lesions)/2+1):nrow(all_lesions),"Wide Peak Limits"])),sep=":")
  lesions<-lesions[!duplicated(lesions),DNA_samples]
  scores<-read.table(gistic2_scores_file,header=T,sep="\t",stringsAsFactors = F,check.names=F)
  colnames(scores)[5]<-"neg_log10_q"
  arms<-read.table(gistic2_arm_file,header=T,sep="\t",stringsAsFactors = F,check.names = F)
  colnames(arms)[1]<-"arm"
  arms<-arms[,c("arm",DNA_samples)]
  arm_significance<-read.csv(gistic2_arm_significance_file,header=T,sep="\t",stringsAsFactors = F,check.names = F)

  #unsupervised lesion heatmap using significanty changed lesions by GISTIC2
  message("Please make sure that annotation should have same sample order as cluster!")
  assertthat::are_equal(colnames(lesions),DNA_samples)
  if(missing(clusters)){
    clusters<-utiltools::estimate_bestNumberofClusters(lesions)$Best.NumberofCluster
  }

  stableclusters_lesions<-utiltools::stableCluster(lesions,scale=scale,method=method,subFeatureSize=subFeatureSize,subSampleSize=subSampleSize,noise=noise,cutFUN=cutFUN,nTimes=nTimes,clusters=clusters,verbose=verbose)
  assertthat::are_equal(colnames(stableclusters_lesions),colnames(lesions))
  stablecluster_lesions_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(lesions))),name="zscore_lesion",cluster_columns = hclust(dist(stableclusters_lesions)),top_annotation = top_anno,column_split = clusters,show_row_names = F)
  if(is.null(sample_order)){
    sample_order=colnames(lesions)[unlist(ComplexHeatmap::column_order(stablecluster_lesions_heatmap))]
  }

  #cluster survival anaysis
  patient_info_<-data.frame(Sample_ID=DNA_samples)
  patient_info_<-merge(patient_info_,obj@experiment@sample_info[,c("Sample_ID","Patient_ID")],by="Sample_ID",all.x=T)
  patient_info_<-merge(patient_info_,obj@experiment@patient_info,by="Patient_ID",all.x=T)
  patient_info_<-patient_info_[match(DNA_samples,patient_info_[["Sample_ID"]]),]
  hCluster_km_plot<-utiltools::hCluster_Surv(heatmap=stablecluster_lesions_heatmap,clinics=patient_info_,survival_col = survival_column,status_col = status_column,palette = palette[1:clusters])

  fb_bin_segments<-lapply(unique(segments[[obj@experiment@cnv_assay@sample_id_column]]),function(samp){data=segments[segments[[obj@experiment@cnv_assay@sample_id_column]]==samp,];
  gr=utiltools::df2granges(df=data,genome=genome,seqlevelsStyle = "NCBI",meta_cols = "segmean");
  bingenome<-utiltools::bingranges(gr);
  return(bingenome$segmean)})
  fb_bin_segments<-do.call(cbind,fb_bin_segments)
  colnames(fb_bin_segments)<-unique(segments[[obj@experiment@cnv_assay@sample_id_column]])

  fb_bin_segments<-fb_bin_segments[,DNA_samples]
  bingenome<-utiltools::bingranges(utiltools::df2granges(df=segments[segments[[obj@experiment@cnv_assay@sample_id_column]]==segments[[obj@experiment@cnv_assay@sample_id_column]][1],],xy=F,genome=genome,seqlevelsStyle = "NCBI",meta_cols = "segmean"))
  seq_df<-data.frame(seq=as.character(GenomeInfoDb::seqnames(bingenome)),seq_name="")
  for(chr in paste("chr",as.character(1:22),sep="")){ seq_df$seq_name[floor((min(which(seq_df$seq==chr))+max(which(seq_df$seq==chr)))/2)]<-gsub("chr","",chr)}
  seq_colors<-rep(c("grey","black"),11)
  names(seq_colors)<-paste("chr",as.character(1:22),sep="")
  ha<-ComplexHeatmap::rowAnnotation(text=ComplexHeatmap::anno_text(seq_df$seq_name,gp=grid::gpar(fontsize=8)),seq=seq_df$seq,show_legend=F,col=list(seq=seq_colors),show_annotation_name=F)
  col_fun<-circlize::colorRamp2(c(-1.5,-0.3,0,0.3,1.5),c("blue","white","white","white","red"))
  cn_heatmap<-ComplexHeatmap::draw(ComplexHeatmap::Heatmap(fb_bin_segments,name="copynumber",col=col_fun,cluster_rows = F,column_order = sample_order,column_dend_reorder = T,top_annotation = top_anno,left_annotation = ha))
  if(!file.exists(mutsig2cv_sig_genes_file)){cat("Please run mutsigcv2 first")} else{
    sig_genes<-read.csv(mutsig2cv_sig_genes_file,header=T,stringsAsFactors = F,check.names = F,sep="\t")
  }

  mutations<-obj@experiment@snp_indel_assay@assay_data
  nonSilent_mutations<-mutations[(!mutations$Variant_Classification %in% c("Silent")) & (mutations$Tumor_Sample_Barcode %in% DNA_samples), ]
  nonSilent_sig_mutations<-nonSilent_mutations[nonSilent_mutations$Hugo_Symbol %in% sig_genes$gene[sig_genes$q<=0.05 & sig_genes$nnon>=2],]


  cancer_gene_info<-data.frame(Symbol=cosmic_tier12_cancer_genes$`Gene Symbol`,
                               Chromosome=gsub(":.*","",cosmic_tier12_cancer_genes$`Genome Location`),
                               Start=as.numeric(gsub("-.*","",gsub(".*:","",cosmic_tier12_cancer_genes$`Genome Location`))),
                               End=as.numeric(gsub(".*-","",gsub(".*:","",cosmic_tier12_cancer_genes$`Genome Location`))))
  cancer_gene_info<-cancer_gene_info[complete.cases(cancer_gene_info),]
  genome_plot<-utiltools::gggenome(segments,sample_col=obj@experiment@cnv_assay@sample_id_column,segments_anno_df = nonSilent_sig_mutations,sample_order = sample_order,scores_df = scores,threshold = threshold,scores_anno_df  = cancer_gene_info,anno_score_text_size  = 2,genome = genome)

  arm_plot<-utiltools::ggarms(arms = arms,arm_significance = arm_significance,sample_order =  sample_order,threshold = threshold,genome = genome,seqlevelsStyle = "NCBI",simplified = T,xy = F)

  obj@output$lesion_heatmap<-list(lesions=lesions,stableclusters_lesions=stableclusters_lesions,stablecluster_lesions_heatmap=stablecluster_lesions_heatmap,hCluster_km_plot=hCluster_km_plot)
  obj@output$cn_heatmap<-cn_heatmap
  obj@output$genome_plot<-list(segments=segments,nonSilent_sig_mutations=nonSilent_sig_mutations,scores=scores,genome_plot=genome_plot)
  obj@output$arm_plot<-list(arms=arms,arm_significance=arm_significance,arm_plot=arm_plot)
  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  if(save_analysis) {
    if(!dir.exists(file.path(output_dir,"lesion_heatmap"))){dir.create(file.path(output_dir,"lesion_heatmap"),recursive = T)}
    write.csv(lesions,file=file.path(output_dir,"lesion_heatmap","lesions.csv"),row.names = F,quote = F)
    write.csv(stableclusters_lesions,file=file.path(output_dir,"lesion_heatmap","stableclusters_lesions.csv"),row.names = F,quote = F)
    svg(filename = file.path(output_dir,"lesion_heatmap","stablecluster_lesions_heatmap.svg"),width = 8,height=8)
    print(stablecluster_lesions_heatmap)
    dev.off()
    svg(filename = file.path(output_dir,"lesion_heatmap","hCluster_km.svg"),width = 8,height=8)
    print(hCluster_km_plot)
    dev.off()

    svg(filename = file.path(output_dir,"cn_heatmap.svg"),width = 8,height=8)
    print(cn_heatmap)
    dev.off()

    if(!dir.exists(file.path(output_dir,"genomeplot"))){dir.create(file.path(output_dir,"genomeplot"),recursive = T)}
    write.csv(segments,file=file.path(output_dir,"genomeplot","segments.csv"),row.names = F,quote = F)
    write.csv(nonSilent_sig_mutations,file=file.path(output_dir,"genomeplot","nonSilent_sig_mutations.csv"),row.names = F,quote = F)
    write.csv(scores,file=file.path(output_dir,"genomeplot","scores.csv"),row.names = F,quote = F)
    svg(filename = file.path(output_dir,"genomeplot","genome_plot.svg"),width = 12,height=8)
    print(genome_plot)
    dev.off()

    if(!dir.exists(file.path(output_dir,"armplot"))){dir.create(file.path(output_dir,"armplot"),recursive = T)}
    write.csv(arms,file=file.path(output_dir,"armplot","arms.csv"),row.names = F,quote = F)
    write.csv(arm_significance,file=file.path(output_dir,"armplot","arm_significance.csv"),row.names = F,quote = F)
    svg(filename = file.path(output_dir,"armplot","arm_plot.svg"),width = 12,height=8)
    print(arm_plot)
    dev.off()

    cat(obj@log,file=file.path(output_dir,"log.txt"))
  }
  return(obj)
}
)


#' runMutationalSigatureAnalysis
#'
#' run analysis about mutational signature (snp)
#'
#' @param obj Analysis. analysis object to be run.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param extension numeric.parameter passed mut_matrix.beside mutation site, number of base pair recruited for mutation signature(pattern) analysis.
#' @param profile_group character. comparison groups or column. parameter passed to plot_profile_heatmap.
#' @param rank numeric. parameter passed to nmf.rank of dimensions.
#' @param nrun numeric. parameter passed to nmf.number of runs.
#' @param cutoff numeric. parameter passed to rename_nmf_signatures.cutoff of correlation with known signatures.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runMutationalSigatureAnalysis",function(obj,top_anno,top_anno_sample_order,extension=2,profile_group,rank = 3, nrun = 10,cutoff = 0.85,save_analysis=F,output_dir,log="run mutational signature analysis.") standardGeneric("runMutationalSigatureAnalysis"))
setMethod("runMutationalSigatureAnalysis","Analysis",function(obj,top_anno,top_anno_sample_order,extension=2,profile_group,rank = 3, nrun = 10,cutoff = 0.85,save_analysis=F,output_dir,log="run mutational signature analysis."){
  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  }
  DNA_samples<-top_anno_sample_order
  mutations=obj@experiment@snp_indel_assay@assay_data
  mutations<-mutations[mutations$Variant_Type=="SNP",]
  mutations<-mutations[mutations$Tumor_Sample_Barcode %in% DNA_samples,]
  mutations<-mutations[!is.na(mutations$Variant_Classification),]
  colnames(mutations)[match(c("Reference_Allele","Tumor_Seq_Allele2"),colnames(mutations))]<-c("REF","ALT")
  sample_info<-merge(obj@experiment@sample_info,obj@experiment@patient_info,by="Patient_ID",all.x=T)
  sample_info<-sample_info[match(DNA_samples,sample_info[["Sample_ID"]]),]
  if(length(profile_group)==1){
    profile_group_<-sample_info[[profile_group]]
    if(any(is.na(profile_group_) | profile_group_=="")){
      excluded_samples<-sample_info[["Sample_ID"]][is.na(profile_group_) | profile_group_==""]
      sample_info=sample_info[-match(excluded_samples,sample_info[["Sample_ID"]]),]
      DNA_samples=DNA_samples[-match(excluded_samples,DNA_samples)]
      top_anno<-top_anno[match(DNA_samples,top_anno_sample_order)]
      assertthat::are_equal(sample_info[["Sample_ID"]],DNA_samples)
      profile_group<-sample_info[[profile_group]]
      mutations=mutations[!mutations$Tumor_Sample_Barcode %in% excluded_samples,]
    } else {
      profile_group<-profile_group_
    }
  }
  mutations$Tumor_Sample_Barcode<-factor(mutations$Tumor_Sample_Barcode,levels=DNA_samples)

  mutations_grl<- utiltools::df2grangelist(mutations,genome="hg19",seqlevelsStyle = "NCBI",simplified = F,xy=T,sample_col = "Tumor_Sample_Barcode",seqnames_col = "Chromosome",start_col = "Start_position",end_col = "End_position",strand_col = "Strand",meta_cols = c("REF","ALT","Hugo_Symbol"))
  genome<-BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
  snv_grl <- MutationalPatterns::get_mut_type(mutations_grl, type = "snv")
  type_occurrences <- MutationalPatterns::mut_type_occurrences(snv_grl, genome)
  spectrum_plot<-MutationalPatterns::plot_spectrum(type_occurrences, CT = TRUE, indv_points = TRUE, legend = T)

  mut_mat_ext_context <- MutationalPatterns::mut_matrix(snv_grl, genome, extension = extension)[,DNA_samples]
  assertthat::are_equal(colnames(mut_mat_ext_context),DNA_samples)

  profile_heatmap<-MutationalPatterns::plot_profile_heatmap(mut_mat_ext_context,by=profile_group)


  mut_mat <- MutationalPatterns::mut_matrix(vcf_list = snv_grl, ref_genome = genome)[,DNA_samples]
  if(any(rowSums(mut_mat)==0)){
    mut_mat[unname(which(rowSums(mut_mat)==0)),sample(1:ncol(mut_mat),1)]<-1
  }
  require(NMF)
  estimate <- NMF::nmf(mut_mat, rank = 2:5, method = "brunet", nrun = 10, seed = 123456, .opt = "v-p")
  estimate_plot<-plot(estimate)


  mut_mat_ <- mut_mat + 0.0001
  nmf_res <- MutationalPatterns::extract_signatures(mut_mat_, rank = rank, nrun = nrun, single_core = TRUE)
  signatures = MutationalPatterns::get_known_signatures()
  nmf_res <- MutationalPatterns::rename_nmf_signatures(nmf_res, signatures, cutoff = cutoff)
  profile96_plot<-MutationalPatterns::plot_96_profile(nmf_res$signatures, condensed = TRUE)
  contribution_barplot<-MutationalPatterns::plot_contribution(nmf_res$contribution, nmf_res$signature,mode = "relative")
  contribution_heatmap<-MutationalPatterns::plot_contribution_heatmap(nmf_res$contribution, cluster_samples = TRUE,cluster_sigs = TRUE)


  rel_contribution<-t(t(nmf_res$contribution)/colSums(nmf_res$contribution))[,DNA_samples]
  col=circlize::colorRamp2(c(0,1),c("white","blue"))
  assertthat::are_equal(colnames(rel_contribution),DNA_samples)
  rel_heatmap<-ComplexHeatmap::Heatmap(rel_contribution,name="relative_contribution",col=col,top_annotation = top_anno,height=grid::unit(4,'cm'))

  cos_sim_samples_signatures <- t(MutationalPatterns::cos_sim_matrix(mut_mat, signatures))[,DNA_samples]
  cos_sim_samples_signature_heatmap<-ComplexHeatmap::Heatmap(cos_sim_samples_signatures,col=col,top_annotation = top_anno)

  obj@output$snv=snv_grl
  obj@output$type_occurrences<-type_occurrences
  obj@output$spectrum_plot=spectrum_plot
  obj@output$mut_mat_ext_context=mut_mat_ext_context
  obj@output$profile_heatmap=profile_heatmap
  obj@output$mut_mat=mut_mat
  obj@output$estimate=estimate
  obj@output$estimate_plot=estimate_plot
  obj@output$nmf_signature<-nmf_res
  obj@output$profile96_plot<-profile96_plot
  obj@output$contribution_barplot<-contribution_barplot
  obj@output$contribution_heatmap<-contribution_heatmap
  obj@output$relative_contribution<-rel_contribution
  obj@output$relative_contribution_heatmap<-rel_heatmap
  obj@output$cos_sim_samples_signatures<-cos_sim_samples_signatures
  obj@output$cos_sim_samples_signature_heatmap<-cos_sim_samples_signature_heatmap
  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  if(save_analysis){
    write.csv(type_occurrences,file.path(output_dir,"type_occurrences.csv"),row.names = F,quote = F)
    write.csv(mut_mat_ext_context,file.path(output_dir,"mut_mat_ext_context.csv"),row.names = F,quote = F)
    write.csv(mut_mat,file.path(output_dir,"mut_mat.csv"),row.names = F,quote = F)
    write.csv(rel_contribution,file.path(output_dir,"rel_contribution.csv"),row.names = F,quote = F)
    write.csv(cos_sim_samples_signatures,file.path(output_dir,"cos_sim_samples_signatures.csv"),row.names = F,quote = F)
    svg(filename = file.path(output_dir,"spectrum_plot.svg"),width = 8,height=8)
    print(spectrum_plot)
    dev.off()
    svg(filename = file.path(output_dir,"profile_heatmap.svg"),width = 8,height=8)
    print(profile_heatmap)
    dev.off()
    svg(filename = file.path(output_dir,"estimate_plot.svg"),width = 8,height=8)
    print(estimate_plot)
    dev.off()
    svg(filename = file.path(output_dir,"profile96_plot.svg"),width = 8,height=8)
    print(profile96_plot)
    dev.off()
    svg(filename = file.path(output_dir,"contribution_barplot.svg"),width = 8,height=8)
    print(contribution_barplot)
    dev.off()
    svg(filename = file.path(output_dir,"contribution_heatmap.svg"),width = 8,height=8)
    print(contribution_heatmap)
    dev.off()
    svg(filename = file.path(output_dir,"rel_heatmap.svg"),width = 8,height=8)
    print(rel_heatmap)
    dev.off()
    svg(filename = file.path(output_dir,"cos_sim_samples_signature_heatmap.svg"),width = 8,height=15)
    print(cos_sim_samples_signature_heatmap)
    dev.off()

    cat(obj@log,file=file.path(output_dir,"log.txt"))
  }
  return(obj)
}
)



#' runMutationalMafAnalysis
#'
#' run analysis about mutation using maftools
#'
#' @param obj Analysis. analysis object to be run.
#' @param maftools_mutation_conversion_dict data.frame. dictionary of cgl mutation to maftools compatible mutation.
#' @param cancer_genes character. cancer gene list.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param sample_order character. oncoplot column order.
#' @param pathways list. pathway list containing pathway name and corresponding genes.
#' @param group_column character. group column needed compare two diffenrent cohorts
#' @param group_terms character. group terms needed for comparisons.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#' @param ... list. parameters passed to oncoplot.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runMutationalMafAnalysis",function(obj,maftools_mutation_conversion_dict=NULL,cancer_genes,top_anno,top_anno_sample_order,sample_order=NULL,pathways,group_column=NULL,group_terms=NULL,save_analysis=F,output_dir,log="run mutational maf analysis.",...) standardGeneric("runMutationalMafAnalysis"))
setMethod("runMutationalMafAnalysis","Analysis",function(obj,maftools_mutation_conversion_dict,cancer_genes,top_anno,top_anno_sample_order,sample_order=NULL,pathways,group_column=NULL,group_terms=NULL,save_analysis=F,output_dir,log="run mutational maf analysis.",...){

  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  } else {
    output_dir<-getwd()
  }
  DNA_samples<-top_anno_sample_order
  sample_info<-merge(merge(data.frame(Sample_ID=DNA_samples),obj@experiment@sample_info[,c("Sample_ID","Patient_ID")],by="Sample_ID",all.x=T),obj@experiment@patient_info,by="Patient_ID",all.x=T)
  sample_info$Tumor_Sample_Barcode<-sample_info$Sample_ID

  write.csv(sample_info,file=file.path(output_dir,"sample_info.tsv"),row.names = F)
  mutations=obj@experiment@snp_indel_assay@assay_data
  if(!is.null(maftools_mutation_conversion_dict)){
    mutations$Variant_Classification<-unname(maftools_mutation_conversion_dict[mutations$Variant_Classification])
    mutations<-mutations[mutations$Variant_Classification %in% maftools_mutation_conversion_dict,]
  }
  mutations<-mutations[!is.na(mutations$Variant_Classification),]
  mutations<-mutations[mutations$Tumor_Sample_Barcode %in% DNA_samples,]
  mutations$Tumor_Sample_Barcode<-factor(mutations$Tumor_Sample_Barcode,levels=DNA_samples)
  write.table(mutations,file=file.path(output_dir,"mutations.maf"),row.names = F,col.names = T,sep="\t",quote = F)
  maf = maftools::read.maf(maf = file.path(output_dir,"mutations.maf"),
                 clinicalData = file.path(output_dir,"sample_info.tsv"),
                 verbose = FALSE)
  if(!is.null(sample_order)){
    top_anno<-top_anno[match(sample_order,DNA_samples),]
    DNA_samples<-sample_order
    mutations$Tumor_Sample_Barcode<-factor(mutations$Tumor_Sample_Barcode,levels=sample_order)
    sample_info<-sample_info[match(DNA_samples,sample_info$Tumor_Sample_Barcode),]
  }


  #oncosummary
  svg(filename = file.path(output_dir,"oncosummary.svg"),width = 8,height=8)
  maftools::plotmafSummary(maf = maf, rmOutlier = TRUE, addStat = 'median', dashboard = TRUE, titvRaw = FALSE)
  obj@output$plotmafSummary=file.path(output_dir,"oncosummary.svg")
  dev.off()

  Mutsig2CV_directory<-obj@experiment@snp_indel_assay@Mutsig2CV_directory
  mutsig2cv_sig_genes_file<-file.path(Mutsig2CV_directory,"sig_genes.txt")
  if(!file.exists(mutsig2cv_sig_genes_file)){cat("Please run mutsigcv2 first")} else{
    sig_genes<-read.csv(mutsig2cv_sig_genes_file,header=T,stringsAsFactors = F,check.names = F,sep="\t")
  }
  selected_genes<-sig_genes$gene[sig_genes$p<=0.05 & sig_genes$nnon>=ifelse(length(unique(mutations$Tumor_Sample_Barcode))>20,2,1)]
  if(!missing(cancer_genes)){
    selected_genes<-selected_genes[selected_genes %in% cancer_genes]
  }
  #oncoplot
  svg(filename = file.path(output_dir,"oncoplot.svg"),width = 8,height=8)
  onco_plot<-maftools::oncoplot(maf,sortByAnnotation = T,genes = selected_genes,removeNonMutated = F,...)
  obj@output$onco_plot=file.path(output_dir,"oncoplot.svg")
  dev.off()

  snp_indel_cnv_col<-c("Missense_Mutation"="#33A02B","Nonsense_Mutation"="#ff0000","In_Frame_Del"="#FEFF99","In_Frame_Ins"="#938038","Frame_Shift_Del"="#1F78B4","Frame_Shift_Ins"="#7b42f5","Splice_Site"="#00ffd5","Nonstop_Mutation"="#0022fc","UTR3"="#707173","UTR5"="#464647","Multi_Hit"="#000000")
  complex_onco_plot<-utiltools::complex_oncoplot(snp_indels=mutations,selected_genes = selected_genes,mutsigCV2_sig_genes = sig_genes,multi_hit=TRUE,cnv_types = c("Loss","Gain","AMP","Del") ,show_heatmap_legend = T,cnv=F,top_filter = 100L,remove_macromolecular_gene = T,col=snp_indel_cnv_col,macromolecular_threshold = 300000,top_annotation = top_anno,sample_order =sample_order,output_dir=output_dir)
  #plot titv summary
  sample_titv = maftools::titv(maf = maf, plot = FALSE, useSyn = TRUE)
  svg(filename = file.path(output_dir,"titv.svg"),width = 8,height=8)
  titv_plot<-maftools::plotTiTv(res = sample_titv)
  obj@output$titv_plot=file.path(output_dir,"titv.svg")
  dev.off()

  svg(filename = file.path(output_dir,"mutload.svg"),width = 8,height=8)
  mutload<-maftools::tcgaCompare(maf = maf, cohortName = '*COHORT', logscale = TRUE, capture_size = 35.8) #no UTRs, non_coding exons
  obj@output$mutload=file.path(output_dir,"mutload.svg")
  dev.off()

  svg(filename = file.path(output_dir,"vaf.svg"),width = 8,height=8)
  maftools::plotVaf(maf = maf, vafCol = 'i_TumorVAF_WU')
  obj@output$plotVaf=file.path(output_dir,"vaf.svg")
  dev.off()

  svg(filename = file.path(output_dir,"mutationcocorrence.svg"),width = 8,height=8)
  somaticinteractions<-maftools::somaticInteractions(maf = maf, top = 25, pvalue = c(0.05, 0.1))
  obj@output$somaticinteractions=file.path(output_dir,"mutationcocorrence.svg")
  dev.off()
  #browser()
  svg(filename = file.path(output_dir,"pfamdomain.svg"),width = 8,height=8)
  pfam = maftools::pfamDomains(maf = maf, AACol = 'Protein_Change', top = 10)
  obj@output$pfam=pfam
  obj@output$pfamplot=file.path(output_dir,"pfamdomain.svg")
  dev.off()

  pathdb <- system.file("extdata", "oncogenic_sig_patwhays.tsv", package = "maftools")
  pathdb = data.table::fread(input = pathdb)
  pathdb_ = split(pathdb$Gene, as.factor(pathdb$Pathway))
  pathdb = split(pathdb, as.factor(pathdb$Pathway))
  if(missing(pathways)){
    pathways=pathdb_
  }

  f<-function(l,multi_hit=TRUE) {
    if(length(l)==1){
      return(as.character(l))
    } else{
      if(multi_hit){
        return ("Multi_Hit")
      } else{
        return(paste(l,collapse=";"))
      }}}
  mutations_<-tidyr::pivot_wider(mutations,id_cols = "Hugo_Symbol",names_from="Tumor_Sample_Barcode",values_from = "Variant_Classification",values_fn = f)
  mutations_<-as.data.frame(mutations_)
  mutations_<-mutations_[,c("Hugo_Symbol",DNA_samples)]
  pathways_<-lapply(pathdb,function(pathway){return(pathway$Gene)})
  oncogenicpathwaytableplot<-utiltools::oncogenicpathway_tableplot(pathways = pathways_,mutations = mutations_,sample_info = sample_info,gene_column = "Hugo_Symbol",sample_id_column = "Sample_ID",top_anno = top_anno,top_anno_sample_order = DNA_samples,sample_order = DNA_samples,show_mutations=T)

  obj@output$mutations<-mutations_
  obj@output$complex_onco_plot<-complex_onco_plot
  obj@output$concogenicpathwaytableplot<-oncogenicpathwaytableplot
  #obj@output$pathwaymuttableplot<-pathwaymuttableplot
  obj@output$maf<-maf

  if(!is.null(group_column)){
    if(is.null(group_terms)){group_terms<-setdiff(na.omit(sample_info[[group_column]]),"")}
    cohort1_term=group_terms[1]
    cohort1_sample_info<-sample_info[sample_info[[group_column]]==cohort1_term,]
    cohort1_samples<-unique(cohort1_sample_info$Sample_ID)
    cohort1_mutations<-mutations[mutations$Tumor_Sample_Barcode %in% cohort1_samples,]
    cohort1_mutations$Tumor_Sample_Barcode<-factor(cohort1_mutations$Tumor_Sample_Barcode,levels=cohort1_samples)
    write.table(cohort1_mutations,file=file.path(tempdir(),"cohort1_mutations.maf"),row.names = F,col.names = T,sep="\t",quote = F)
    cohort1_maf = maftools::read.maf(maf = file.path(tempdir(),"cohort1_mutations.maf"))
    file.remove(file.path(tempdir(),"cohort1_mutations.maf"))
    cohort2_term=group_terms[2]
    cohort2_sample_info<-sample_info[sample_info[[group_column]]==cohort2_term,]
    cohort2_samples<-unique(cohort2_sample_info$Sample_ID)
    cohort2_mutations<-mutations[mutations$Tumor_Sample_Barcode %in% cohort2_samples,]
    cohort2_mutations$Tumor_Sample_Barcode<-factor(cohort2_mutations$Tumor_Sample_Barcode,levels=cohort2_samples)
    write.table(cohort2_mutations,file=file.path(tempdir(),"cohort2_mutations.maf"),row.names = F,col.names = T,sep="\t",quote = F)
    cohort2_maf = maftools::read.maf(maf = file.path(tempdir(),"cohort2_mutations.maf"))
    file.remove(file.path(tempdir(),"cohort2_mutations.maf"))
    cohort1_hits<-cohort1_mutations %>% dplyr::filter(Variant_Classification!="Silent") %>% dplyr::group_by(Hugo_Symbol) %>% dplyr::summarise(cohort1_hits=dplyr::n_distinct(Tumor_Sample_Barcode))
    cohort2_hits<-cohort2_mutations %>% dplyr::filter(Variant_Classification!="Silent") %>% dplyr::group_by(Hugo_Symbol) %>% dplyr::summarise(cohort2_hits=dplyr::n_distinct(Tumor_Sample_Barcode))
    OR_table<-merge(cohort1_hits,cohort2_hits,by="Hugo_Symbol",all=T)
    OR_table[is.na(OR_table)]<-0
    OR_table<-as.data.frame(OR_table)
    OR_table$cohort1_nonhits<-length(unique(cohort1_mutations$Tumor_Sample_Barcode))-OR_table$cohort1_hits
    OR_table$cohort2_nonhits<-length(unique(cohort2_mutations$Tumor_Sample_Barcode))-OR_table$cohort2_hits
    pvalues<-apply(OR_table,1,function(row_data){mat<-matrix(as.numeric(c(row_data["cohort1_hits"],row_data["cohort2_hits"],row_data["cohort1_nonhits"],row_data["cohort2_nonhits"])),nrow=2,byrow = T);test=stats::fisher.test(mat);return(test$p.value)})
    OR_table$OR_test=pvalues
    OR_table<- OR_table[order(OR_table$OR_test,decreasing = F),]
    top_genes<-OR_table$Hugo_Symbol[1:max(sum(OR_table$OR_test<=0.25),20)]
    obj@output$OR_table<-OR_table
    mafcomparison<-maftools::mafCompare(m1=cohort1_maf,m2=cohort2_maf,m1Name = cohort1_term,m2Name = cohort2_term,minMut = 2)
    svg(filename = file.path(output_dir,"mafcompare_forestplot.svg"),width = 8,height=8)
    maftools::forestPlot(mafCompareRes = mafcomparison, pVal = 0.1)
    dev.off()
    svg(filename = file.path(output_dir,"mafcompare_cobarplot.svg"),width = 8,height=12)
    maftools::coBarplot(m1=cohort1_maf,m2=cohort2_maf,m1Name = cohort1_term,m2Name = cohort2_term,genes = top_genes)
    dev.off()
  }

  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  if(save_analysis){
    svg(filename = file.path(output_dir,"oncogenicpathwaytableplot.svg"),width = 10,height=8)
    print(oncogenicpathwaytableplot)
    dev.off()
    cat(obj@log,file=file.path(output_dir,"log.txt"))
  }
  return(obj)
}
)



#' runMDLAnalysis
#' run analysis from pyNBS output
#'
#' @param obj Analysis. analysis object to be run.
#' @param mdl_info data.frame. mdl sample information.
#' @param convert2maftools logic. whether convert wes mutation into maftools compatible.
#' @param maftools_mutation_conversion_dict character. a named vector mapping cgl mutation and maftools mutation.
#' @param mutation_panel_info data.frame. mdl mutation panel information.
#' @param with_oncokb logic. whether annotate with oncokb mutation.
#' @param oncokb_cancer_genes_info data.frame. oncokb database.
#' @param with_cosmic logic. whether annotate with cosmic mutation.
#' @param cosmic_mutations data.frame. cosmic mutation database.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param cancer_genes character. cancer gene list.
#' @param sample_order character. sample order for complex onco plot.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#' @param ... list. parameters passed to complex_onco_plot
#'
#' @return Analysis.
#' @export
#'
setGeneric("runMDLAnalysis",function(obj,mdl_info,convert2maftools,maftools_mutation_conversion_dict,mutation_panel_info,with_oncokb=T,oncokb_cancer_genes_info,with_cosmic=T,cosmic_mutations,top_anno,top_anno_sample_order,cancer_genes,sample_order=NULL,save_analysis=F,output_dir,log="run mutational pynbs analysis.",...) standardGeneric("runMDLAnalysis"))
setMethod("runMDLAnalysis","Analysis",function(obj,mdl_info,convert2maftools,maftools_mutation_conversion_dict,mutation_panel_info,with_oncokb=T,oncokb_cancer_genes_info,with_cosmic=T,cosmic_mutations,top_anno,top_anno_sample_order,cancer_genes,sample_order=NULL,save_analysis=F,output_dir,log="run mutational pynbs analysis.",...){

  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  }
  mdl_mutations<-obj@experiment@mdl_assay@assay_data
  #remove wildtype, annotate mdl_mutations with Variant_Classification and Variant_Type
  mdl_mutations<-mdl_mutations[mdl_mutations$indication_name!="wildtype",]
  mdl_mutations <- mdl_mutations %>% dplyr::group_by(Sample_ID,gene_name) %>% dplyr::mutate(n_annotator=dplyr::n_distinct(annotator_version_number)) %>% dplyr::filter(!(n_annotator>1 & annotator_version_number=="mdl_1.1"))  %>% dplyr::select(-n_annotator)#multiple annotated genes with old version?
  mdl_mutations <- mdl_mutations %>% dplyr::group_by(Patient_ID,Sample_ID,gene_name,document_subtype_description) %>% dplyr::mutate(conflicthits=dplyr::n())
  mdl_mutations<-mdl_mutations[!(mdl_mutations$conflicthits>1 & mdl_mutations$indication_name=="wildtype"),]# remove conflict wildtype (mutation and wildtype found in same sample and same gene),total 433 mutations(germline,mutation,variant)
  mdl_mutations$func.knowngene<-"exonic"
  mdl_mutations$del_nucleotide_start<-gsub("c\\.(\\d+).*", "\\1", mdl_mutations$nucleotide_description)
  mdl_mutations$del_nucleotide_end<-gsub("c\\..*_(\\d+).*", "\\1", mdl_mutations$nucleotide_description)
  mdl_mutations$del_nucleotide_end[!grepl("_",mdl_mutations$nucleotide_description)]<-mdl_mutations$del_nucleotide_start[!grepl("_",mdl_mutations$nucleotide_description)]
  mdl_mutations$del_nucleotide_length<-as.numeric(mdl_mutations$del_nucleotide_end)-as.numeric(mdl_mutations$del_nucleotide_start)+1
  mdl_mutations$ins_length<-NA
  mdl_mutations$ins_length[grepl("(delins)|(dup)",mdl_mutations$nucleotide_description)]<-nchar(gsub("(.*delins)|(.*dup)","",mdl_mutations$nucleotide_description[grepl("(delins)|(dup)",mdl_mutations$nucleotide_description)]))
  mdl_mutations$Variant_Classification<-NA
  mdl_mutations$Variant_Classification[grepl("del$",mdl_mutations$nucleotide_description,perl=T) & (grepl("Deletion",mdl_mutations$evidence)) & (mdl_mutations$del_nucleotide_length %% 3 !=0)]<-"Frame_Shift_Del"
  mdl_mutations$Variant_Classification[grepl("del$",mdl_mutations$nucleotide_description,perl=T) & (grepl("Deletion",mdl_mutations$evidence)) & (mdl_mutations$del_nucleotide_length %% 3 ==0)]<-"In_Frame_Del"
  mdl_mutations$Variant_Classification[grepl("dup",mdl_mutations$nucleotide_description,perl=T) & (mdl_mutations$ins_length %% 3 ==0)]<-"In_Frame_Ins"
  mdl_mutations$Variant_Classification[grepl("dup",mdl_mutations$nucleotide_description,perl=T) & (mdl_mutations$ins_length %% 3 !=0)]<-"Frame_Shift_Ins"
  mdl_mutations$Variant_Classification[grepl("delins",mdl_mutations$evidence) & (mdl_mutations$ins_length>mdl_mutations$del_nucleotide_length) & (abs(mdl_mutations$ins_length-mdl_mutations$del_nucleotide_length) %% 3==0)]<-"In_Frame_Ins"
  mdl_mutations$Variant_Classification[grepl("delins",mdl_mutations$evidence) & (mdl_mutations$ins_length>mdl_mutations$del_nucleotide_length) & (abs(mdl_mutations$ins_length-mdl_mutations$del_nucleotide_length) %% 3 !=0)]<-"Frame_Shift_Ins"
  mdl_mutations$Variant_Classification[grepl("delins",mdl_mutations$evidence) & (mdl_mutations$ins_length<mdl_mutations$del_nucleotide_length) & (abs(mdl_mutations$ins_length-mdl_mutations$del_nucleotide_length) %% 3 ==0)]<-"In_Frame_Del"
  mdl_mutations$Variant_Classification[grepl("delins",mdl_mutations$evidence) & (mdl_mutations$ins_length<mdl_mutations$del_nucleotide_length) & (abs(mdl_mutations$ins_length-mdl_mutations$del_nucleotide_length) %% 3 !=0)]<-"Frame_Shift_Del"
  mdl_mutations$Variant_Classification[grepl("Missense",mdl_mutations$evidence)]<-"Missense_Mutation"
  mdl_mutations$Variant_Classification[grepl("Nonsense",mdl_mutations$evidence)]<-"Nonsense_Mutation"
  mdl_mutations$Variant_Type<-"SNP"
  mdl_mutations$Variant_Type[grepl("Del",mdl_mutations$Variant_Classification)]<-"DEL"
  mdl_mutations$Variant_Type[grepl("Ins",mdl_mutations$Variant_Classification)]<-"INS"

  mdl_info<-mdl_info[mdl_info$Patient_ID %in% mdl_mutations$Patient_ID,]

  sample_info<-obj@experiment@sample_info
  wes_sample_info<-sample_info[sample_info$Assay=="DNA",]

  wes_mutations<-obj@experiment@snp_indel_assay@assay_data
  if(convert2maftools){
    wes_mutations<-utiltools::convert_mutations(wes_mutations,target = "maftools",maftools_mutation_conversion_dict = maftools_mutation_conversion_dict,)
  }

  mdl_wes_mutations<-utiltools::harmonize_mutations(wes_sample_info=wes_sample_info,wes_sample_id_col="Sample_ID",wes_patient_id_col="Patient_ID",wes_timepoint_col="Timepoint",wes_sample_collect_col="Collectdate",
                                         wes_mutations=wes_mutations,wes_mutation_sample_id_col="Tumor_Sample_Barcode",wes_gene_col="Hugo_Symbol",wes_chromosome_col="Chromosome",wes_start_col="Start_position",wes_end_col="End_position",wes_keep_columns=c("Variant_Classification","Variant_Type","Reference_Allele","Tumor_Seq_Allele1","Tumor_Seq_Allele2","Protein_Change","i_TumorVAF_WU","cdna"),
                                         mdl_sample_info=mdl_info,mdl_sample_id_col="Sample_ID",mdl_patient_id_col="Patient_ID",mdl_timepoint_col="Timepoint",mdl_sample_collect_time_col="Collected_At",
                                         mdl_mutations=mdl_mutations,mdl_mutation_sample_id_col="Sample_ID",mdl_mutation_patient_id_col="Patient_ID",mdl_mutation_collect_time_col="document_service_at",mdl_mutation_gene_col="gene_name",mdl_mutation_indication_col="indication_name",mdl_mutations_aa_col="protein_description",mdl_mutation_panel_col="document_subtype_description",mdl_mutation_keep_columns=c("indication_name","aa_change","reference_nucleotide","alternate_nucleotide","ref_seq_type","nucleotide_change" ,"nucleotide_change_type" ,"nucleotide_change_subtype"),
                                         panel_info=mutation_panel_info,panel_id_col="document_type",panel_chromosome_col="chromosome",panel_start_col="start_pos",panel_end_col="end_pos",panel_gene_col="gene",
                                         with_oncokb=T,oncokb_db=oncokb_cancer_genes_info,oncokb_gene_col="Hugo Symbol",with_cosmic=T,cosmic_db=cosmic_mutations,cosmic_gene_col="GENE_NAME",cosmic_aa_col="Mutation AA",
                                         combined=TRUE)


  mdl_sample_info<-sample_info[sample_info$Assay=="MDL",c("Patient_ID","Timepoint","Sample_ID")]
  DNA_sample_info<-sample_info[sample_info$Assay=="DNA",c("Patient_ID","Timepoint","Sample_ID")]
  mdl_sample2DNA_sample_<-merge(mdl_sample_info,DNA_sample_info,by=c("Patient_ID","Timepoint"),suffixes=c("_MDL","_DNA"))
  mdl_sample2DNA_sample<-mdl_sample2DNA_sample_$Sample_ID_DNA
  mdl_sample2DNA_sample<-setNames(mdl_sample2DNA_sample,nm=mdl_sample2DNA_sample_$Sample_ID_MDL)
  mdl_wes_mutations$Sample_ID<-sapply(mdl_wes_mutations$Sample_ID,function(item){items<-unlist(strsplit(item,";")); last_item=items[length(items)];ifelse(last_item %in% names(mdl_sample2DNA_sample),mdl_sample2DNA_sample[last_item],last_item)})
  DNA_consistent_mdl_wes_mutations<-mdl_wes_mutations[!mdl_wes_mutations$Sample_ID %in% mdl_sample_info$Sample_ID,]
  mutations<-DNA_consistent_mdl_wes_mutations
  mutations$Hugo_Symbol<-mutations$Gene
  mutations$Tumor_Sample_Barcode<-mutations[["Sample_ID"]]


  DNA_samples<-top_anno_sample_order
  #sample_info<-merge(merge(data.frame(Sample_ID=DNA_samples),obj@experiment@sample_info[,c("Sample_ID","Patient_ID")],by="Sample_ID",all.x=T),obj@experiment@patient_info,by="Patient_ID",all.x=T)
  #sample_info$Tumor_Sample_Barcode<-sample_info$Sample_ID
  mutations$Tumor_Sample_Barcode<-factor(mutations$Tumor_Sample_Barcode,levels=DNA_samples)

  Mutsig2CV_directory<-obj@experiment@snp_indel_assay@Mutsig2CV_directory
  mutsig2cv_sig_genes_file<-file.path(Mutsig2CV_directory,"sig_genes.txt")
  if(!file.exists(mutsig2cv_sig_genes_file)){cat("Please run mutsigcv2 first")} else{
    sig_genes<-read.csv(mutsig2cv_sig_genes_file,header=T,stringsAsFactors = F,check.names = F,sep="\t")
  }
  selected_genes<-sig_genes$gene[sig_genes$p<=0.05 & sig_genes$nnon>=ifelse(length(unique(mutations$Tumor_Sample_Barcode))>20,2,1)]
  if(!missing(cancer_genes)){
    selected_genes<-selected_genes[selected_genes %in% cancer_genes]
  }

  if(!is.null(sample_order)){
    top_anno<-top_anno[match(sample_order,DNA_samples),]
    DNA_samples<-sample_order
    mutations$Tumor_Sample_Barcode<-factor(mutations$Tumor_Sample_Barcode,levels=sample_order)
    sample_info<-sample_info[match(DNA_samples,sample_info$Tumor_Sample_Barcode),]
  }
  sequence_methods<-rep("WES",length(DNA_samples))
  sequence_methods<-setNames(sequence_methods,DNA_samples)
  sequence_methods[unname(mdl_sample2DNA_sample)]<-"MDL;WES"
  sequence_methods_anno<-ComplexHeatmap::HeatmapAnnotation(Method=factor(unname(sequence_methods)),annotation_name_side = "left",col=list(Method=c("WES"="grey","MDL;WES"="blue")))

  top_anno<-c(top_anno,sequence_methods_anno)
  snp_indel_cnv_col<-c("Missense_Mutation"="#33A02B","Nonsense_Mutation"="#ff0000","In_Frame_Del"="#FEFF99","In_Frame_Ins"="#938038","Frame_Shift_Del"="#1F78B4","Frame_Shift_Ins"="#7b42f5","Splice_Site"="#00ffd5","Nonstop_Mutation"="#0022fc","UTR3"="#707173","UTR5"="#464647","Multi_Hit"="#a68a6c")
  complex_onco_plot<-utiltools::complex_oncoplot(snp_indels=mutations,selected_genes = selected_genes,mutsigCV2_sig_genes = sig_genes,multi_hit=TRUE,cnv_types = c("Loss","Gain","AMP","Del") ,show_heatmap_legend = T,cnv=F,top_filter = 100L,remove_macromolecular_gene = T,col=snp_indel_cnv_col,multi_hit_col = "#a68a6c",macromolecular_threshold = 300000,top_annotation = top_anno,sample_order =sample_order,annotate_source = T,source_column = "Support",output_dir=output_dir,...)

  obj@output$DNA_consistent_mdl_wes_mutations<-DNA_consistent_mdl_wes_mutations
  obj@output$complex_onco_plot<-complex_onco_plot
  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  if(save_analysis){
    write.csv(DNA_consistent_mdl_wes_mutations,file.path(output_dir,"DNA_consistent_mdl_wes_mutations.csv"))
    cat(obj@log,file=file.path(output_dir,"log.txt"))
  }
  return(obj)
}
)



#' runpyNBSAnalysis
#' run analysis from pyNBS output
#'
#' @param obj Analysis. analysis object to be run.
#' @param pyNBS_dir character. path to directory of pyNBS output.
#' @param survival_column character. column for survival days.
#' @param status_column character. column for survival status.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param palette character. colors for KM curve.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runpyNBSAnalysis",function(obj,pyNBS_dir,survival_column,status_column,top_anno,top_anno_sample_order,palette,save_analysis=F,output_dir,log="run mutational pynbs analysis.") standardGeneric("runpyNBSAnalysis"))
setMethod("runpyNBSAnalysis","Analysis",function(obj,pyNBS_dir,survival_column,status_column,top_anno,top_anno_sample_order,palette,save_analysis=F,output_dir,log="run mutational pynbs analysis."){

  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  }
browser()
  DNA_samples<-top_anno_sample_order
  sample_info<-merge(merge(data.frame(Sample_ID=DNA_samples),obj@experiment@sample_info,by="Sample_ID",all.x=T),obj@experiment@patient_info,by="Patient_ID",all.x=T)
  sample_info$Tumor_Sample_Barcode<-sample_info$Sample_ID
  colnames(sample_info)[match(c(survival_column,status_column),colnames(sample_info))]<-c("Survival","Status")

  #pyNBS cluster
  stopifnot("No pyNBS directory found. Please run pyNBS first."=dir.exists(pyNBS_dir))
  NBS_cc_table_file<-list.files(pyNBS_dir,pattern = "cc_matrix")
  NBS_cluster_assignment_file<-list.files(pyNBS_dir,pattern = "cluster_assignments")
  pynbs_cc<-read.csv(file.path(pyNBS_dir,NBS_cc_table_file),header=T,row.names = 1,check.names = F)
  pynbs_cc<-pynbs_cc[DNA_samples,DNA_samples]
  obj@output$pynbs_cc<-pynbs_cc

  pynbs_cluster<-read.csv(file.path(pyNBS_dir,NBS_cluster_assignment_file),header=F,row.names = 1)
  colnames(pynbs_cluster)<-"Cluster"
  pynbs_cluster[["pyNBS_Cluster"]]<-paste("Cluster_",pynbs_cluster[["Cluster"]],sep="")
  obj@output$pynbs_cluster<-pynbs_cluster
  clusters<-length(unique(pynbs_cluster[["pyNBS_Cluster"]]))
  pyNBS_Cluster_col<-RColorBrewer::brewer.pal(n=clusters,name="Paired")[1:clusters]
  names(pyNBS_Cluster_col)<-unique(pynbs_cluster[["pyNBS_Cluster"]])
  pynbs_cluster_anno<-ComplexHeatmap::HeatmapAnnotation(df=pynbs_cluster[DNA_samples,"pyNBS_Cluster",drop=F],col=list(pyNBS_Cluster=pyNBS_Cluster_col))
  pynbs_cc_heatmap<-ComplexHeatmap::Heatmap(pynbs_cc,top_annotation = c(top_anno,pynbs_cluster_anno),name="coef")
  sample_info$pyNBS_Cluster[match(rownames(pynbs_cluster),sample_info$Tumor_Sample_Barcode)]<-pynbs_cluster[['pyNBS_Cluster']]
  pyNBSOS<-sample_info %>% survivalAnalysis::analyse_survival(dplyr::vars(Survival,Status),by=pyNBS_Cluster)

  if(missing(palette)){
    palette<-RColorBrewer::brewer.pal(n=clusters,name="Paired")[1:clusters]
  }
  pynbs_km_plot<-survivalAnalysis::kaplan_meier_plot(pyNBSOS,
                                   break.time.by="breakByMonthYear",
                                   xlab="Survival (month)",
                                   legend.title="pyNBS_Cluster",
                                   hazard.ratio=TRUE,
                                   risk.table=TRUE,
                                   table.layout="clean",
                                   ggtheme=ggplot2::theme_bw(10),
                                   palette=palette,
                                   legend=c(0.8,0.8))

  obj@output$pynbs_cc_heatmap<-pynbs_cc_heatmap
  obj@output$pynbs_km_plot<-pynbs_km_plot
  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  if(save_analysis){
    svg(filename = file.path(output_dir,"pynbs_cc_heatmap.svg"),width = 8,height=8)
    print(pynbs_cc_heatmap)
    dev.off()
    svg(filename = file.path(output_dir,"pynbs_km_plot.svg"),width = 8,height=8)
    print(pynbs_km_plot)
    dev.off()
    cat(obj@log,file=file.path(output_dir,"log.txt"))
  }
  return(obj)
}
)


#' runCloneunSupervisedAnalysis
#' run unsupervised analysis of Canopy derived clone deconvolution
#'
#' @param obj Analysis. analysis object to be run.
#' @param survival_column character. column for survival days.
#' @param status_column character. column for survival status.
#' @param stablecluster_method character. parameter passed to stableCluster
#' @param cutFUN function. parameter passed to stableCluster
#' @param clusters numeric.parameter passed to stableCluster
#' @param nTimes numeric. parameter passed to stableCluster
#' @param subSampleSize numeric. parameter passed to stableCluster
#' @param subFeatureSize numeric. parameter passed to stableCluster
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param palette character. colors for KM curve.
#' @param pval_cutoff numeric. p value cutoff for significance.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runCloneunSupervisedAnalysis",function(obj,
                                                   survival_column,status_column,
                                                   stablecluster_method="combine",cutFUN=c(ClassDiscovery::cutHclust),clusters,nTimes=100,subSampleSize = 0.9,subFeatureSize=0.9,
                                                   top_anno,top_anno_sample_order,palette,
                                                   pval_cutoff = 0.05,
                                                   save_analysis = T,output_dir,
                                                   log="run canopy clone fraction unsuperved analysis.") standardGeneric("runCloneunSupervisedAnalysis"))
setMethod("runCloneunSupervisedAnalysis","Analysis",function(obj,
                                                             survival_column,status_column,
                                                             stablecluster_method="combine",cutFUN=c(ClassDiscovery::cutHclust),clusters,nTimes=100,subSampleSize = 0.9,subFeatureSize=0.9,
                                                             top_anno,top_anno_sample_order,palette,
                                                             pval_cutoff = 0.05,
                                                             save_analysis = T,output_dir,
                                                             log="run canopy clone fraction unsupervised analysis."){

  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  } else {
    output_dir<-getwd()
  }
  #browser()
  DNA_samples<-top_anno_sample_order
  sample_info<-merge(merge(data.frame(Sample_ID=DNA_samples),obj@experiment@sample_info,by="Sample_ID",all.x=T),obj@experiment@patient_info,by="Patient_ID",all.x=T)
  sample_info<-sample_info[match(DNA_samples, sample_info$Sample_ID),]

  ccf<-as.data.frame(obj@experiment@clone_assay@assay_data[,DNA_samples])


  if(missing(clusters)){
    clusters<-utiltools::estimate_bestNumberofClusters(ccf,max.nc = as.integer(ncol(ccf)/2))$Best.NumberofCluster
  }

  stableclusters<-utiltools::stableCluster(ccf,method = stablecluster_method,cutFUN = cutFUN,clusters = clusters,nTimes = nTimes,subSampleSize = subSampleSize,subFeatureSize = subFeatureSize)
  assertthat::are_equal(colnames(stableclusters),colnames(ccf))

  stablecluster_ccf_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(ccf))),name="ccf",cluster_columns = hclust(dist(stableclusters)),top_annotation = top_anno,column_split = clusters)
  stablecluster_heatmap_data<-ccf[ComplexHeatmap::row_order(stablecluster_ccf_heatmap),unlist(ComplexHeatmap::column_order(stablecluster_ccf_heatmap))]

  #stablecluster survival anaysis
  if(missing(palette)){
    palette=RColorBrewer::brewer.pal(n=clusters,name="Paired")[1:clusters]
  }
  cluster_km_plot<-utiltools::hCluster_Surv(heatmap=stablecluster_ccf_heatmap,clinics=sample_info,survival_col = survival_column,status_col = status_column,palette = palette,legend=c(0.8,0.2))

  #differential comparison between clusters
  sample_info$CCF_Cluster[unlist(ComplexHeatmap::column_order(stablecluster_ccf_heatmap))]<-rep(paste("Cluster_",1:clusters,sep=""),times=sapply(ComplexHeatmap::column_order(stablecluster_ccf_heatmap),length))
  if(length(unique(sample_info$CCF_Cluster))>2){
    cluster_statistics=apply(ccf,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[["CCF_Cluster"]]);test=aov(value~group,data=df);return(summary(test)[[1]][["Pr(>F)"]][1])})
  } else {
    cluster_statistics=apply(ccf,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[["CCF_Cluster"]]);test=t.test(value~group,data=df);return(test$p.value)})
  }

  neg_log_pvalue=-log10(cluster_statistics)
  row_anno<-ComplexHeatmap::rowAnnotation(`-log10(p)`=ComplexHeatmap::row_anno_barplot(neg_log_pvalue,gp=grid::gpar(fill=ifelse(neg_log_pvalue>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)),width=grid::unit(3,"cm"),annotation_name_side="bottom",annotation_name_rot=0)

  fontsize=rep(8,nrow(ccf))
  fontsize[neg_log_pvalue>=(-log10(pval_cutoff))]<-10
  fontface=rep("plain",nrow(ccf))
  fontface[neg_log_pvalue>=(-log10(pval_cutoff))]<-"bold"

  stablecluster_ccf_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(ccf))),name="ccf",cluster_columns = hclust(dist(stableclusters)),top_annotation = top_anno,column_split = clusters,right_annotation = row_anno,row_names_gp = grid::gpar(fontsize=fontsize,fontface=fontface))

  cluster_mean_<-list()
  for (ccf_cluster in unique(sample_info$CCF_Cluster)){
    cluster_mean_[[paste(ccf_cluster,"_mean",sep="")]]<-rowMeans(ccf[,sample_info[["CCF_Cluster"]]==ccf_cluster],na.rm=T)
  }
  cluster_mean_<-as.data.frame(cluster_mean_)
  cluster_results<-data.frame(cluster_mean_[names(cluster_statistics),],P.Value=cluster_statistics)
  cluster_results[["adj.P.Val"]]<-p.adjust(cluster_results[["P.Value"]],method='BH')

  if(save_analysis){
    write.csv(sample_info,file.path(output_dir,"sample_info.csv"))
    write.csv(ccf,file.path(output_dir,"ccf.csv"))
    if(!dir.exists(file.path(output_dir,"unsupervised_cluster"))){dir.create(file.path(output_dir,"unsupervised_cluster"),recursive = T)}
    write.csv(cluster_results,file.path(output_dir,"unsupervised_cluster","stablecluster_statistics.csv"))
    svg(filename = file.path(output_dir,"unsupervised_cluster","stablecluster_ccf_heatmap"),width = 8,height=10)
    ComplexHeatmap::draw(stablecluster_ccf_heatmap)
    if(any(neg_log_pvalue>(-(log10(pval_cutoff))))){
      xpos<-0.95/(max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
      ComplexHeatmap::decorate_annotation("-log10(p)",{
        grid::grid.lines(c(xpos,xpos),c(0,1),gp = grid::gpar(lty = 2,col="red"))
      })
    }
    dev.off()
    svg(filename = file.path(output_dir,"unsupervised_cluster","cluster_km_plot"),width = 8,height=8)
    print(cluster_km_plot)
    dev.off()
    cat(obj@log,file=file.path(output_dir,"log.txt"))
  }


  obj@output$sample_info<-sample_info
  obj@output$ccf<-ccf
  obj@output$stableclusters=stableclusters
  obj@output$unsupervised_cluster$stablecluster_ccf_heatmap<-stablecluster_ccf_heatmap
  obj@output$unsupervised_cluster$stablecluster_statistics<-cluster_results
  obj@output$unsupervised_cluster$cluster_km_plot<-cluster_km_plot

  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")

  return(obj)
}
)



#' runCCFDifAnalysis
#'
#' run differentiate analysis of canopy derived cell fraction
#'
#' @param obj Analysis. analysis object to be run.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param palette character. colors for KM curve.
#' @param group character. comparison group or column. parameter passed to aov or ttest
#' @param block character. comparison block used for paired ttest.
#' @param contrasts character. contrasts comparison post aov
#' @param pval_cutoff numeric. p value cutoff for significance.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runCCFDifAnalysis",function(obj,
                                        top_anno,top_anno_sample_order,palette,
                                        group,block,contrasts,
                                        pval_cutoff = 0.05,
                                        save_analysis = T,output_dir,
                                        log="run canopy clone fraction differential analysis.") standardGeneric("runCCFDifAnalysis"))
setMethod("runCCFDifAnalysis","Analysis",function(obj,
                                                  top_anno,top_anno_sample_order,palette,
                                                  group,block,contrasts,
                                                  pval_cutoff = 0.05,
                                                  save_analysis = T,output_dir,
                                                  log="run canopy clone fraction differential analysis."){

  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  } else {
    output_dir<-getwd()
  }
  #browser()
  DNA_samples<-top_anno_sample_order
  sample_info<-merge(merge(data.frame(Sample_ID=DNA_samples),obj@experiment@sample_info,by="Sample_ID",all.x=T),obj@experiment@patient_info,by="Patient_ID",all.x=T)
  sample_info<-sample_info[match(DNA_samples, sample_info$Sample_ID),]

  ccf<-as.data.frame(obj@experiment@clone_assay@assay_data[,DNA_samples])

  group_mean_<-list()
  for (ccf_group in unique(sample_info[[group]])){
    group_mean_[[paste(ccf_group,"_mean",sep="")]]<-rowMeans(ccf[,sample_info[[group]]==ccf_group],na.rm=T)
  }
  group_results<-as.data.frame(group_mean_)
  if(missing(block)){
    if(missing(contrasts)){
      if(length(unique(sample_info[[group]]))>2){
        group_statistics=apply(ccf,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[[group]]);test=aov(value~group,data=df);return(summary(test)[[1]][["Pr(>F)"]][1])})
      } else {
        group_statistics=apply(ccf,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[[group]]);test=t.test(value~group,data=df);return(test$p.value)})
      }
      group_results<-data.frame(group_results,P.Value=group_statistics)
      neg_log_pvalue=-log10(group_statistics)
      row_anno<-ComplexHeatmap::rowAnnotation(`-log10(p)`=ComplexHeatmap::row_anno_barplot(neg_log_pvalue,gp=grid::gpar(fill=ifelse(neg_log_pvalue>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)),width=grid::unit(3,"cm"),annotation_name_side="bottom",annotation_name_rot=0)
      fontsize=rep(8,nrow(ccf))
      fontsize[neg_log_pvalue>=(-log10(pval_cutoff))]<-10
      fontface=rep("plain",nrow(ccf))
      fontface[neg_log_pvalue>=(-log10(pval_cutoff))]<-"bold"
      ccf_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(ccf))),name="zscores_ccf",cluster_columns =ComplexHeatmap::cluster_between_groups(t(scale(t(ccf))),factor = sample_info[[group]]),top_annotation = top_anno,column_split = length(unique(sample_info[[group]])),right_annotation = row_anno,row_names_gp = grid::gpar(fontsize=fontsize,fontface=fontface))
      obj@output$contrasts$ccf_heatmap<-ccf_heatmap
      obj@output$contrasts$statistics<-group_results
      if(save_analysis){
        if(!dir.exists(file.path(output_dir,"contrasts"))){dir.create(file.path(output_dir,"contrasts"),recursive = T)}
        write.csv(group_results,file.path(output_dir,"contrasts","statistics.csv"))
        svg(filename = file.path(output_dir,"contrasts","ccf_heatmap.svg"),width = 8,height=8)
        ComplexHeatmap::draw(ccf_heatmap)
        if(any(neg_log_pvalue>(-(log10(pval_cutoff))))){
          xpos<-0.95/(max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
          ComplexHeatmap::decorate_annotation("-log10(p)",{
            grid::grid.lines(c(xpos,xpos),c(0,1),gp = grid::gpar(lty = 2,col="red"))
          })
        }
        dev.off()
      }
    } else {
      for(contrast in contrasts){
        subset_sample_info<-sample_info[sample_info[[group]] %in% unlist(strsplit(contrast,"-")),]
        subset_ccf<-ccf[,subset_sample_info$Sample_ID]
        contrast_group_statistics=apply(subset_ccf,1,function(scf){df<-data.frame(value=as.numeric(scf),group=subset_sample_info[[group]]);test=t.test(value~group,data=df);return(test$p.value)})
        contrast_group_results[[paste(contrast,"P.Value",sep=":")]]<-contrast_group_statistics
        contrast_top_anno<-top_anno[match(subset_sample_info$Sample_ID,DNA_samples),]
        neg_log_pvalue=-log10(group_statistics)
        row_anno<-ComplexHeatmap::rowAnnotation(`-log10(p)`=ComplexHeatmap::row_anno_barplot(neg_log_pvalue,gp=grid::gpar(fill=ifelse(neg_log_pvalue>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)),width=grid::unit(3,"cm"),annotation_name_side="bottom",annotation_name_rot=0)
        fontsize=rep(8,nrow(ccf))
        fontsize[neg_log_pvalue>=(-log10(pval_cutoff))]<-10
        fontface=rep("plain",nrow(ccf))
        fontface[neg_log_pvalue>=(-log10(pval_cutoff))]<-"bold"
        contrast_ccf_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(subset_ccf))),name="zscores_ccf",cluster_columns =ComplexHeatmap::cluster_between_groups(t(scale(t(subset_ccf))),factor = subset_sample_info[[group]]),top_annotation = contrast_top_anno,column_split = length(unique(subset_sample_info[[group]])),right_annotation = row_anno,row_names_gp = grid::gpar(fontsize=fontsize,fontface=fontface))

        obj@output$contrasts[[contrast]]$sample_info<-subset_sample_info
        obj@output$contrasts[[contrast]]$ccf<-subset_ccf
        obj@output$contrasts[[contrast]]$ccf_heatmap<-contrast_ccf_heatmap
        obj@output$contrasts[[contrast]]$statistics<-contrast_group_results
        if(save_analysis){
          if(!dir.exists(file.path(output_dir,"contrasts",contrast))){dir.create(file.path(output_dir,"contrasts",contrast),recursive = T)}
          write.csv(contrast_ccf_heatmap,file.path(output_dir,"contrasts",contrast,"statistics.csv"))
          write.csv(subset_sample_info,file.path(output_dir,"contrasts",contrast,"sample_info.csv"))
          write.csv(subset_ccf,file.path(output_dir,"contrasts",contrast,"subset_ccf.csv"))
          svg(filename = file.path(output_dir,"contrasts",contrast,"contrast_ccf_heatmap.svg"),width = 8,height=8)
          ComplexHeatmap::draw(contrast_ccf_heatmap)
          if(any(neg_log_pvalue>(-(log10(pval_cutoff))))){
            xpos<-0.95/(max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
            ComplexHeatmap::decorate_annotation("-log10(p)",{
              grid::grid.lines(c(xpos,xpos),c(0,1),gp = grid::gpar(lty = 2,col="red"))
            })
          }
          dev.off()
        }
      }
    }
  } else {
    sample_info<-sample_info[order(sample_info[[block]],sample_info[[group]]),]
    ccf<-ccf[,sample_info$Sample_ID]

    if(missing(contrasts)){
      if(length(unique(sample_info[[group]]))>2){
        group_statistics=apply(ccf,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[[group]]);test=aov(value~group,data=df);return(summary(test)[[1]][["Pr(>F)"]][1])})
      } else {
        group_statistics=apply(ccf,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[[group]]);test=t.test(value~group,data=df,paired=T);return(test$p.value)})
      }
      group_results<-data.frame(group_results,P.Value=group_statistics)
      neg_log_pvalue=-log10(group_statistics)
      row_anno<-ComplexHeatmap::rowAnnotation(`-log10(p)`=ComplexHeatmap::row_anno_barplot(neg_log_pvalue,gp=grid::gpar(fill=ifelse(neg_log_pvalue>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)),width=grid::unit(3,"cm"),annotation_name_side="bottom",annotation_name_rot=0)
      fontsize=rep(8,nrow(ccf))
      fontsize[neg_log_pvalue>=(-log10(pval_cutoff))]<-10
      fontface=rep("plain",nrow(ccf))
      fontface[neg_log_pvalue>=(-log10(pval_cutoff))]<-"bold"
      ccf_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(ccf))),name="zscores_ccf",cluster_columns=ComplexHeatmap::cluster_between_groups(t(scale(t(ccf))),factor = sample_info[[group]]),top_annotation = top_anno,column_split = length(unique(sample_info[[group]])),right_annotation = row_anno,row_names_gp = grid::gpar(fontsize=fontsize,fontface=fontface))

      obj@output$contrasts$ccf_heatmap<-ccf_heatmap
      obj@output$contrasts$statistics<-group_results
      if(save_analysis){
        if(!dir.exists(file.path(output_dir,"contrasts"))){dir.create(file.path(output_dir,"contrasts"),recursive = T)}
        write.csv(group_results,file.path(output_dir,"contrasts","statistics.csv"))
        svg(filename = file.path(output_dir,"contrasts","ccf_heatmap.svg"),width = 8,height=8)
        ComplexHeatmap::draw(ccf_heatmap)
        if(any(neg_log_pvalue>(-(log10(pval_cutoff))))){
          xpos<-0.95/(max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
          ComplexHeatmap::decorate_annotation("-log10(p)",{
            grid::grid.lines(c(xpos,xpos),c(0,1),gp = grid::gpar(lty = 2,col="red"))
          })
        }
        dev.off()
      }
    } else {
      for(contrast in contrasts){
        subset_sample_info<-sample_info[sample_info[[group]] %in% unlist(strsplit(contrast,"-")),]
        subset_ccf<-ccf[,subset_sample_info$Sample_ID]
        contrast_group_statistics=apply(subset_ccf,1,function(scf){df<-data.frame(value=as.numeric(scf),group=subset_sample_info[[group]]);test=t.test(value~group,data=df,paired=T);return(test$p.value)})
        contrast_group_results[[paste(contrast,"P.Value",sep=":")]]<-contrast_group_statistics
        contrast_top_anno<-top_anno[match(subset_sample_info$Sample_ID,DNA_samples),]
        neg_log_pvalue=-log10(group_statistics)
        row_anno<-ComplexHeatmap::rowAnnotation(`-log10(p)`=ComplexHeatmap::row_anno_barplot(neg_log_pvalue,gp=grid::gpar(fill=ifelse(neg_log_pvalue>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)),width=grid::unit(3,"cm"),annotation_name_side="bottom",annotation_name_rot=0)
        fontsize=rep(8,nrow(ccf))
        fontsize[neg_log_pvalue>=(-log10(pval_cutoff))]<-10
        fontface=rep("plain",nrow(ccf))
        fontface[neg_log_pvalue>=(-log10(pval_cutoff))]<-"bold"
        contrast_ccf_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(subset_ccf))),name="zscores_ccf",cluster_columns=ComplexHeatmap::cluster_between_groups(t(scale(t(subset_ccf))),factor = subset_sample_info[[group]]),top_annotation = contrast_top_anno,column_split = length(unique(subset_sample_info[[group]])),right_annotation = row_anno,row_names_gp = grid::gpar(fontsize=fontsize,fontface=fontface))

        obj@output$contrasts[[contrast]]$sample_info<-subset_sample_info
        obj@output$contrasts[[contrast]]$ccf<-subset_ccf
        obj@output$contrasts[[contrast]]$ccf_heatmap<-contrast_ccf_heatmap
        obj@output$contrasts[[contrast]]$statistics<-contrast_group_results
        if(save_analysis){
          if(!dir.exists(file.path(output_dir,"contrasts",contrast))){dir.create(file.path(output_dir,"contrasts",contrast),recursive = T)}
          write.csv(contrast_ccf_heatmap,file.path(output_dir,"contrasts",contrast,"statistics.csv"))
          write.csv(subset_sample_info,file.path(output_dir,"contrasts",contrast,"sample_info.csv"))
          write.csv(subset_ccf,file.path(output_dir,"contrasts",contrast,"subset_ccf.csv"))
          svg(filename = file.path(output_dir,"contrasts",contrast,"contrast_ccf_heatmap.svg"),width = 8,height=8)
          ComplexHeatmap::draw(contrast_ccf_heatmap)
          if(any(neg_log_pvalue>(-(log10(pval_cutoff))))){
            xpos<-0.95/(max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
            ComplexHeatmap::decorate_annotation("-log10(p)",{
              grid::grid.lines(c(xpos,xpos),c(0,1),gp = grid::gpar(lty = 2,col="red"))
            })
          }
          dev.off()
        }
      }
    }
  }
  obj@output$sample_info<-sample_info
  obj@output$ccf<-ccf

  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")

  return(obj)
}
)



#' deltaAnalysis
#'
#' delta assay anlysis
#'
#' @param assay data.frame. assay data.
#' @param scale logic. whether to scale assay.
#' @param assay_name character. assay name.
#' @param sample_info data.frame. sample information.
#' @param patient_id_column character. column for patient id.
#' @param sample_id_column character. column for sample id.
#' @param timepoint_column character. column for timepoint.
#' @param timepoint_terms character. timepoint terms pair as delta.
#' @param group_column character. column for group.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_patient_order character. patient orders for top annotation.
#' @param test_method character. test method, such as "ttest" or "limma".
#' @param contrasts character. specific contrast.
#' @param pval_cutoff numeric. p value cutoff for significance.
#' @param only_sig_genes logic. whether only significant genes will be used.
#' @param only_show_sig_genes logic. whether only significant genes will be showed in heatmap.
#' @param palette character. colors for KM curve.
#' @param output_dir character. output directory.
#'
#' @return list. contain combined_assay,combined_assay_heapmap,combined_statistics,group_terms,groups.
#' @export
#'

deltaAnalysis<-function(assay,scale=F,assay_name="assay",
                        sample_info,patient_id_column="Patient_ID",sample_id_column="Sample_ID",timepoint_column,timepoint_terms=c("baseline"="Baseline","tp2"="TP2"),group_column,
                        top_anno,top_anno_patient_order,
                        test_method=c("ttest","limma"),contrasts,pval_cutoff = 0.05,
                        only_sig_genes=F,only_show_sig_genes=F,palette,
                        output_dir){
  if(scale){
    assay<-as.data.frame(t(scale(t(assay))))
    heatmap_name=paste("zscore",assay_name,sep="_")
  } else {
    heatmap_name=assay_name
  }
  test_method=match.arg(test_method)
  baseline_samples_info<-sample_info[sample_info[[timepoint_column]]==timepoint_terms["baseline"],]
  baseline_sample_info<-baseline_samples_info[match(top_anno_patient_order,baseline_samples_info[[patient_id_column]]),]
  baseline_samples<-baseline_sample_info[[sample_id_column]]
  tp2_samples_info<-sample_info[sample_info[[timepoint_column]]==timepoint_terms["tp2"],]
  tp2_samples_info<-tp2_samples_info[match(top_anno_patient_order,tp2_samples_info[[patient_id_column]]),]
  tp2_samples<-tp2_samples_info[[sample_id_column]]
  assertthat::are_equal(sample_info[[patient_id_column]][match(baseline_samples,sample_info[[sample_id_column]])],sample_info[[patient_id_column]][match(tp2_samples,sample_info[[sample_id_column]])])
  groups<-baseline_sample_info[[group_column]]

  baseline_assay<-assay[,baseline_samples]
  colnames(baseline_assay)<-baseline_sample_info[[patient_id_column]]
  tp2_assay<-assay[,tp2_samples]
  colnames(tp2_assay)<-tp2_samples_info[[patient_id_column]]
  assertthat::are_equal(colnames(baseline_assay),colnames(tp2_assay))
  delta_assay<-tp2_assay-baseline_assay
  assertthat::are_equal(colnames(baseline_assay),colnames(delta_assay))


  delta_statistics<-apply(delta_assay,1,function(r_data){data=data.frame(Value=as.numeric(r_data),Group=groups);test=stats::t.test(Value~Group,data=data);return(test$p.value)})

  if(only_sig_genes) {
    sig_genes<-names(delta_statistics)[delta_statistics<=pval_cutoff]
    delta_assay_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(delta_assay[sig_genes,]))),name=paste("zscore",assay_name,sep="_"),column_split = groups)
    feature_order<-sig_genes[ComplexHeatmap::row_order(delta_assay_heatmap)]
  } else {
    delta_assay_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(delta_assay))),name=paste("zscore",assay_name,sep="_"),column_split = groups)
    feature_order<-rownames(delta_assay)[ComplexHeatmap::row_order(delta_assay_heatmap)]
  }


  baseline_assay<-baseline_assay[feature_order,]
  tp2_assay<-tp2_assay[feature_order,]
  delta_assay<-delta_assay[feature_order,]
  combined_assay<-rbind(as.matrix(baseline_assay),as.matrix(tp2_assay),as.matrix(delta_assay))
  assertthat::are_equal(colnames(combined_assay),colnames(delta_assay))
  assertthat::are_equal(rownames(combined_assay),rep(feature_order,3))

  group_terms<-rev(sort(unique(as.character(groups))))
  means_<-data.frame(term1=rowMeans(combined_assay[,groups==group_terms[1]],na.rm=T),term2=rowMeans(combined_assay[,groups==group_terms[2]],na.rm=T))
  colnames(means_)<-group_terms
  if(test_method=="ttest"){
    baseline_statistics<-apply(baseline_assay,1,function(r_data){data=data.frame(Value=as.numeric(r_data),Group=groups);test=stats::t.test(Value~Group,data=data);return(test$p.value)})
    baseline_statistics<-baseline_statistics[feature_order]
    tp2_statistics<-apply(tp2_assay,1,function(r_data){data=data.frame(Value=as.numeric(r_data),Group=groups);test=stats::t.test(Value~Group,data=data);return(test$p.value)})
    tp2_statistics<-tp2_statistics[feature_order]
    delta_statistics<-delta_statistics[feature_order]
  }
  if(test_method=="limma"){
    baseline_limma<-utiltools::dge_limma(baseline_assay,is_rawcount = F,is_logged = T,normalize = F,sample_frequency_threshold = 0.5,clinic_info =data.frame(Sample_ID=colnames(baseline_assay),Group=groups),ID_col = "Sample_ID",group_col = "Group",contrasts = contrasts,method ="limma_trend")
    baseline_statistics_<-baseline_statistics$statistics[[colnames(baseline_statistics$statistics)[grepl("p.value",colnames(baseline_statistics$statistics),ignore.case=TRUE)]]]
    baseline_statistics<-stats::setNames(baseline_statistics_,rownames(baseline_statistics$statistics))
    baseline_statistics<-baseline_statistics[feature_order]

    tp2_limma<-utiltools::dge_limma(tp2_assay,is_rawcount = F,is_logged = T,normalize = F,sample_frequency_threshold = 0.5,clinic_info =data.frame(Sample_ID=colnames(tp2_assay),Group=groups),ID_col = "Sample_ID",group_col = "Group",contrasts = contrasts,method ="limma_trend")
    tp2_statistics_<-tp2_statistics$statistics[[colnames(tp2_statistics$statistics)[grepl("p.value",colnames(tp2_statistics$statistics),ignore.case=TRUE)]]]
    tp2_statistics<-stats::setNames(tp2_statistics_,rownames(tp2_statistics$statistics))
    tp2_statistics<-tp2_statistics[feature_order]

    delta_statistics<-apply(delta_assay,1,function(r_data){data=data.frame(Value=as.numeric(r_data),Group=groups);test=stats::t.test(Value~Group,data=data);return(test$p.value)})
    delta_statistics<-delta_statistics[feature_order]
  }

  combined_statistics<-c(baseline_statistics,tp2_statistics,delta_statistics)

  feature_splits<-factor(c(rep("Baseline",nrow(baseline_assay)),rep("TP2",nrow(tp2_assay)),rep("Delta",nrow(delta_assay))),levels=c("Baseline","TP2","Delta"),labels=c("Baseline","TP2","Delta"))


  combined_neg_log_statistics=-log10(combined_statistics)
  if(missing(palette)){palette=c("blue","white","red")}
  col=circlize::colorRamp2(breaks=c(min(combined_assay),0,max(combined_assay)),colors = palette)

  if(only_show_sig_genes){
    row_anno<-ComplexHeatmap::rowAnnotation(`-log10(p)`=ComplexHeatmap::anno_barplot(combined_neg_log_statistics,gp=grid::gpar(fill=ifelse(combined_neg_log_statistics>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(combined_neg_log_statistics,-log10(pval_cutoff)))+0.5)),
                                            sig_genes=ComplexHeatmap::anno_mark(at=unname(which(combined_neg_log_statistics>=(-log10(pval_cutoff)))),labels=names(which(combined_neg_log_statistics>=(-log10(pval_cutoff)))),labels_gp=grid::gpar(fontsize=5)),
                                            width=grid::unit(5,"cm"),show_annotation_name = T,annotation_name_side="bottom",annotation_name_rot=0)
    combined_assay_heatmap<-ComplexHeatmap::Heatmap(combined_assay,name=heatmap_name,col = col,cluster_rows = F,top_annotation = top_anno,column_split = groups,right_annotation = row_anno,show_row_names = F,row_split = feature_splits)

  } else {
    row_anno<-ComplexHeatmap::rowAnnotation(`-log10(p)`=ComplexHeatmap::anno_barplot(combined_neg_log_statistics,gp=grid::gpar(fill=ifelse(combined_neg_log_statistics>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(combined_neg_log_statistics,-log10(pval_cutoff)))+0.5)),
                                            width=grid::unit(3,"cm"),show_annotation_name = T,annotation_name_side="bottom",annotation_name_rot=0)
    fontsize=rep(6,nrow(combined_assay))
    fontsize[combined_neg_log_statistics>=(-log10(pval_cutoff))]<-10
    fontface=rep("plain",nrow(combined_assay))
    fontface[combined_neg_log_statistics>=(-log10(pval_cutoff))]<-"bold"
    #combined_activity_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(combined_activity))),name="zscores_activity",cluster_columns =cluster_between_groups(t(scale(t(delta_activity))),factor = patient_info[[group]]),cluster_rows = F,top_annotation = top_anno,column_split = length(unique(patient_info[[group]])),right_annotation = row_anno,row_names_gp = gpar(fontsize=fontsize,fontface=fontface),row_split = feature_splits)
    combined_assay_heatmap<-ComplexHeatmap::Heatmap(combined_assay,name=heatmap_name,col = col,cluster_rows = F,top_annotation = top_anno,column_split = groups,right_annotation = row_anno,row_names_gp = grid::gpar(fontsize=fontsize,fontface=fontface),row_split = feature_splits)
  }
  if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
  svg(filename = file.path(output_dir,"combined_assay_heatmap.svg"),width = 8,height=10)
  ComplexHeatmap::draw(combined_assay_heatmap)
  if(any(combined_neg_log_statistics>(-(log10(pval_cutoff))))){
    xpos<-0.95/(max(c(combined_neg_log_statistics,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
    for(i in 1:length(unique(feature_splits))){
      ComplexHeatmap::decorate_annotation("-log10(p)",{
        grid::grid.lines(c(xpos,xpos),c(0,1),gp = grid::gpar(lty = 2,col="red"))
      },slice=i)
    }
  }
  dev.off()
  combined_statistics=data.frame(Group=feature_splits,Feature=rep(feature_order,3),means_,pvalue=combined_statistics)
  return(list(combined_assay=combined_assay,combined_assay_heapmap=combined_assay_heatmap,combined_statistics=combined_statistics,group_terms=group_terms,groups=groups))
}



#' runAssayDeltaAnalysis
#'
#' run delta analysis about various assay
#'
#' @param obj Analysis. analysis object to be run.
#' @param assay_name character. assay name.
#' @param assay_type character. assay type, such as "DNA" or "RNA"
#' @param scale logic. whether to scale assay.
#' @param patient_id_column character. column for patient id.
#' @param sample_id_column character. column for sample id.
#' @param timepoint_column character. column for timepoint.
#' @param timepoint_terms character. timepoint terms pair as delta.
#' @param group_column character. column for group.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_patient_order character. patient orders for top annotation.
#' @param palette character. colors for KM curve.
#' @param test_method character. test method, such as "ttest" or "limma".
#' @param contrasts character. specific contrast.
#' @param pval_cutoff numeric. p value cutoff for significance.
#' @param run_go logic. whether run go analysis.
#' @param run_gsea logic. whether run gsea analysis.
#' @param logFC_cutoff numeric. log fold change cutoff for signigicant genes.
#' @param padj_cutoff numeric. adjusted p value cutoff for significance.
#' @param cancergenes character. cancer gene list.
#' @param pathwaylists list. pathway list.
#' @param specialpathwaylists list. special pathway(geneset) list.
#' @param only_sig_genes logic. whether only significant genes will be used.
#' @param only_show_sig_genes logic. whether only significant genes will be showed in heatmap.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runAssayDeltaAnalysis",function(obj,assay_name,assay_type="RNA",scale,
                                           patient_id_column="Patient_ID",sample_id_column="Sample_ID",timepoint_column="Timepoint",timepoint_terms=c("baseline"="Baseline","tp2"="TP2"),group_column,
                                           top_anno,top_anno_patient_order,palette,
                                           test_method=c("ttest","limma"),contrasts="SHORT-LONG",pval_cutoff = 0.05,
                                           run_go=FALSE,run_gsea=FALSE,
                                           logFC_cutoff = 1,padj_cutoff = 0.05,
                                           cancergenes=cancergenes, pathwaylists=pathwaylists,specialpathwaylists=specialpathwaylists,
                                           only_sig_genes=F,only_show_sig_genes=F,
                                           save_analysis = T,output_dir,
                                           log="run assay delta analysis.") standardGeneric("runAssayDeltaAnalysis"))
setMethod("runAssayDeltaAnalysis","Analysis",function(obj,assay_name,assay_type="RNA",scale,
                                                     patient_id_column="Patient_ID",sample_id_column="Sample_ID",timepoint_column="Timepoint",timepoint_terms=c("baseline"="Baseline","tp2"="TP2"),group_column,
                                                     top_anno,top_anno_patient_order,palette,
                                                     test_method=c("ttest","limma"),contrasts,pval_cutoff = 0.05,
                                                     run_go=FALSE,run_gsea=FALSE,
                                                     logFC_cutoff = 1,padj_cutoff = 0.05,
                                                     cancergenes=cancergenes, pathwaylists=pathwaylists,specialpathwaylists=specialpathwaylists,
                                                     only_sig_genes=F,only_show_sig_genes=F,
                                                     save_analysis = T,output_dir,
                                                     log="run assay delta analysis."){
  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  } else {
    output_dir<-getwd()
  }
  if(missing(scale)){scale=F}
  patients<-top_anno_patient_order
  patient_info<-merge(data.frame(Patient_ID=patients),obj@experiment@patient_info,by="Patient_ID",all.x=T)
  patient_info<-patient_info[match(patients,patient_info$Patient_ID),]
  sample_info<-obj@experiment@sample_info
  sample_info<-sample_info[sample_info$Assay==assay_type,]
  sample_info<-merge(sample_info,patient_info,by="Patient_ID",all.x=T)
  assay_data=methods::slot(obj@experiment,assay_name)@assay_data
  if(missing(palette)){palette=c("blue","white","red")}
  if(missing(contrasts)){paste(rev(sort(unique(sample_info[[group_column]]))),collapse = "-")}
  delta_analysis<-deltaAnalysis(assay=assay_data,scale=scale,
                                sample_info = sample_info,patient_id_column = patient_id_column,sample_id_column = sample_id_column,timepoint_column = timepoint_column,timepoint_terms = timepoint_terms,group_column = group_column,
                                top_anno = top_anno,top_anno_patient_order = top_anno_patient_order,palette = palette,
                                test_method = test_method,contrasts = contrasts,pval_cutoff = pval_cutoff,
                                only_sig_genes = only_sig_genes,only_show_sig_genes = only_show_sig_genes,
                                output_dir = output_dir)
  obj@output$combined_assay<-delta_analysis$combined_assay
  obj@output$combined_assay_heatmap<-delta_analysis$combined_assay_heatmap
  obj@output$combined_statistics<-delta_analysis$combined_statistics

  delta_statistics<-delta_analysis$combined_statistics[delta_analysis$combined_statistics$Group=="Delta",-1]
  delta_statistics[["logFC"]]<-delta_statistics[[delta_analysis[["group_terms"]][1]]]-delta_statistics[[delta_analysis[["group_terms"]][2]]]
  delta_statistics<-as.data.frame(delta_statistics)
  delta_statistics<-utiltools::set_column_as_rownames(delta_statistics,"Feature")
  delta_assay<-delta_analysis$combined_assay[delta_analysis$combined_statistics$Group=="Delta",]
  sig_delta_assay<-delta_assay[delta_statistics$pvalue<=0.05,]
  groups=delta_analysis$groups

  logFC_col="logFC"
  pval_col="pvalue"

  if(run_go) {
    go_output_dir<-file.path(output_dir,"GO")
    if(!dir.exists(go_output_dir)){dir.create(go_output_dir,recursive = T)}
    go<-utiltools::goEnrich(statistics = delta_statistics,pval_col = pval_col,logFC_col = logFC_col,pval_cutoff = pval_cutoff,logFC_cutoff = logFC_cutoff,padj_cutoff = padj_cutoff,output_dir = go_output_dir)
    obj@output$go<-go
  }
  if(run_gsea){
    for(pathwayset in names(pathwaylists)){
      pathways<-pathwaylists[[pathwayset]]
      specialpathways<-specialpathwaylists[[pathwayset]]
      pathwayset_output_dir<-file.path(output_dir,pathwayset)
      if(!dir.exists(pathwayset_output_dir)){dir.create(pathwayset_output_dir,recursive = T)}
      gsea_res<-utiltools::gsea(delta_statistics,pval_cutoff = pval_cutoff,FC_cutoff = logFC_cutoff,output_dir = pathwayset_output_dir,logFC_col = logFC_col,pval_col = pval_col,pathways = pathways)
      obj@output[[pathwayset]][["gsea"]]<-gsea_res
      selected_pathways<-pathways[intersect(names(gsea_res$sig_pathways),specialpathways)]
      if(length(selected_pathways)!=0){
        pathway_colors<-RColorBrewer::brewer.pal(length(selected_pathways),"Paired")[1:length(selected_pathways)]
        names(pathway_colors)<-names(selected_pathways)
        if(save_analysis){
          svg(filename = file.path(pathwayset_output_dir,"volcano_plot.svg"),width = 8,height=8)
          volcano_plot<-utiltools::ggvolcano(delta_statistics,x_col=logFC_col,y_col = pval_col,pathways =selected_pathways,pathway_colors = pathway_colors,cancergenes = cancergenes,ylab = "pvalue",FC_Cutoff = 1)
          obj@output[[pathwayset]][["volcano_plot"]]=volcano_plot
          dev.off()
        }

        pathway_gene_table_<-gsea_res$pathway_gene_table
        if(save_analysis){
          svg(filename = file.path(pathwayset_output_dir,"gseaheatmap.svg"),width = 10,height=8)
          gseaheatmap<-utiltools::gsea_heatmap(sig_expressions = sig_delta_assay,name="delta assay",scale = F,pathway_gene_table = pathway_gene_table_,pathway_font_size = 3,pathway_gene_heatmap_width = grid::unit(8,'cm'),top_annotation=top_anno,column_split=groups,show_column_dend=F,show_row_names=F,use_raster=T)
          obj@output[[pathwayset]][["gseaheatmap"]]=gseaheatmap
          dev.off()
        }
        if(save_analysis){
          svg(filename = file.path(pathwayset_output_dir,"sankeyheatmap.svg"),width = 10,height=8)
          sankeyheatmap<-utiltools::sankey_heatmap(sig_expressions = sig_delta_assay,name="delta assay",scale = F,keep_other = F,pathways = selected_pathways,pathway_colors = pathway_colors,top_annotation=top_anno,column_split=groups,show_column_dend=F,show_row_names=F,line_size = 3,text_size = 8)
          obj@output[[pathwayset]][["sankeyheatmap"]]=sankeyheatmap
          dev.off()
        }
      } else {
        if(save_analysis){
          svg(filename = file.path(pathwayset_output_dir,"volcano_plot.svg"),width = 8,height=8)
          volcano_plot<-utiltools::ggvolcano(delta_statistics,x_col=logFC_col,y_col = pval_col,cancergenes = cancergenes,ylab = "pvalue",FC_Cutoff = 1)
          obj@output[[pathwayset]][["volcano_plot"]]=volcano_plot
          dev.off()
        }
      }
    }
  }

  if(save_analysis){
    write.csv(obj@output$combined_assay,file.path(output_dir,"combined_assay.csv"),row.names = F)
    write.csv(obj@output$combined_statistics,file.path(output_dir,"combined_statistics.csv"),row.names = F)
  }

  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")

  return(obj)
}
)



#' difAnalysis
#'
#' differential assay anlysis
#'
#' @param assay data.frame. assay data.
#' @param scale logic. whether to scale assay.
#' @param assay_name character. assay name.
#' @param sample_info data.frame. sample information.
#' @param patient_id_column character. column for patient id.
#' @param sample_id_column character. column for sample id.
#' @param group_column character. column for group.
#' @param block_column character. column for block factor.
#' @param contrasts character. specific contrast.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param test_method character. test method, such as "ttest" or "limma".
#' @param pval_cutoff numeric. p value cutoff for significance.
#' @param only_sig_genes logic. whether only significant genes will be used.
#' @param palette character. colors used for heatmap.
#' @param output_dir character. output directory.
#' @param ... list. parameter passed tp ComplexHeatmap::Heatmap
#'
#' @return list. contain various contrast result.
#' @export
#'

difAnalysis<-function(assay,scale=F,assay_name="assay",
                      sample_info,patient_id_column="Patient_ID",sample_id_column="Sample_ID",group_column,block_column=NULL,contrasts=NULL,
                      top_anno,top_anno_sample_order,
                      test_method=c("ttest","limma"),pval_cutoff = 0.05,
                      only_sig_genes=F,palette=NULL,
                      output_dir,...){
  output<-list()
  if(scale){
    assay<-as.data.frame(t(scale(t(assay))))
  }
  test_method=match.arg(test_method)

  assay_samples<-top_anno_sample_order
  sample_info<-merge(data.frame(Sample_ID=assay_samples),sample_info,by=sample_id_column,all.x=T)
  sample_info<-sample_info[match(assay_samples, sample_info$Sample_ID),]

  assay<-as.data.frame(assay[,assay_samples])
  feature_order=rownames(assay)

  group_terms<-sort(unique(sample_info[[group_column]]))
  group_mean_<-list()
  for (group_term in group_terms){
    group_mean_[[paste(group_term,"_mean",sep="")]]<-rowMeans(assay[,sample_info[[group_column]]==group_term,drop=F],na.rm=T)
  }
  group_results<-as.data.frame(group_mean_)[feature_order,]

  if(is.null(contrasts)){
    contrasts=combinat::combn(group_terms,2,fun = function(items){paste(items,collapse="-")})
  }
  if(is.null(palette)){palette=c("blue","white","red")}
  col=circlize::colorRamp2(breaks=c(min(t(scale(t(assay)))),0,max(t(scale(t(assay))))),colors = palette)

  if(length(group_terms)>2){
    group_statistics=apply(assay,1,function(rowdata){df<-data.frame(value=as.numeric(rowdata),group=sample_info[[group_column]]);test=aov(value~group,data=df);return(summary(test)[[1]][["Pr(>F)"]][1])})
    group_statistics<-group_statistics[feature_order]
    group_results[["F_pvalue"]]<-group_statistics
    if(only_sig_genes){
      sig_genes<-feature_order[group_statistics<=pval_cutoff]
      f_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(assay[sig_genes,]))),name=paste("zscore",assay_name,sep="_"),col=col,column_split = sample_info[[group_column]],top_annotation = top_anno,...)
      if(!is.null(output_dir)){
        if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
        write.csv(group_results[sig_genes,],file=file.path(output_dir,"sig_F_Statistics.csv"))
        write.csv(assay[sig_genes,],file=file.path(output_dir,"sig_assay.csv"))
        svg(filename = file.path(output_dir,"F_heatmap.svg"),width = 20,height=12)
        ComplexHeatmap::draw(f_heatmap)
        dev.off()
      }
    } else {
      neg_log_statistics=-log10(group_statistics)
      row_anno<-ComplexHeatmap::rowAnnotation(`-log10(p)`=ComplexHeatmap::anno_barplot(neg_log_statistics,gp=grid::gpar(fill=ifelse(neg_log_statistics>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_statistics,-log10(pval_cutoff)))+0.5)),
                                              width=grid::unit(3,"cm"),show_annotation_name = T,annotation_name_side="bottom",annotation_name_rot=0)
      fontsize=rep(6,nrow(assay))
      fontsize[neg_log_statistics>=(-log10(pval_cutoff))]<-10
      fontface=rep("plain",nrow(assay))
      fontface[neg_log_statistics>=(-log10(pval_cutoff))]<-"bold"
      f_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(assay))),name=paste("zscore",assay_name,sep="_"),col = col,top_annotation = top_anno,column_split = sample_info[[group_column]],right_annotation = row_anno,row_names_gp = grid::gpar(fontsize=fontsize,fontface=fontface),...)
      if(!is.null(output_dir)){
        if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
        write.csv(group_results,file=file.path(output_dir,"F_Statistics.csv"))
        write.csv(assay,file=file.path(output_dir,"assay.csv"))
        svg(filename = file.path(output_dir,"F_heatmap.svg"),width = 20,height=12)
        ComplexHeatmap::draw(f_heatmap)
        if(any(neg_log_statistics>(-(log10(pval_cutoff))))){
          xpos<-0.95/(max(c(neg_log_statistics,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
          ComplexHeatmap::decorate_annotation("-log10(p)",{
            grid::grid.lines(c(xpos,xpos),c(0,1),gp = grid::gpar(lty = 2,col="red"))
          })
        }
        dev.off()
      }
    }
    output[["f_test"]]<-group_results
    output[["f_heatmamp"]]<-f_heatmap
  }
  for(contrast in contrasts){
    contrast_terms<-unlist(strsplit(contrast,"-"))
    contrast_results<-group_results[,paste(contrast_terms,"_mean",sep="")]
    contrast_results[["logFC"]]<-contrast_results[[paste(contrast_terms[1],"_mean",sep="")]]-contrast_results[[paste(contrast_terms[2],"_mean",sep="")]]
    subset_sample_info<-sample_info[sample_info[[group_column]] %in% contrast_terms,]
    subset_assay<-assay[,subset_sample_info[[sample_id_column]]]
    subset_top_anno<-top_anno[match(subset_sample_info[[sample_id_column]],assay_samples),]
    if(is.null(block_column)){
      if(test_method=="ttest"){
        subset_statistics<-apply(subset_assay,1,function(rowdata){data=data.frame(Value=as.numeric(rowdata),Group=subset_sample_info[[group_column]]);test=stats::t.test(Value~Group,data=data);return(test$p.value)})
        subset_statistics<-subset_statistics[feature_order]
      }
      if(test_method=="limma"){
        subset_statistics<-utiltools::dge_limma(subset_assay,is_rawcount = F,is_logged = T,normalize = F,sample_frequency_threshold = 0.5,clinic_info =data.frame(Sample_ID=colnames(subset_assay),Group=subset_sample_info[[group_column]]),ID_col = "Sample_ID",group_col = "Group",contrasts = contrast,method ="limma_trend")
        subset_statistics_<-subset_statistics$statistics[[colnames(subset_statistics$statistics)[grepl("p.value",colnames(subset_statistics$statistics),ignore.case=TRUE)]]]
        subset_statistics<-stats::setNames(subset_statistics_,rownames(subset_statistics$statistics))
        subset_statistics<-subset_statistics[feature_order]
      }
    } else {
      if(test_method=="ttest"){
        subset_statistics<-apply(subset_assay,1,function(rowdata){data=data.frame(Value=as.numeric(rowdata),Group=subset_sample_info[[group_column]],Block=subset_sample_info[[block_column]]);data=data[order(data$Group,data$Block),];test=stats::t.test(Value~Group,paired=T,data=data);return(test$p.value)})
        subset_statistics<-subset_statistics[feature_order]
      }
      if(test_method=="limma"){
        subset_statistics<-utiltools::dge_limma(subset_assay,is_rawcount = F,is_logged = T,normalize = F,sample_frequency_threshold = 0.5,clinic_info =data.frame(Sample_ID=colnames(subset_assay),Group=subset_sample_info[[group_column]],Block=subset_sample_info[[block_column]]),ID_col = "Sample_ID",group_col = "Group",block_col = "Block",contrasts = contrast,method ="limma_trend")
        ubset_statistics_<-subset_statistics$statistics[[colnames(subset_statistics$statistics)[grepl("p.value",colnames(subset_statistics$statistics),ignore.case=TRUE)]]]
        subset_statistics<-stats::setNames(subset_statistics_,rownames(subset_statistics$statistics))
        subset_statistics<-subset_statistics[feature_order]
      }
    }

    if(only_sig_genes){
      sig_genes<-feature_order[subset_statistics<=pval_cutoff]
      contrast_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(subset_assay[sig_genes,]))),name=paste("zscore",assay_name,sep="_"),column_split = subset_sample_info[[group_column]],top_annotation = subset_top_anno,...)
      if(!is.null(output_dir)){
        if(!dir.exists(file.path(output_dir,contrast))){dir.create(file.path(output_dir,contrast),recursive = T)}
        svg(filename = file.path(output_dir,contrast,paste(contrast,"_heatmap.svg",sep="")),width = 20,height=12)
        ComplexHeatmap::draw(contrast_heatmap)
        dev.off()
      }
    } else {
      neg_log_statistics=-log10(subset_statistics)
      row_anno<-ComplexHeatmap::rowAnnotation(`-log10(p)`=ComplexHeatmap::anno_barplot(neg_log_statistics,gp=grid::gpar(fill=ifelse(neg_log_statistics>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_statistics,-log10(pval_cutoff)))+0.5)),
                                              width=grid::unit(3,"cm"),show_annotation_name = T,annotation_name_side="bottom",annotation_name_rot=0)
      fontsize=rep(6,nrow(subset_assay))
      fontsize[neg_log_statistics>=(-log10(pval_cutoff))]<-10
      fontface=rep("plain",nrow(subset_assay))
      fontface[neg_log_statistics>=(-log10(pval_cutoff))]<-"bold"
      contrast_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(subset_assay))),name=paste("zscore",assay_name,sep="_"),col = col,top_annotation = subset_top_anno,column_split = subset_sample_info[[group_column]],right_annotation = row_anno,row_names_gp = grid::gpar(fontsize=fontsize,fontface=fontface),...)
      if(!is.null(output_dir)){
        if(!dir.exists(file.path(output_dir,contrast))){dir.create(file.path(output_dir,contrast),recursive = T)}
        svg(filename = file.path(output_dir,contrast,paste(contrast,"_heatmap.svg",sep="")),width = 20,height=12)
        ComplexHeatmap::draw(contrast_heatmap)
        if(any(neg_log_statistics>(-(log10(pval_cutoff))))){
          xpos<-0.95/(max(c(neg_log_statistics,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
          ComplexHeatmap::decorate_annotation("-log10(p)",{
            grid::grid.lines(c(xpos,xpos),c(0,1),gp = grid::gpar(lty = 2,col="red"))
          })
        }
        dev.off()
      }
    }
    contrast_assay=subset_assay
    contrast_group=subset_sample_info[[group_column]]
    contrast_results[["pvalue"]]<-subset_statistics
    contrast_results[["adj_pvalue"]]<-stats::p.adjust(subset_statistics,method = "BH")

    output[[contrast]]<-list(contrast_assay=contrast_assay,contrast_group=contrast_group,contrast_results=contrast_results,contrast_heatmap=contrast_heatmap)
  }

  return(output)
}


#' runAssayDifAnalysis
#'
#' run delta analysis about various assay
#'
#' @param obj Analysis. analysis object to be run.
#' @param assay_name character. assay name.
#' @param assay_type character. assay type, such as "DNA" or "RNA"
#' @param scale logic. whether to scale assay.
#' @param patient_id_column character. column for patient id.
#' @param sample_id_column character. column for sample id.
#' @param group_column character. column for group.
#' @param block_column character. column for block.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. patient orders for top annotation.
#' @param palette character. colors for KM curve.
#' @param test_method character. test method, such as "ttest" or "limma".
#' @param contrasts character. specific contrast.
#' @param pval_cutoff numeric. p value cutoff for significance.
#' @param run_go logic. whether run go analysis.
#' @param run_gsea logic. whether run gsea analysis.
#' @param logFC_cutoff numeric. log fold change cutoff for signigicant genes.
#' @param padj_cutoff numeric. adjusted p value cutoff for significance.
#' @param cancergenes character. cancer gene list.
#' @param pathwaylists list. pathway list.
#' @param specialpathwaylists list. special pathway(geneset) list.
#' @param only_sig_genes logic. whether only significant genes will be used.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if is.null, use project/analysis_name as directory
#' @param log character. any comments.
#' @param ... list. parameter passed to difAnalysis.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runAssayDifAnalysis",function(obj,assay_name,assay_type="RNA",scale,
                                          patient_id_column="Patient_ID",sample_id_column="Sample_ID",group_column,block_column=NULL,
                                          top_anno,top_anno_sample_order,palette=NULL,
                                          test_method=c("ttest","limma"),contrasts=NULL,pval_cutoff = 0.05,
                                          run_go=FALSE,run_gsea=FALSE,
                                          logFC_cutoff = 1,padj_cutoff = 0.05,
                                          cancergenes=cancergenes, pathwaylists=pathwaylists,specialpathwaylists=specialpathwaylists,
                                          only_sig_genes=F,
                                          save_analysis = T,output_dir,
                                          log="run assay differential analysis.",...) standardGeneric("runAssayDifAnalysis"))
setMethod("runAssayDifAnalysis","Analysis",function(obj,assay_name,assay_type="RNA",scale,
                                                    patient_id_column="Patient_ID",sample_id_column="Sample_ID",group_column,block_column=NULL,
                                                    top_anno,top_anno_sample_order,palette=NULL,
                                                    test_method=c("ttest","limma"),contrasts=NULL,pval_cutoff = 0.05,
                                                    run_go=FALSE,run_gsea=FALSE,
                                                    logFC_cutoff = 1,padj_cutoff = 0.05,
                                                    cancergenes=cancergenes, pathwaylists=pathwaylists,specialpathwaylists=specialpathwaylists,
                                                    only_sig_genes=F,
                                                    save_analysis = T,output_dir,
                                                    log="run assay differential analysis.",...){
  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  } else {
    output_dir<-getwd()
  }
  if(missing(scale)){scale=F}
  assay_samples<-top_anno_sample_order
  sample_info<-obj@experiment@sample_info
  sample_info<-sample_info[sample_info$Assay==assay_type,]
  sample_info<-merge(sample_info,obj@experiment@patient_info,by=patient_id_column,all.x=T)
  sample_info<-sample_info[match(assay_samples,sample_info[[sample_id_column]]),]

  assay_data=methods::slot(obj@experiment,assay_name)@assay_data[,assay_samples]

  if(length(group_column)==1){
    profile_group_<-sample_info[[group_column]]
    if(any(is.na(profile_group_) | profile_group_=="")){
      excluded_samples<-sample_info[["Sample_ID"]][is.na(profile_group_) | profile_group_==""]
      sample_info=sample_info[-match(excluded_samples,sample_info[["Sample_ID"]]),]
      assay_samples=assay_samples[-match(excluded_samples,assay_samples)]
      top_anno<-top_anno[match(assay_samples,top_anno_sample_order)]
      assertthat::are_equal(sample_info[["Sample_ID"]],assay_samples)
      assay_data=assay_data[,assay_samples]
    }
  }
  if(missing(palette)){palette=c("blue","white","red")}

  if(is.null(contrasts)){
    group_terms<-sort(unique(sample_info[[group_column]]))
    contrasts=combinat::combn(group_terms,2,fun = function(items){paste(items,collapse="-")})
  }

  dif_analysis<-difAnalysis(assay=assay_data,scale=scale,assay_name=assay_name,
                            sample_info=sample_info,patient_id_column=patient_id_column,sample_id_column=sample_id_column,group_column=group_column,block_column=block_column,contrasts=contrasts,
                            top_anno=top_anno,top_anno_sample_order=assay_samples,
                            test_method=test_method,pval_cutoff = pval_cutoff,
                            only_sig_genes=only_sig_genes,palette=palette,
                            output_dir=output_dir,...)

  obj@output$dif_analysis<-dif_analysis


  logFC_col="logFC"
  pval_col="pvalue"
  for(contrast in contrasts){

    contrast_output_dir<-file.path(output_dir,contrast)
    contrast_assay=dif_analysis[[contrast]][["contrast_assay"]]
    groups=dif_analysis[[contrast]][["contrast_group"]]
    contrast_statistics=dif_analysis[[contrast]][["contrast_results"]]
    sig_assay=contrast_assay[contrast_statistics[[pval_col]]<=pval_cutoff,]
    top_anno_<-top_anno[match(colnames(sig_assay),assay_samples)]
    obj@output[[contrast]][["contrast_group"]]=groups
    obj@output[[contrast]][["contrast_statistics"]]=contrast_statistics
    obj@output[[contrast]][["sig_assay"]]=sig_assay
    if(run_go) {
      go_output_dir<-file.path(contrast_output_dir,"GO")
      if(!dir.exists(go_output_dir)){dir.create(go_output_dir,recursive = T)}
      go<-utiltools::goEnrich(statistics = contrast_statistics,pval_col = pval_col,logFC_col = logFC_col,pval_cutoff = pval_cutoff,logFC_cutoff = logFC_cutoff,padj_cutoff = padj_cutoff,output_dir = go_output_dir)
      obj@output[[contrast]]$go<-go
    }

    if(run_gsea){
      for(pathwayset in names(pathwaylists)){
        pathways<-pathwaylists[[pathwayset]]
        specialpathways<-specialpathwaylists[[pathwayset]]
        pathwayset_output_dir<-file.path(contrast_output_dir,pathwayset)
        if(!dir.exists(pathwayset_output_dir)){dir.create(pathwayset_output_dir,recursive = T)}
        gsea_res<-utiltools::gsea(contrast_statistics,pval_cutoff = pval_cutoff,FC_cutoff = logFC_cutoff,output_dir = pathwayset_output_dir,logFC_col = logFC_col,pval_col = pval_col,pathways = pathways)
        obj@output[[contrast]][[pathwayset]][["gsea"]]<-gsea_res
        selected_pathways<-pathways[intersect(names(gsea_res$sig_pathways),specialpathways)]
        if(length(selected_pathways)!=0){
          pathway_colors<-RColorBrewer::brewer.pal(length(selected_pathways),"Paired")[1:length(selected_pathways)]
          names(pathway_colors)<-names(selected_pathways)
          if(save_analysis){
            svg(filename = file.path(pathwayset_output_dir,"volcano_plot.svg"),width = 10,height=10)
            volcano_plot<-utiltools::ggvolcano(contrast_statistics,x_col=logFC_col,y_col = pval_col,pathways =selected_pathways,pathway_colors = pathway_colors,cancergenes = cancergenes,ylab = "pvalue",FC_Cutoff = 1)
            obj@output[[contrast]][[pathwayset]][["volcano_plot"]]=volcano_plot
            dev.off()
          }

          pathway_gene_table_<-gsea_res$pathway_gene_table
          if(save_analysis){
            svg(filename = file.path(pathwayset_output_dir,"gseaheatmap.svg"),width = 20,height=12)
            gseaheatmap<-utiltools::gsea_heatmap(sig_expressions = t(scale(t(sig_assay))),name="differential assay",scale = F,pathway_gene_table = pathway_gene_table_,pathway_font_size = 3,pathway_gene_heatmap_width = grid::unit(8,'cm'),top_annotation=top_anno_,column_split=groups,show_column_dend=F,show_row_names=F,use_raster=T)
            obj@output[[contrast]][[pathwayset]][["gseaheatmap"]]=gseaheatmap
            dev.off()
          }
          if(save_analysis){
            svg(filename = file.path(pathwayset_output_dir,"sankeyheatmap.svg"),width = 20,height=12)
            sankeyheatmap<-utiltools::sankey_heatmap(sig_expressions = t(scale(t(sig_assay))),name="differential assay",scale = F,keep_other = F,pathways = selected_pathways,pathway_colors = pathway_colors,top_annotation=top_anno_,column_split=groups,show_column_dend=F,show_row_names=F,line_size = 3,text_size = 8)
            obj@output[[contrast]][[pathwayset]][["sankeyheatmap"]]=sankeyheatmap
            dev.off()
          }
        } else {
          if(save_analysis){
            svg(filename = file.path(pathwayset_output_dir,"volcano_plot.svg"),width = 10,height=10)
            volcano_plot<-utiltools::ggvolcano(contrast_statistics,x_col=logFC_col,y_col = pval_col,cancergenes = cancergenes,ylab = "pvalue",FC_Cutoff = 1)
            obj@output[[contrast]][[pathwayset]][["volcano_plot"]]=volcano_plot
            dev.off()
          }
        }
      }
    }
    if(save_analysis){
      if(!dir.exists(contrast_output_dir)){dir.create(contrast_output_dir,recursive = T)}
      write.csv(sig_assay,file.path(contrast_output_dir,"sig_assay.csv"),row.names = T)
      write.csv(contrast_statistics,file.path(contrast_output_dir,"combined_statistics.csv"),row.names = T)
    }
  }

  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")

  return(obj)
}
)


#' runAssayunSupervisedAnalysis
#'
#' run unsupervised assay anlysis
#'
#' @param obj analysis. an analysis object.
#' @param assay_name character. assay name.
#' @param assay_type character. type of assay
#' @param scale logic. whether to scale assay.
#' @param stablecluster_method character. stable cluster method.
#' @param cutFUN function. cluster function.
#' @param clusters numeric. number of cluters.
#' @param nTimes numeric. number of permutation.
#' @param subSampleSize numeric. ratio of samples for each permutation.
#' @param subFeatureSize numetic. ratio of features for each permutation.
#' @param patient_id_column character. column for patient id.
#' @param sample_id_column character. column for sample id.
#' @param survival_column character. column for survival.
#' @param status_column character. column for status.
#' @param survival_type character. type of survival analysis.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. patient orders for top annotation.
#' @param palette character. colors used for heatmap.
#' @param test_method character. test method, such as "ttest" or "limma".
#' @param contrasts character. contrast for test.
#' @param pval_cutoff numeric. p value cutoff for significance.
#' @param run_go logic. whether run go analysis.
#' @param run_gsea logic. whether run gsea analysis.
#' @param logFC_cutoff numeric. log fold change cutoff for signigicant genes.
#' @param padj_cutoff numeric. adjusted p value cutoff for significance.
#' @param cancergenes character. cancer gene list.
#' @param pathwaylists list. pathway list.
#' @param specialpathwaylists list. special pathway(geneset) list.
#' @param only_sig_genes logic. whether only significant genes will be used.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory.
#' @param log chracter. any comments.
#' @param ... list. parameter passed tp difAnalysis
#'
#' @return analysis. contain various contrast result.
#' @export
#'

setGeneric("runAssayunSupervisedAnalysis",function(obj,assay_name,assay_type="RNA",scale=F,
                                                   stablecluster_method="combine",cutFUN=c(ClassDiscovery::cutHclust),clusters=NULL,nTimes=100,subSampleSize = 0.9,subFeatureSize=0.9,
                                                   patient_id_column="Patient_ID",sample_id_column="Sample_ID",survival_column=NULL,status_column=NULL,survival_type=c("OS","PFS"),
                                                   top_anno,top_anno_sample_order,palette=NULL,
                                                   test_method=c("ttest","limma"),contrasts=NULL,pval_cutoff = 0.05,
                                                   run_go=FALSE,run_gsea=FALSE,
                                                   logFC_cutoff = 1,padj_cutoff = 0.05,
                                                   cancergenes=cancergenes, pathwaylists=pathwaylists,specialpathwaylists=specialpathwaylists,
                                                   only_sig_genes=F,
                                                   save_analysis = T,output_dir,log="run unsupervised assay analysis.",...) standardGeneric("runAssayunSupervisedAnalysis"))
setMethod("runAssayunSupervisedAnalysis","Analysis",function(obj,assay_name,assay_type="RNA",scale=F,
                                                             stablecluster_method="combine",cutFUN=c(ClassDiscovery::cutHclust),clusters=NULL,nTimes=100,subSampleSize = 0.9,subFeatureSize=0.9,
                                                             patient_id_column="Patient_ID",sample_id_column="Sample_ID",survival_column=NULL,status_column=NULL,survival_type=c("OS","PFS"),
                                                             top_anno,top_anno_sample_order,palette=NULL,
                                                             test_method=c("ttest","limma"),contrasts=NULL,pval_cutoff = 0.05,
                                                             run_go=FALSE,run_gsea=FALSE,
                                                             logFC_cutoff = 1,padj_cutoff = 0.05,
                                                             cancergenes=cancergenes, pathwaylists=pathwaylists,specialpathwaylists=specialpathwaylists,
                                                             only_sig_genes=F,
                                                             save_analysis = T,output_dir,log="run unsupervised assay analysis.",...){

  test_method=match.arg(test_method)
  survival_type=match.arg(survival_type)

  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  } else {
    output_dir<-getwd()
  }

  assay_samples<-top_anno_sample_order
  sample_info<-obj@experiment@sample_info
  sample_info<-sample_info[sample_info$Assay==assay_type,]
  sample_info<-merge(sample_info,obj@experiment@patient_info,by=patient_id_column,all.x=T)
  if((!is.null(survival_column)) & (!is.null(status_column))){
    colnames(sample_info)[match(c(survival_column,status_column),colnames(sample_info))]<-c(survival_type,"Status")
  }
  sample_info<-sample_info[match(assay_samples,sample_info[[sample_id_column]]),]

  assay_data=methods::slot(obj@experiment,assay_name)@assay_data[,assay_samples]
  if(scale){
    assay_data<-as.data.frame(t(scale(t(assay_data))))
    if(is.null(palette)){palette=c("blue","white","red")}
    col=circlize::colorRamp2(breaks=c(min(assay_data),0,max(assay_data)),colors = palette)
  } else {
    if(is.null(palette)){palette=c("blue","white","red")}
    col=circlize::colorRamp2(breaks=c(min(t(scale(t(assay_data)))),0,max(t(scale(t(assay_data))))),colors = palette)

  }

  if(nrow(assay_data)>=300){
    unsupervised_assay<-utiltools::prepare_unsupervised_data(assay_data,method="GUMBEL",gumbel_p_cutoff = 0.1,remove_outlier = T)
  } else {
    unsupervised_assay<-assay_data
  }

  if(is.null(clusters)){
    clusters<-utiltools::estimate_bestNumberofClusters(unsupervised_assay,max.nc = as.integer(ncol(unsupervised_assay)/2))$Best.NumberofCluster
  }

  stableclusters<-utiltools::stableCluster(unsupervised_assay,method = stablecluster_method,cutFUN = cutFUN,clusters = clusters,nTimes = nTimes,subSampleSize = subSampleSize,subFeatureSize = subFeatureSize)
  assertthat::are_equal(colnames(stableclusters),colnames(unsupervised_assay))

  stablecluster_assay_heatmap<-ComplexHeatmap::draw(ComplexHeatmap::Heatmap(t(scale(t(unsupervised_assay))),name=paste("zscores_",assay_name,sep=""),col=col,cluster_columns = stats::hclust(stats::dist(stableclusters)),top_annotation = top_anno,column_split = clusters,...))
  clusters_=ComplexHeatmap::column_order(stablecluster_assay_heatmap)
  names(clusters_)<-paste("Cluster_",1:length(clusters_),sep="")
  clusters_<-utils::stack(clusters_)
  clusters_<-clusters_[order(clusters_[["values"]]),]
  sample_info[["Cluster"]]<-as.character(clusters_[["ind"]])


  svg(filename = file.path(output_dir,"stablecluster_assay_heatmap.svg"),width = 20,height=12)
  ComplexHeatmap::draw(stablecluster_assay_heatmap)
  dev.off()


  stablecluster_assay_data<-unsupervised_assay[ComplexHeatmap::row_order(stablecluster_assay_heatmap),unlist(ComplexHeatmap::column_order(stablecluster_assay_heatmap))]
  obj@output[["stablecluster_assay_heatmap"]]=stablecluster_assay_heatmap
  obj@output[["stablecluster_assay_data"]]=stablecluster_assay_data

  #stablecluster survival anaysis

  km_palette=RColorBrewer::brewer.pal(n=max(c(3,clusters)),name="Paired")[1:clusters]
  assertthat::are_equal(sample_info[[sample_id_column]],colnames(unsupervised_assay))
  if((!is.null(survival_column)) & (!is.null(status_column))){
    cluster_km_plot<-utiltools::hCluster_Surv(heatmap=stablecluster_assay_heatmap,clinics=sample_info,survival_col = survival_type,status_col = "Status",palette = km_palette,legend=c(0.8,0.2))
    svg(filename = file.path(output_dir,"cluster_km_plot.svg"),width = 10,height=10)
    print(cluster_km_plot)
    dev.off()
    obj@output[["cluster_km_plotp"]]=cluster_km_plot
  }

  group_column="Cluster"
  if(missing(palette)){palette=c("blue","white","red")}
  col=circlize::colorRamp2(breaks=c(min(assay_data),0,max(assay_data)),colors = palette)
  group_terms<-sort(unique(sample_info[[group_column]]))
  contrasts=combinat::combn(group_terms,2,fun = function(items){paste(items,collapse="-")})

  dif_analysis<-difAnalysis(assay=assay_data,scale=scale,assay_name=assay_name,
                            sample_info=sample_info,patient_id_column=patient_id_column,sample_id_column=sample_id_column,group_column=group_column,contrasts = contrasts,
                            top_anno=top_anno,top_anno_sample_order=top_anno_sample_order,
                            test_method=test_method,pval_cutoff = pval_cutoff,
                            only_sig_genes=only_sig_genes,palette=palette,
                            output_dir=output_dir,...)

  obj@output[["dif_analysis"]]<-dif_analysis


  logFC_col="logFC"
  pval_col="pvalue"
  for(contrast in contrasts){
    contrast_output_dir<-file.path(output_dir,contrast)
    if(!dir.exists(contrast_output_dir)){dir.create(contrast_output_dir)}
    contrast_assay=dif_analysis[[contrast]][["contrast_assay"]]
    groups=dif_analysis[[contrast]][["contrast_group"]]
    contrast_statistics=dif_analysis[[contrast]][["contrast_results"]]
    sig_assay=contrast_assay[contrast_statistics[[pval_col]]<=pval_cutoff,]
    obj@output[[contrast]][["contrast_group"]]=groups
    obj@output[[contrast]][["contrast_statistics"]]=contrast_statistics
    obj@output[[contrast]][["sig_assay"]]=sig_assay

    if(save_analysis){
      write.csv(contrast_statistics,file.path(contrast_output_dir,"statistics.csv"))
      write.csv(sig_assay,file.path(contrast_output_dir,"sig_assay.csv"))
    }
    if(run_go) {
      go_output_dir<-file.path(contrast_output_dir,"GO")
      if(!dir.exists(go_output_dir)){dir.create(go_output_dir,recursive = T)}
      go<-utiltools::goEnrich(statistics = contrast_statistics,pval_col = pval_col,logFC_col = logFC_col,pval_cutoff = pval_cutoff,logFC_cutoff = logFC_cutoff,padj_cutoff = padj_cutoff,output_dir = go_output_dir)
      obj@output[[contrast]]$go<-go
    }
    if(run_gsea){
      for(pathwayset in names(pathwaylists)){
        pathways<-pathwaylists[[pathwayset]]
        specialpathways<-specialpathwaylists[[pathwayset]]
        pathwayset_output_dir<-file.path(contrast_output_dir,pathwayset)
        if(!dir.exists(pathwayset_output_dir)){dir.create(pathwayset_output_dir,recursive = T)}
        gsea_res<-utiltools::gsea(contrast_statistics,pval_cutoff = pval_cutoff,FC_cutoff = logFC_cutoff,output_dir = pathwayset_output_dir,logFC_col = logFC_col,pval_col = pval_col,pathways = pathways)
        obj@output[[contrast]][[pathwayset]][["gsea"]]<-gsea_res
        selected_pathways<-pathways[intersect(names(gsea_res$sig_pathways),specialpathways)]
        if(length(selected_pathways)!=0){
          pathway_colors<-RColorBrewer::brewer.pal(length(selected_pathways),"Paired")[1:length(selected_pathways)]
          names(pathway_colors)<-names(selected_pathways)
          if(save_analysis){
            svg(filename = file.path(pathwayset_output_dir,"volcano_plot.svg"),width = 10,height=8)
            volcano_plot<-utiltools::ggvolcano(contrast_statistics,x_col=logFC_col,y_col = pval_col,pathways =selected_pathways,pathway_colors = pathway_colors,cancergenes = cancergenes,ylab = "pvalue",FC_Cutoff = 1)
            obj@output[[contrast]][[pathwayset]][["volcano_plot"]]=volcano_plot
            dev.off()
          }

          pathway_gene_table_<-gsea_res$pathway_gene_table
          if(save_analysis){
            svg(filename = file.path(pathwayset_output_dir,"gseaheatmap.svg"),width = 20,height=12)
            gseaheatmap<-utiltools::gsea_heatmap(sig_expressions = t(scale(t(sig_assay))),name="differential assay",scale = F,pathway_gene_table = pathway_gene_table_,pathway_font_size = 3,pathway_gene_heatmap_width = grid::unit(8,'cm'),top_annotation=top_anno,column_split=groups,show_column_dend=F,show_row_names=F,use_raster=T)
            obj@output[[contrast]][[pathwayset]][["gseaheatmap"]]=gseaheatmap
            dev.off()
          }
          if(save_analysis){
            svg(filename = file.path(pathwayset_output_dir,"sankeyheatmap.svg"),width = 20,height=12)
            sankeyheatmap<-utiltools::sankey_heatmap(sig_expressions = t(scale(t(sig_assay))),name="differential assay",scale = F,keep_other = F,pathways = selected_pathways,pathway_colors = pathway_colors,top_annotation=top_anno,column_split=groups,show_column_dend=F,show_row_names=F,line_size = 3,text_size = 8)
            obj@output[[contrast]][[pathwayset]][["sankeyheatmap"]]=sankeyheatmap
            dev.off()
          }
        } else {
          if(save_analysis){
            svg(filename = file.path(pathwayset_output_dir,"volcano_plot.svg"),width = 10,height=8)
            volcano_plot<-utiltools::ggvolcano(contrast_statistics,x_col=logFC_col,y_col = pval_col,cancergenes = cancergenes,ylab = "pvalue",FC_Cutoff = 1)
            obj@output[[contrast]][[pathwayset]][["volcano_plot"]]=volcano_plot
            dev.off()
          }
        }
      }
    }

  }

  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  if(save_analysis){
    write.csv(sample_info,file.path(output_dir,"sample_info.csv"))
    cat(obj@log,file=file.path(output_dir,"log.txt"))
  }
  return(obj)
}
)


#' runFusionAnalysis
#' run analysis of fusions
#'
#' @param obj analysis. an analysis object.
#' @param top_fusions numeric. top fusions number.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. patient orders for top annotation.
#' @param type_col character. colors for fusion type.
#' @param sample_order character.parameter passed to oncoprint.
#' @param annotate_frequence logic. whether mark number of fusion detection method.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory.
#' @param log chracter. any comments.
#' @param ... list. parameter passed oncoprint.
#'
#' @return analysis. contain various contrast result.
#' @export
#'
setGeneric("runFusionAnalysis",function(obj,top_fusions=30L,top_anno,top_anno_sample_order,type_col=NULL,sample_order=NULL,annotate_frequence=T,save_analysis=F,output_dir,log="run fusion analysis.",...) standardGeneric("runFusionAnalysis"))
setMethod("runFusionAnalysis","Analysis",function(obj,top_fusions=30L,top_anno,top_anno_sample_order,type_col=NULL,sample_order=NULL,annotate_frequence=T,save_analysis=F,output_dir,log="run fusion analysis.",...){
  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  }
  RNA_samples<-top_anno_sample_order
  sample_info<-merge(merge(data.frame(Sample_ID=RNA_samples),obj@experiment@sample_info,by="Sample_ID",all.x=T),obj@experiment@patient_info,by="Patient_ID",all.x=T)
  sample_info$Tumor_Sample_Barcode<-sample_info$Sample_ID

  fusions<-as.data.frame(obj@experiment@fusion_assay@assay_data)
  fusions<-fusions[fusions[["Sample_ID"]] %in% RNA_samples,]
  fusion_samples<-unique(fusions$Sample_ID)
  RNA_samples<-fusion_samples
  top_anno<-top_anno[match(RNA_samples,top_anno_sample_order)]

  if(!is.null(sample_order)){
    assertthat::are_equal(sort(RNA_samples),sort(sample_order))
  }
  fp<-utiltools::fusionPlot(fusions=fusions,top_filter=top_fusions,top_anno = top_anno,top_anno_sample_order = RNA_samples,type_col=type_col,sample_order=sample_order,annotate_frequence=annotate_frequence,output_dir=output_dir,...)

  obj@output$fusion_plot<-fp
  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  if(save_analysis){
    write.csv(fp$fusions_svtype_matrix,file.path(output_dir,"fusions_svtype_matrix.csv"))
    write.csv(fp$fusions_frequence_matrix,file.path(output_dir,"fusions_frequence_matrix.csv"))
    write.csv(fp$fusions_caller_matrix,file.path(output_dir,"fusions_caller_matrix.csv"))
    write.csv(fp$ffusions_genefunction_matrix,file.path(output_dir,"fusions_genefunction_matrix.csv"))
    cat(obj@log,file=file.path(output_dir,"log.txt"))
  }
  return(obj)
}
)


#' runCloneAnalysis
#' run pyclone Analysis
#'
#' @param obj analysis. an analysis object.
#' @param mutect_snps data.frame. snps derived from mutect.
#' @param sequenza_snps data.frame. snps derived from sequenza.
#' @param sequenza_segments data.frame. segments derived from sequenza.
#' @param purities data.frame. tumor purity derived from sequenza.
#' @param overlapping_ref_counts_cutoff numeric. ref count cutoff for overlapping mutations(Baseline and TP2).
#' @param overlapping_alt_counts_cutoff numeric. alt count cutoff for overlapping mutations(Baseline and TP2).
#' @param nonoverlapping_ref_counts_cutoff numeric. ref count cutoff for nonoverlapping mutations(Baseline and TP2).
#' @param nonoverlapping_alt_counts_cutoff numeric. alt count cutoff for nonoverlapping mutations(Baseline and TP2).
#' @param gene_locs data.frame. gene infomation
#' @param cancer_genes character. cancer gene list.
#' @param smg_genes character. SMG gene list.
#' @param level character. arm level, cytoband level.
#' @param pyclone_path character. path to pyclone
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory.
#' @param log chracter. any comments.
#'
#' @return analysis. contain various contrast result.
#' @export
#'
setGeneric("runCloneAnalysis",function(obj,mutect_snps,sequenza_snps=NULL,sequenza_segments,purities,overlapping_ref_counts_cutoff=20,overlapping_alt_counts_cutoff=3,nonoverlapping_ref_counts_cutoff=100,nonoverlapping_alt_counts_cutoff=10,gene_locs,cancer_genes,smg_genes,level="arm",pyclone_path="/home/harryjerry/anaconda3/bin/pyclone-vi",save_analysis=F,output_dir,log="run sequenza clone analysis.") standardGeneric("runCloneAnalysis"))
setMethod("runCloneAnalysis","Analysis",function(obj,mutect_snps,sequenza_snps=NULL,sequenza_segments,purities,overlapping_ref_counts_cutoff=20,overlapping_alt_counts_cutoff=3,nonoverlapping_ref_counts_cutoff=100,nonoverlapping_alt_counts_cutoff=10,gene_locs,cancer_genes,smg_genes,level="arm",pyclone_path="/home/harryjerry/anaconda3/bin/pyclone-vi",save_analysis=F,output_dir,log="run sequenza clone analysis."){
  f<-function(l) {
    if(length(l)==1){
      return(as.character(l))
    } else{
      return(paste(l,collapse=";"))
    }}
  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(getwd(),obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  }

  if(!missing(gene_locs)){
    gene_locs<-utiltools::df2granges(gene_locs,xy=T,seqnames_col = "chromosome",meta_cols=c("Gene_name"))
  }

  if(missing(purities)){ purities<-obj@experiment@purity_assay@assay_data}
  sample_timepoint_info<-obj@experiment@sample_info %>% dplyr::filter(Assay=='DNA') %>% dplyr::select(c("Sample_ID","Patient_ID","Timepoint")) %>% tidyr::pivot_wider(id_cols = Patient_ID,names_from = Timepoint,values_from=Sample_ID,values_fn = f)
  sample_timepoint_info<-sample_timepoint_info[rowSums(!is.na(sample_timepoint_info[,-match("Patient_ID",colnames(sample_timepoint_info))]))>=2,]
  if(!is.null(sequenza_snps)){
    sequenza_snps<-utiltoools::df2granges(sequenza_snps,xy=T,seqnames_col = "chrom",meta_cols = c("Sample_ID","ref_counts","alt_counts"))
    sequenza_snps<-plyranges::join_overlap_left(sequenza_snps,gene_locs)
    sequenza_snps<-as.data.frame(sequenza_snps)
    sequenza_snps<-with(sequenza_snps,data.frame(Sample_ID=Sample_ID,
                                                 #Patient_ID=gsub("_[a-zA-Z0-9]+-.*","",Sample_ID),
                                                 Mutation_ID=paste(seqnames,start,sep=":"),
                                                 chrom=seqnames,
                                                 start=start,
                                                 end=end,
                                                 gene=Gene_name,
                                                 ref_counts=ref_counts,
                                                 alt_counts=alt_counts))
    sequenza_snps<-sequenza_snps[!duplicated(sequenza_snps[,c("Sample_ID","Mutation_ID")]),]
    #sequenza_snps_<-sequenza_snps_ %>% dplyr::group_by(Patient_ID,Mutation_ID) %>% dplyr::mutate(n_sample_inpatient=n(),patient_alt_counts=sum(alt_counts))

    #n_sample_inpatient<-pivot_wider(sequenza_snps_,id_cols ="Mutation_ID",names_from = "Patient_ID",values_from = "n_sample_inpatient",values_fn = function(i){unique(i)},values_fill = 0)

    #patient_alt_counts<-pivot_wider(sequenza_snps_,id_cols ="Mutation_ID",names_from = "Patient_ID",values_from = "alt_counts",values_fn = function(i){sum(i)},values_fill = 0)
    #assertthat::are_equal(n_sample_inpatient$Mutation_ID,patient_alt_counts$Mutation_ID)

    #sequenza_snps_wh_genename<-sequenza_snps_[!is.na(sequenza_snps_$gene),]
    #sequenza_snps_wh_genename<- sequenza_snps_wh_genename %>% dplyr::group_by(Patient_ID,gene) %>% filter(patient_alt_counts==max(patient_alt_counts))
    #sequenza_snps_wo_genename<-sequenza_snps_[is.na(sequenza_snps_$gene),]

  }
  for(i in 1:nrow(sample_timepoint_info)){
    patient=as.character(sample_timepoint_info[i,"Patient_ID"])
    cat("\n\nstart", patient,"clonal evolution analysis...\n\n")
    patient_output_dir<-file.path(output_dir,patient)
    if(!dir.exists(patient_output_dir)){dir.create(patient_output_dir,recursive = T)}
    samples<-as.character(sample_timepoint_info[i,-match("Patient_ID",colnames(sample_timepoint_info))])
    names(samples)<-colnames(sample_timepoint_info)[-match("Patient_ID",colnames(sample_timepoint_info))]
    samples<-samples[!is.na(samples)]

    if(!is.null(sequenza_snps)){
      patient_sequenza_snps<-sequenza_snps[sequenza_snps$Sample_ID %in% samples,]
      patient_sequenza_snps<-patient_sequenza_snps %>% dplyr::group_by(Mutation_ID) %>% dplyr::mutate(n_sample_inpatient=dplyr::n()) %>% dplyr::filter(n_sample_inpatient>1)
      patient_sequenza_snps<-patient_sequenza_snps[!is.na(patient_sequenza_snps$gene) & patient_sequenza_snps$n_sample_inpatient>1,]
      patient_sequenza_snps<-patient_sequenza_snps %>% dplyr::group_by(Mutation_ID) %>% dplyr::mutate(patient_refalt_counts=sum(ref_counts)+sum(alt_counts))
      patient_sequenza_snps<-patient_sequenza_snps %>% dplyr::group_by(gene) %>% dplyr::filter(patient_refalt_counts==max(patient_refalt_counts))
      patient_sequenza_snps<-patient_sequenza_snps[,c("Sample_ID","Mutation_ID","chrom","start","end","gene","ref_counts","alt_counts")]
    }

    patient_mutect_snps<-mutect_snps[mutect_snps$Sample_ID %in% samples,c("Sample_ID","chrom","start","t_ref_count","t_alt_count","gene")]
    patient_mutect_snps<-with(patient_mutect_snps,data.frame(Sample_ID=Sample_ID,
                                                             Mutation_ID=paste(paste("chr",chrom,sep=""),start,sep=":"),
                                                             chrom=paste("chr",chrom,sep=""),
                                                             start=start,
                                                             end=start,
                                                             gene=gene,
                                                             ref_counts=t_ref_count,
                                                             alt_counts=t_alt_count))
    ol_patient_mutect_snps<- patient_mutect_snps %>% dplyr::group_by(Mutation_ID) %>% dplyr::mutate(n_sample_inpatient=dplyr::n()) %>% dplyr::filter(n_sample_inpatient>1)
    nol_patient_mutect_snps<- patient_mutect_snps %>% dplyr::group_by(Mutation_ID) %>% dplyr::mutate(n_sample_inpatient=dplyr::n()) %>% dplyr::filter(n_sample_inpatient==1 & ref_counts>=nonoverlapping_ref_counts_cutoff & alt_counts>=nonoverlapping_alt_counts_cutoff)
    patient_mutect_snps<-rbind(ol_patient_mutect_snps,nol_patient_mutect_snps)
    patient_mutect_snps<-patient_mutect_snps[,c("Sample_ID","Mutation_ID","chrom","start","end","gene","ref_counts","alt_counts")]

    ref_counts_<-tidyr::pivot_wider(patient_mutect_snps,id_cols="Mutation_ID",names_from = "Sample_ID",values_from = "ref_counts",values_fn = function(i){mean(i,na.rm=T)})
    ref_counts_[,-1]<-ref_counts_[,-1] %>% dplyr::mutate_if(is.numeric,~ifelse(is.na(.x), as.integer(mean(.x, na.rm = TRUE)), .x))
    ref_counts<-tidyr::pivot_longer(ref_counts_,cols = unname(samples),names_to = "Sample_ID",values_to = "ref_counts")
    alt_counts_<-tidyr::pivot_wider(patient_mutect_snps,id_cols="Mutation_ID",names_from = "Sample_ID",values_from = "alt_counts",values_fill = 0)
    alt_counts<-tidyr::pivot_longer(alt_counts_,cols = unname(samples),names_to = "Sample_ID",values_to = "alt_counts")
    ref_alt_counts<-merge(ref_counts,alt_counts,by=c("Mutation_ID","Sample_ID"))
    patient_mutect_snps_<-data.frame(Sample_ID=ref_alt_counts$Sample_ID,
                                     Mutation_ID=ref_alt_counts$Mutation_ID,
                                     chrom=gsub(":.*","",ref_alt_counts$Mutation_ID),
                                     start=as.integer(gsub(".*:","",ref_alt_counts$Mutation_ID)),
                                     end=as.integer(gsub(".*:","",ref_alt_counts$Mutation_ID)),
                                     ref_counts=ref_alt_counts$ref_counts,
                                     alt_counts=ref_alt_counts$alt_counts
    )

    patient_mutect_snps<-merge(patient_mutect_snps_,patient_mutect_snps[,c("Mutation_ID","gene")][!duplicated(patient_mutect_snps[,c("Mutation_ID")]),],by=c("Mutation_ID"),all.x=T)
    patient_mutect_snps<-patient_mutect_snps[,c("Sample_ID","Mutation_ID","chrom","start","end","gene","ref_counts","alt_counts")]
    if(!is.null(sequenza_snps)){
      patient_snps<-rbind(patient_sequenza_snps[!patient_sequenza_snps$Mutation_ID %in% patient_mutect_snps$Mutation_ID,],patient_mutect_snps)
    } else {patient_snps<-patient_mutect_snps}
    patient_snps[["is_cancer_gene"]]<-ifelse(patient_snps[["gene"]] %in% cancer_genes,TRUE,FALSE)
    patient_snps[["is_smg_gene"]]<-ifelse(patient_snps[["gene"]] %in% smg_genes,TRUE,FALSE)

    patient_snps_<-patient_snps[,c("Sample_ID","chrom","start","end","ref_counts","alt_counts")]
    patient_snps_$chrom<-gsub("chr","",patient_snps_$chrom)
    patient_segments<-sequenza_segments[sequenza_segments$Sample_ID %in% samples,]
    patient_purities<-purities[purities$Sample_ID %in% samples,]

    patient_snp_segments<-utiltools::prepare_snpsegments_for_pyclone_VI(snps = patient_snps_,segments = patient_segments,purities = patient_purities)
    patient_snp_segments<-patient_snp_segments[patient_snp_segments$mutation_id %in% patient_snp_segments$mutation_id[duplicated(patient_snp_segments$mutation_id)],]
    patient_snp_segments$sample_id<-names(samples)[match(patient_snp_segments$sample_id,samples)]
    patient_snp_segments<-patient_snp_segments[complete.cases(patient_snp_segments),]
    patient_snp_segments_<-patient_snp_segments
    #patient_snp_segments$mutation_id<-gsub("chr","",patient_snp_segments$mutation_id)
    patient_snp_segments_[['vaf']]<-with(patient_snp_segments_,round(alt_counts/(ref_counts*(tumour_content*major_cn)/((1-tumour_content)*normal_cn+(tumour_content*major_cn))+alt_counts)*100,2))
    #mutation_vaf<-pivot_wider(patient_snp_segments_,id_cols = "mutation_id",names_from = "sample_id",values_from = "vaf",values_fill = 0)
    #vaf_plot<-ggplot(mutation_vaf,aes(x=Baseline,y=TP2))+geom_point()

    #browser()
    patient_snp_segments_file<-file.path(patient_output_dir,"snp_segemnt.tsv")
    write.table(patient_snp_segments,patient_snp_segments_file,row.names = F,quote = F,sep="\t")
    tracerx_file<-file.path(patient_output_dir,"tracerx.h5")
    fit_command<-paste(pyclone_path,"fit -i",patient_snp_segments_file,"-o",tracerx_file,"-c 10 -d beta-binomial -r 10",sep=" ")
    result_file<-file.path(patient_output_dir,"tracerx.tsv")
    write_command<-paste(pyclone_path,"write-results-file -i",tracerx_file, "-o", result_file,sep=" ")
    system(fit_command,intern = F)
    system(write_command,intern = F)
    prevalence<-read.csv(result_file,header=T,stringsAsFactors = F,check.names = F,sep='\t')
    prevalence$cluster<-prevalence$cluster_id+1L
    prevalence$cellular_prevalence<-prevalence$cellular_prevalence*100
    founding.cluster<-prevalence %>% dplyr::filter(sample_id=="Baseline") %>% dplyr::group_by(cluster) %>% dplyr::summarise(mean_ccf=mean(cellular_prevalence,na.rm=T))
    founding.cluster<-founding.cluster$cluster[which.max(founding.cluster$mean_ccf)]
    prevalence<-merge(prevalence,patient_snp_segments_,by=c("mutation_id","sample_id"),all.x=T)
    prevalence<-utiltools::addccf(prevalence)
    prevalence$ccf<-prevalence$ccf*100
    mutation_anno<-patient_snps[,c("Mutation_ID","gene","is_cancer_gene","is_smg_gene")]
    mutation_anno<-mutation_anno[!duplicated(mutation_anno),]
    prevalence<-merge(prevalence,mutation_anno,by.x="mutation_id",by.y="Mutation_ID",all.x=T)

    mutation_vaf<-tidyr::pivot_wider(prevalence,id_cols = "mutation_id",names_from = "sample_id",values_from = "vaf",values_fill = 0)
    mutation_vaf<-merge(mutation_vaf,prevalence[,c("mutation_id","cluster")][!duplicated(prevalence[,c("mutation_id","cluster")]),],by="mutation_id")
    vaf_plot<-ggplot2::ggplot(mutation_vaf,ggplot2::aes(x=Baseline,y=TP2,color=factor(cluster)))+ggplot2::geom_point()

    ccf<-tidyr::pivot_wider(prevalence,id_cols=c('mutation_id',"cluster"),names_from = c("sample_id"),values_from = c("vaf"))
    ccf<-merge(ccf,mutation_anno,by.x="mutation_id",by.y="Mutation_ID",all.x=T)
    ccf$driver_gene_rank<-ccf$is_cancer_gene*2+ccf$is_smg_gene
    ccf<-ccf[order(ccf$driver_gene_rank,decreasing = T),]
    ccf<-ccf %>% dplyr::group_by(cluster) %>% dplyr::mutate(fake_rank=dplyr::row_number())
    ccf[["is_driver_gene"]]<-ifelse(ccf$is_cancer_gene & ccf$fake_rank<=3,TRUE,FALSE)
    ccf[["cluster"]]<-as.character(ccf[["cluster"]])
    ccf<-as.data.frame(ccf)
    #browser()
    vaf.col.names<-names(samples)
    sample.groups<-names(samples)
    names(sample.groups)<-names(samples)
    clone.colors<-RColorBrewer::brewer.pal(length(unique(ccf[["cluster"]])),"Paired")[1:length(unique(ccf[["cluster"]]))]
    svg(filename = file.path(patient_output_dir,"variantcluster_plot.svg"),width = 8,height=8)
    variantcluster_plot <- clonevol::plot.variant.clusters(ccf,
                                                           cluster.col.name = 'cluster',
                                                           show.cluster.size = FALSE,
                                                           cluster.size.text.color = 'blue',
                                                           vaf.col.names = vaf.col.names,
                                                           vaf.limits = 70,
                                                           sample.title.size = 20,
                                                           violin = FALSE,
                                                           box = FALSE,
                                                           jitter = TRUE,
                                                           jitter.shape = 1,
                                                           jitter.color = clone.colors,
                                                           jitter.size = 3,
                                                           jitter.alpha = 1,
                                                           jitter.center.method = 'median',
                                                           jitter.center.size = 1,
                                                           jitter.center.color = 'darkgray',
                                                           jitter.center.display.value = 'none',
                                                           highlight = 'is_driver_gene',
                                                           highlight.shape = 21,
                                                           highlight.color = 'blue',
                                                           highlight.fill.color = 'green',
                                                           highlight.note.col.name = 'gene',
                                                           highlight.note.size = 2,
                                                           order.by.total.vaf = FALSE)
    dev.off()

    clusterflow_plot <-clonevol::plot.cluster.flow(ccf, vaf.col.names = vaf.col.names,
                                                   sample.names = sample.groups,
                                                   colors = clone.colors)
    if(save_analysis){
      svg(filename = file.path(patient_output_dir,"vaf_plot.svg"),width = 8,height=8)
      print(vaf_plot)
      dev.off()
      svg(filename = file.path(patient_output_dir,"clusterflow_plot.svg"),width = 8,height=8)
      print(clusterflow_plot)
      dev.off()
      write.csv(patient_snp_segments,file.path(patient_output_dir,"patient_snp_segments.csv"),row.names = F,quote=F)
      write.csv(prevalence,file.path(patient_output_dir,"prevalence.csv"),row.names = F,quote=F)
      write.csv(ccf,file.path(patient_output_dir,"ccf.csv"),row.names = F,quote=F)
      cat(obj@log,file=file.path(output_dir,"log.txt"))
    }
    obj@output[[patient]][["patient_snp_segments"]]<-patient_snp_segments
    obj@output[[patient]][["prevalence"]]<-prevalence
    obj@output[[patient]][["ccf"]]<-ccf
    model = clonevol::infer.clonal.models(variants = ccf,
                                          cluster.col.name = "cluster",
                                          vaf.col.names=vaf.col.names,
                                          sample.groups = sample.groups,
                                          cancer.initiation.model='monoclonal',
                                          subclonal.test = 'bootstrap',
                                          subclonal.test.model = 'non-parametric',
                                          num.boots = 1000,
                                          founding.cluster = founding.cluster,
                                          cluster.center = 'median',
                                          ignore.clusters = NULL,
                                          clone.colors = clone.colors,
                                          min.cluster.vaf = 0.01,
                                          # min probability that CCF(clone) is non-negative
                                          sum.p = 0.05,
                                          # alpha level in confidence interval estimate for CCF(clone)
                                          alpha = 0.05)

    if(model$num.matched.models==0){
      cat("No clonal model found for",patient,"!")
      ignore.clusters_<-ccf %>% dplyr::group_by(cluster) %>% dplyr::summarize(n_gene=dplyr::n())
      ignore.clusters<-ignore.clusters_$cluster[which.min(ignore.clusters_$n_gene)]
      model = clonevol::infer.clonal.models(variants = ccf,
                                            cluster.col.name = "cluster",
                                            vaf.col.names=vaf.col.names,
                                            sample.groups = sample.groups,
                                            cancer.initiation.model='monoclonal',
                                            subclonal.test = 'bootstrap',
                                            subclonal.test.model = 'non-parametric',
                                            num.boots = 1000,
                                            founding.cluster = founding.cluster,
                                            cluster.center = 'median',
                                            ignore.clusters = ignore.clusters,
                                            clone.colors = clone.colors,
                                            min.cluster.vaf = 0.01,
                                            # min probability that CCF(clone) is non-negative
                                            sum.p = 0.05,
                                            # alpha level in confidence interval estimate for CCF(clone)
                                            alpha = 0.05)

    }
    model<-clonevol::transfer.events.to.consensus.trees(model,ccf[ccf$is_driver_gene,],cluster.col.name = "cluster",event.col.name = "gene")
    model <- clonevol::convert.consensus.tree.clone.to.branch(model, branch.scale = 'sqrt')
    for(i in seq_len(length(model$matched$merged.trees))){ # guarantee Y chunk appear first
      branches<-model$matched$merged.trees[[i]]$branches
      branch_ord<-c(which(branches=="Y"),setdiff(order(branches),which(branches=="Y")))
      model$matched$merged.trees[[i]]<-model$matched$merged.trees[[i]][branch_ord,]
    }
    #browser()
    obj@output[[patient]][["model"]]<-model

    modelclone_plot<-clonevol::plot.clonal.models(model,
                                                  # box plot parameters
                                                  box.plot = TRUE,
                                                  fancy.boxplot = TRUE,
                                                  fancy.variant.boxplot.highlight = 'is_driver_gene',
                                                  fancy.variant.boxplot.highlight.shape = 21,
                                                  fancy.variant.boxplot.highlight.fill.color = 'red',
                                                  fancy.variant.boxplot.highlight.color = 'black',
                                                  fancy.variant.boxplot.highlight.note.col.name = 'gene',
                                                  fancy.variant.boxplot.highlight.note.color = 'blue',
                                                  fancy.variant.boxplot.highlight.note.size = 2,
                                                  fancy.variant.boxplot.jitter.alpha = 1,
                                                  fancy.variant.boxplot.jitter.center.color = 'grey50',
                                                  fancy.variant.boxplot.base_size = 12,
                                                  fancy.variant.boxplot.plot.margin = 1,
                                                  fancy.variant.boxplot.vaf.suffix = '.VAF',
                                                  # bell plot parameters
                                                  clone.shape = 'bell',
                                                  bell.event = TRUE,
                                                  bell.event.label.color = 'blue',
                                                  bell.event.label.angle = 60,
                                                  clone.time.step.scale = 1,
                                                  bell.curve.step = 2,
                                                  # node-based consensus tree parameters
                                                  merged.tree.plot = TRUE,
                                                  tree.node.label.split.character = NULL,
                                                  tree.node.shape = 'circle',
                                                  tree.node.size = 30,
                                                  tree.node.text.size = 0.5,
                                                  merged.tree.node.size.scale = 1.25,
                                                  merged.tree.node.text.size.scale = 2.5,
                                                  merged.tree.cell.frac.ci = FALSE,
                                                  # branch-based consensus tree parameters
                                                  merged.tree.clone.as.branch = TRUE,
                                                  mtcab.event.sep.char = ',',
                                                  mtcab.branch.text.size = 1,
                                                  mtcab.branch.width = 0.75,
                                                  mtcab.node.size = 3,
                                                  mtcab.node.label.size = 1,
                                                  mtcab.node.text.size = 1.5,
                                                  # cellular population parameters
                                                  cell.plot = TRUE,
                                                  num.cells = 100,
                                                  cell.border.size = 0.25,
                                                  cell.border.color = 'black',
                                                  clone.grouping = 'horizontal',
                                                  #meta-parameters
                                                  scale.monoclonal.cell.frac = TRUE,
                                                  show.score = FALSE,
                                                  cell.frac.ci = TRUE,
                                                  disable.cell.frac = FALSE,
                                                  # output figure parameters
                                                  out.dir = patient_output_dir,
                                                  out.format = 'pdf',
                                                  overwrite.output = TRUE,
                                                  width = 12,
                                                  height = 4,
                                                  # vector of width scales for each panel from left to right
                                                  panel.widths = c(3,4,2,4,2))
  }
  obj@output[["sample_timepoint_info"]]<-sample_timepoint_info
  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")
  return(obj)
}
)


#' runICGDifAnalysis
#' run immune checkpoint gene analysis
#' @param obj analysis. an analysis object.
#' @param icg_info data.frame. immune checkpoint gene information
#' @param gene_id_column character. gene ic column from icg_info
#' @score_method character. method passed to utiltools::geneset_activity, used for calculation of activity
#' @param survival_column character. suvival days column from patient_info
#' @param status_column character. suvival status column from patient_info
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param palette character. colors used for heatmap.
#' @param group_column character. group column from patient_info/sample_info
#' @param block_column character. block column from patient_info/sample_info
#' @param contrasts character. contrast for test.
#' @param pval_cutoff numeric. p value cutoff for significance.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory.
#' @param log chracter. any comments.
#'
#' @return analysis. contain various contrast result.
#' @export
#'
setGeneric("runICGDifAnalysis",function(obj,
                                        icg_info,gene_id_column="Hgnc_Symbol",
                                        score_method="zscore",
                                        survival_column,status_column,
                                        top_anno,top_anno_sample_order,palette,
                                        group_column,block_column,contrasts,
                                        pval_cutoff = 0.05,
                                        save_analysis = T,output_dir,
                                        log="run immune checkpoint differential gene expression analysis.") standardGeneric("runICGDifAnalysis"))
setMethod("runICGDifAnalysis","Analysis",function(obj,
                                                  icg_info,gene_id_column="Hgnc_Symbol",
                                                  score_method="zscore",
                                                  survival_column,status_column,
                                                  top_anno,top_anno_sample_order,palette,
                                                  group_column,block_column,contrasts,
                                                  pval_cutoff = 0.05,
                                                  save_analysis = T,output_dir,
                                                  log="run immune checkpoint differential gene expression analysis."){

  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  } else {
    output_dir<-getwd()
  }

  RNA_samples<-top_anno_sample_order
  sample_info<-merge(merge(data.frame(Sample_ID=RNA_samples),obj@experiment@sample_info,by="Sample_ID",all.x=T),obj@experiment@patient_info,by="Patient_ID",all.x=T)
  sample_info<-sample_info[match(RNA_samples, sample_info$Sample_ID),]

  icg<-as.data.frame(obj@experiment@mRNA_assay@assay_data[na.omit(intersect(rownames(obj@experiment@mRNA_assay@assay_data),icg_info[[gene_id_column]])),RNA_samples])
  icg<-icg[rowMeans(icg)>=2,]

  if(length(group_column)==1){
    profile_group_<-sample_info[[group_column]]
    if(any(is.na(profile_group_) | profile_group_=="")){
      excluded_samples<-sample_info[["Sample_ID"]][is.na(profile_group_) | profile_group_==""]
      sample_info=sample_info[-match(excluded_samples,sample_info[["Sample_ID"]]),]
      RNA_samples=RNA_samples[-match(excluded_samples,RNA_samples)]
      top_anno<-top_anno[match(RNA_samples,top_anno_sample_order)]
      assertthat::are_equal(sample_info[["Sample_ID"]],RNA_samples)
      icg=icg[,RNA_samples]
    }
  }

  icg_info_<-icg_info[match(rownames(icg),icg_info[[gene_id_column]]),]
  icg_function_anno<-ComplexHeatmap::rowAnnotation(Function=icg_info_[["Immune_Checkpoint"]],col=list(Function=c("Inhibitory"="red","Stimulatory"="blue","MHC"="green")))

  signatures<-list(Inhibitory=icg_info_$Hgnc_Symbol[icg_info_$Immune_Checkpoint=="Inhibitory"],
                   Stimulatory=icg_info_$Hgnc_Symbol[icg_info_$Immune_Checkpoint=="Stimulatory"],
                   MHC=icg_info_$Hgnc_Symbol[icg_info_$Immune_Checkpoint=="MHC"])

  scores<-utiltools::geneset_activity(icg,genesets = signatures,methods = score_method)
  ic_scores<-scores[["agg_activity"]]
  if(length(unique(sample_info[[group_column]]))==2){
    ic_score_statistics<-apply(rbind(ic_scores,ic_scores["Inhibitory",]-ic_scores["Stimulatory",]),1,function(r_data){data=data.frame(Value=as.numeric(r_data),Group=sample_info[[group_column]]);test=t.test(Value~Group,data=data);return(test$p.value)})  #ic_scores<-as.data.frame(t(ic_scores))
  } else {
    ic_score_statistics<-apply(rbind(ic_scores,ic_scores["Inhibitory",]-ic_scores["Stimulatory",]),1,function(r_data){data=data.frame(Value=as.numeric(r_data),Group=sample_info[[group_column]]);test=stats::aov(Value~Group,data=data);return(summary(test)[[1]][["Pr(>F)"]][1])})  #ic_scores<-as.data.frame(t(ic_scores))
  }
    assertthat::are_equal(colnames(ic_scores),colnames(icg))
  #browser()
  group_mean_<-list()
  for (icg_group in unique(sample_info[[group_column]])){
    group_mean_[[paste(icg_group,"_mean",sep="")]]<-rowMeans(icg[,sample_info[[group_column]]==icg_group,drop=F],na.rm=T)
  }
  group_results<-as.data.frame(group_mean_)
  panel_fun<-function(index,nm){
    df<-data.frame(Value=as.numeric(ic_scores[nm,]),Group=factor(sample_info[[group_column]]))
    colors<-c("blue","red","green")[1:length(levels(df$Group))]
    names(colors)<-levels(df$Group)
    xscale=c(0,max(as.numeric(df$Group))+1)
    yscale=range(df$Value)+c(-0.3,0.3)*abs(range(df$Value))
    pushViewport(viewport(xscale=xscale,yscale=yscale))
    grid.rect()
    grid.xaxis(gp=gpar(fontsize=6),at=1:max(as.numeric(df$Group)),label = levels(df$Group))
    grid.yaxis(gp=gpar(fontsize=6),main=T)
    for(g in levels(df$Group)){
      grid.boxplot(df$Value[df$Group==g],pos=which(levels(df$Group)==g),direction = "vertical",gp=gpar(fill=colors[g]))
    }
    grid.text(paste("pval=",round(ic_score_statistics[nm],3),sep=" "),x = 0.3,y=0.9,gp=gpar(fontsize=8))
    popViewport()
  }
  if(missing(block_column)){
    if(missing(contrasts)){
      if(length(unique(sample_info[[group_column]]))>2){
        group_statistics=apply(icg,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[[group_column]]);test=aov(value~group,data=df);return(summary(test)[[1]][["Pr(>F)"]][1])})
      } else {
        group_statistics=apply(icg,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[[group_column]]);test=t.test(value~group,data=df);return(test$p.value)})
      }
      group_results<-data.frame(group_results,P.Value=group_statistics)
      neg_log_pvalue=-log10(group_statistics)
      icg_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(icg))),name="zscores_icg",top_annotation = top_anno,row_split=icg_info_$Immune_Checkpoint,column_split = sample_info[[group_column]])
      row_anno<-rowAnnotation(Function=anno_block(gp=gpar(fill=c("red","blue","green")),labels=names(ComplexHeatmap::row_order(icg_heatmap))),
                              `-log10(p)`=ComplexHeatmap::row_anno_barplot(neg_log_pvalue,gp=gpar(fill=ifelse(neg_log_pvalue>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)),
                              Score=anno_link(align_to=icg_info_$Immune_Checkpoint,which="row",panel_fun=panel_fun,size=unit(3,"cm"),width=unit(4,"cm")),
                              gap=unit(2,"points"),annotation_name_side="top",annotation_name_rot=0)
      fontsize=rep(6,nrow(icg))
      fontsize[neg_log_pvalue>=(-log10(pval_cutoff))]<-8
      fontface=rep("plain",nrow(icg))
      fontface[neg_log_pvalue>=(-log10(pval_cutoff))]<-"bold"
      icg_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(icg))),name="zscores_icg",top_annotation = top_anno,row_split=icg_info_$Immune_Checkpoint,right_annotation = row_anno,column_split = sample_info[[group_column]],row_names_gp = gpar(fontsize=fontsize,fontface=fontface),row_title_gp = gpar(fontsize=0),row_names_side = "left")
      obj@output$contrasts$icg_heatmap<-icg_heatmap
      obj@output$contrasts$statistics<-group_results
      if(save_analysis){
        if(!dir.exists(file.path(output_dir,"contrasts"))){dir.create(file.path(output_dir,"contrasts"),recursive = T)}
        write.csv(group_results,file.path(output_dir,"contrasts","statistics.csv"))
        svg(filename = file.path(output_dir,"contrasts","icg_heatmap.svg"),width = 20,height=15)
        ComplexHeatmap::draw(icg_heatmap,merge_legend=T)
        if(any(neg_log_pvalue>(-(log10(pval_cutoff))))){
          xpos<-0.95/(max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
          for(i in 1:length(unique(icg_info_$Immune_Checkpoint))){
            decorate_annotation("-log10(p)",{
              grid.lines(c(xpos,xpos),c(0,1),gp = gpar(lty = 2,col="red"))
            },slice=i)
          }

        }
        dev.off()
      }
    } else {
      for(contrast in contrasts){
        subset_sample_info<-sample_info[sample_info[[group_column]] %in% unlist(strsplit(contrast,"-")),]
        subset_icg<-icg[,subset_sample_info$Sample_ID]
        contrast_group_statistics=apply(subset_icg,1,function(scf){df<-data.frame(value=as.numeric(scf),group=subset_sample_info[[group_column]]);test=t.test(value~group,data=df);return(test$p.value)})
        contrast_group_results[[paste(contrast,"P.Value",sep=":")]]<-contrast_group_statistics
        contrast_top_anno<-top_anno[match(subset_sample_info$Sample_ID,RNA_samples),]
        neg_log_pvalue=-log10(group_statistics)
        icg_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(icg))),name="zscores_icg",top_annotation = top_anno,row_split=icg_info_$Immune_Checkpoint,column_split = sample_info[[group_column]])
        row_anno<-rowAnnotation(Function=anno_block(gp=gpar(fill=c("red","blue","green")),labels=names(ComplexHeatmap::row_order(icg_heatmap))),
                                `-log10(p)`=ComplexHeatmap::row_anno_barplot(neg_log_pvalue,gp=gpar(fill=ifelse(neg_log_pvalue>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)),
                                Score=anno_link(align_to=icg_info_$Immune_Checkpoint,which="row",panel_fun=panel_fun,size=unit(3,"cm"),width=unit(4,"cm")),
                                gap=unit(2,"points"),annotation_name_side="top",annotation_name_rot=0)
        fontsize=rep(6,nrow(icg))
        fontsize[neg_log_pvalue>=(-log10(pval_cutoff))]<-8
        fontface=rep("plain",nrow(icg))
        fontface[neg_log_pvalue>=(-log10(pval_cutoff))]<-"bold"
        contrast_icg_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(subset_icg))),name="zscores_icg",top_annotation = contrast_top_anno,column_split = subset_sample_info[[group_column]],row_split=icg_info_$Immune_Checkpoint,right_annotation = row_anno,row_names_gp = gpar(fontsize=fontsize,fontface=fontface))

        obj@output$contrasts[[contrast]]$sample_info<-subset_sample_info
        obj@output$contrasts[[contrast]]$icg<-subset_icg
        obj@output$contrasts[[contrast]]$icg_heatmap<-contrast_icg_heatmap
        obj@output$contrasts[[contrast]]$statistics<-contrast_group_results
        if(save_analysis){
          if(!dir.exists(file.path(output_dir,"contrasts",contrast))){dir.create(file.path(output_dir,"contrasts",contrast),recursive = T)}
          write.csv(contrast_icg_heatmap,file.path(output_dir,"contrasts",contrast,"statistics.csv"))
          write.csv(subset_sample_info,file.path(output_dir,"contrasts",contrast,"sample_info.csv"))
          write.csv(subset_icg,file.path(output_dir,"contrasts",contrast,"subset_icg.csv"))
          svg(filename = file.path(output_dir,"contrasts",contrast,"contrast_icg_heatmap.svg"),width = 20,height=15)
          ComplexHeatmap::draw(contrast_icg_heatmap,merge_legend=T)
          if(any(neg_log_pvalue>(-(log10(pval_cutoff))))){
            xpos<-0.95/(max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
            for(i in 1:length(unique(icg_info_$`Immune Checkpoint`))){
              decorate_annotation("-log10(p)",{
                grid.lines(c(xpos,xpos),c(0,1),gp = gpar(lty = 2,col="red"))
              },slice=i)
            }
          }
          dev.off()
        }
      }
    }
  } else {
    sample_info<-sample_info[order(sample_info[[block_column]],sample_info[[group_column]]),]
    icg<-icg[,sample_info$Sample_ID]
    if(missing(contrasts)){
      if(length(unique(sample_info[[group_column]]))>2){
        group_statistics=apply(icg,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[[group_column]]);test=aov(value~group,data=df);return(summary(test)[[1]][["Pr(>F)"]][1])})
      } else {
        group_statistics=apply(icg,1,function(scf){df<-data.frame(value=as.numeric(scf),group=sample_info[[group_column]]);test=t.test(value~group,data=df,paired=T);return(test$p.value)})
      }
      group_results<-data.frame(group_results,P.Value=group_statistics)
      neg_log_pvalue=-log10(group_statistics)
      icg_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(icg))),name="zscores_icg",top_annotation = top_anno,row_split=icg_info_$Immune_Checkpoint,column_split = sample_info[[group_column]])
      row_anno<-rowAnnotation(Function=anno_block(gp=gpar(fill=c("red","blue","green")),labels=names(ComplexHeatmap::row_order(icg_heatmap))),
                              `-log10(p)`=ComplexHeatmap::row_anno_barplot(neg_log_pvalue,gp=gpar(fill=ifelse(neg_log_pvalue>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)),
                              Score=anno_link(align_to=icg_info_$Immune_Checkpoint,which="row",panel_fun=panel_fun,size=unit(3,"cm"),width=unit(4,"cm")),
                              gap=unit(2,"points"),annotation_name_side="top",annotation_name_rot=0)
      fontsize=rep(6,nrow(icg))
      fontsize[neg_log_pvalue>=(-log10(pval_cutoff))]<-8
      fontface=rep("plain",nrow(icg))
      fontface[neg_log_pvalue>=(-log10(pval_cutoff))]<-"bold"
      icg_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(icg))),name="zscores_icg",top_annotation = top_anno,column_split = sample_info[[group_column]],row_split=icg_info_$Immune_Checkpoint,right_annotation = row_anno,row_names_gp = gpar(fontsize=fontsize,fontface=fontface))

      obj@output$contrasts$icg_heatmap<-icg_heatmap
      obj@output$contrasts$statistics<-group_results
      if(save_analysis){
        if(!dir.exists(file.path(output_dir,"contrasts"))){dir.create(file.path(output_dir,"contrasts"),recursive = T)}
        write.csv(group_results,file.path(output_dir,"contrasts","statistics.csv"))
        svg(filename = file.path(output_dir,"contrasts","icg_heatmap.svg"),width = 20,height=15)
        ComplexHeatmap::draw(icg_heatmap)
        if(any(neg_log_pvalue>(-(log10(pval_cutoff))))){
          xpos<-0.95/(max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
          for(i in 1:length(unique(icg_info_$`Immune Checkpoint`))){
            decorate_annotation("-log10(p)",{
              grid.lines(c(xpos,xpos),c(0,1),gp = gpar(lty = 2,col="red"))
            },slice=i)
          }
        }
        dev.off()
      }
    } else {
      for(contrast in contrasts){
        subset_sample_info<-sample_info[sample_info[[group_column]] %in% unlist(strsplit(contrast,"-")),]
        subset_icg<-icg[,subset_sample_info$Sample_ID]
        contrast_group_statistics=apply(subset_icg,1,function(scf){df<-data.frame(value=as.numeric(scf),group=subset_sample_info[[group_column]]);test=t.test(value~group,data=df,paired=T);return(test$p.value)})
        contrast_group_results[[paste(contrast,"P.Value",sep=":")]]<-contrast_group_statistics
        contrast_top_anno<-top_anno[match(subset_sample_info$Sample_ID,RNA_samples),]
        neg_log_pvalue=-log10(group_statistics)
        icg_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(icg))),name="zscores_icg",top_annotation = top_anno,row_split=icg_info_$Immune_Checkpoint,column_split = sample_info[[group_column]])
        row_anno<-rowAnnotation(Function=anno_block(gp=gpar(fill=c("red","blue","green")),labels=names(ComplexHeatmap::row_order(icg_heatmap))),
                                `-log10(p)`=ComplexHeatmap::row_anno_barplot(neg_log_pvalue,gp=gpar(fill=ifelse(neg_log_pvalue>(-log10(pval_cutoff)),"red","gray")),ylim=c(0,max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)),
                                Score=anno_link(align_to=icg_info_$Immune_Checkpoint,which="row",panel_fun=panel_fun,size=unit(3,"cm"),width=unit(4,"cm")),
                                gap=unit(2,"points"),annotation_name_side="top",annotation_name_rot=0)
        fontsize=rep(6,nrow(icg))
        fontsize[neg_log_pvalue>=(-log10(pval_cutoff))]<-8
        fontface=rep("plain",nrow(icg))
        fontface[neg_log_pvalue>=(-log10(pval_cutoff))]<-"bold"
        contrast_icg_heatmap<-ComplexHeatmap::Heatmap(t(scale(t(subset_icg))),name="zscores_icg",top_annotation = contrast_top_anno,column_split = subset_sample_info[[group_column]],row_split=icg_info_$Immune_Checkpoint,right_annotation = row_anno,row_names_gp = gpar(fontsize=fontsize,fontface=fontface))

        obj@output$contrasts[[contrast]]$sample_info<-subset_sample_info
        obj@output$contrasts[[contrast]]$icg<-subset_icg
        obj@output$contrasts[[contrast]]$icg_heatmap<-contrast_icg_heatmap
        obj@output$contrasts[[contrast]]$statistics<-contrast_group_results
        if(save_analysis){
          if(!dir.exists(file.path(output_dir,"contrasts",contrast))){dir.create(file.path(output_dir,"contrasts",contrast),recursive = T)}
          write.csv(contrast_icg_heatmap,file.path(output_dir,"contrasts",contrast,"statistics.csv"))
          write.csv(subset_sample_info,file.path(output_dir,"contrasts",contrast,"sample_info.csv"))
          write.csv(subset_icg,file.path(output_dir,"contrasts",contrast,"subset_icg.csv"))
          svg(filename = file.path(output_dir,"contrasts",contrast,"contrast_icg_heatmap.svg"),width = 20,height=15)
          ComplexHeatmap::draw(contrast_icg_heatmap)
          if(any(neg_log_pvalue>(-(log10(pval_cutoff))))){
            xpos<-0.95/(max(c(neg_log_pvalue,-log10(pval_cutoff)))+0.5)*(-log10(pval_cutoff))
            for(i in 1:length(unique(icg_info_$`Immune Checkpoint`))){
              decorate_annotation("-log10(p)",{
                grid.lines(c(xpos,xpos),c(0,1),gp = gpar(lty = 2,col="red"))
              },slice=i)
            }
          }
          dev.off()
        }
      }
    }
  }
  obj@output$sample_info<-sample_info
  obj@output$icg<-icg

  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")

  return(obj)
}
)


#' runAssaySurvivalSignatureAnalysis
#'
#' run assay survival signature analysis
#'
#' @param obj Analysis. analysis object to be run.
#' @param assay_name character. assay name.
#' @param assay_type character. assay type, such as "DNA" or "RNA"
#' @param scale logic. whether to scale assay.
#' @param patient_id_column character. column for patient id.
#' @param sample_id_column character. column for sample id.
#' @param survival_column character. column for survival days.
#' @param status_column character. column for survival status.
#' @param suvival_type character. type of survival analysis, such as OS or PFS.
#' @param top_anno HeatmapAnnotation. top annotation, passed to heatmap.
#' @param top_anno_sample_order character. sample orders for top annotation.
#' @param palette character. colors for KM curve.
#' @param uni_cox_pvalue_cutoff numeric. p value cutoff for significance of unicox test.
#' @param beta_filter numeric. p value cutoff for beta coefficient of unicox test.
#' @param save_analysis logic. whether to save results locally.
#' @param output_dir character. output directory, if missing, use project/analysis_name as directory
#' @param log character. any comments.
#'
#' @return Analysis.
#' @export
#'
setGeneric("runAssaySurvivalSignatureAnalysis",function(obj,
                                                        assay_name,assay_type="RNA",scale,
                                                        patient_id_column="Patient_ID",sample_id_column="Sample_ID",survival_column,status_column,survival_type=c("OS"),
                                                        top_anno,top_anno_sample_order,palette=NULL,
                                                        uni_cox_pvalue_cutoff = 0.01,beta_filter = 0.001,
                                                        save_analysis = T,output_dir,
                                                        log="run assay survival signature analysis.") standardGeneric("runAssaySurvivalSignatureAnalysis"))
setMethod("runAssaySurvivalSignatureAnalysis","Analysis",function(obj,
                                                                  assay_name,assay_type="RNA",scale,
                                                                  patient_id_column="Patient_ID",sample_id_column="Sample_ID",survival_column,status_column,survival_type=c("OS"),
                                                                  top_anno,top_anno_sample_order,palette=NULL,
                                                                  uni_cox_pvalue_cutoff = 0.01,beta_filter = 0.001,
                                                                  save_analysis = T,output_dir,
                                                                  log="run assay survival signature analysis."){

  if(save_analysis){
    if(missing(output_dir)){
      output_dir<-file.path(obj@project,obj@analysis_name)
      if(!dir.exists(output_dir)){dir.create(output_dir,recursive = T)}
    }
  } else {
    output_dir<-getwd()
  }
  if(missing(scale)){scale=F}

  sample_info<-obj@experiment@sample_info
  sample_info<-sample_info[sample_info$Assay==assay_type,]
  sample_info<-merge(sample_info,obj@experiment@patient_info,by=patient_id_column,all.x=T)
  excluded_samples<-sample_info[[sample_id_column]][is.na(sample_info[[survival_column]]) | is.na(sample_info[[status_column]])]

  assay_samples<-setdiff(top_anno_sample_order,excluded_samples)

  sample_info<-sample_info[match(assay_samples,sample_info[[sample_id_column]]),]
  assay_data=methods::slot(obj@experiment,assay_name)@assay_data[,assay_samples]
  top_anno=top_anno[match(assay_samples,top_anno_sample_order)]
  assertthat::are_equal(sample_info[["Sample_ID"]],assay_samples)

  if(missing(palette)){palette=c("blue","white","red")}

  signatures<-survivalSignatures(sample_info[[survival_column]],sample_info[[status_column]],assay_data,uni_cox_pvalue_cutoff = 0.01,beta_filter = 0.001,scale=scale)
  signature_expressions<-signatures$selected_covariates[,assay_samples]
  riskscore<-signatures$risk_scores[assay_samples]
  betas<-signatures$selected_betas[rownames(signature_expressions)]
  optimal_cutoff<-signatures$optimal_cut
  riskscore_anno=ComplexHeatmap::HeatmapAnnotation(RiskScore=anno_barplot(riskscore,gp=gpar(col=ifelse(riskscore>=optimal_cutoff,"green","red"))),annotation_name_side = "left")
  anno_signature<-c(top_anno,riskscore_anno)
  beta_anno<-ComplexHeatmap::rowAnnotation(Beta=anno_barplot(betas,gp=gpar(col=ifelse(betas>=0,"green","red"))))

  saveRDS(signatures,file=file.path(output_dir,"signature.rds"))
  write.csv(signature_expressions,file=file.path(output_dir,"signature.csv"))
  tiff(filename = file.path(output_dir,"signature_heatmap.tiff"),width = 20,height=20,units = "in",res=300,compression = "lzw")
  signature_heatmap<-Heatmap(t(scale(t(signature_expressions))),name="zscore_signature",top_annotation = anno_signature,right_annotation = beta_anno,column_order = names(sort(riskscore)),row_order = names(sort(betas)),row_names_gp = gpar(fontsize=5))
  print(signature_heatmap)
  dev.off()

  obj@output$signature<-signature
  obj@output$sample_info<-sample_info
  obj@output$signature_heatmap<-signature_heatmap

  obj@log=paste(obj@log,"\n\n",Sys.time(),";\n\t",log,sep="")

  return(obj)
}
)
