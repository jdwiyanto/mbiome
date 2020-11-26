#' Runs PLSDA and SPLSDA for the visualisation, feature-selection and validation of a microbiome dataset
#'
#' @param ps A phyloseq object -- raw read count
#' @param var A variable(s) for feature selection
#' @param valmodel A validation model. To be improved. Defaults to 'all'
#' @param no.comp Component numbers for PLSDA and SPLSDA to analyse
#' @param no.repeat No. of repeat to be performed by the model
#' @return output files in a new folder listing the compiled details and plots of SPLSDA and SPLSA along with the feature selection plot
#' @export
jd_splsda <- function(ps, var, valmodel = 'all', no.comp, no.repeat) {

  for (var in var) {

    sampledata = phyloseq::sample_data(ps)
    sampledata = subset(sampledata, !is.na(sampledata[[var]]))
    ps = phyloseq::phyloseq(sampledata, otu_table(ps), tax_table(ps))

    for (valmodel in valmodel) {

      path2 <- paste0('splsda_', var, '_', valmodel)

      setwd(path)
      dir.create(path2)
      setwd(paste0(path, path2))

      x3 <- ( phyloseq::otu_table(ps) %>% as.data.frame() ) + 1                        # create otu table

      set.seed(5)

      train <- sample(nrow(phyloseq::sample_data(ps)), 0.7*nrow(phyloseq::sample_data(ps)), replace = F) # option 1 - split ps 70:30 for training and validation
      trainset.all <- phyloseq::sample_data(ps)[train,]
      validset.all <- phyloseq::sample_data(ps)[-train,]

      trainset.msia <- phyloseq::subset_samples(ps, country != 'malaysia') %>% phyloseq::sample_data  # option 2 - split ps into immigrant and mainland communities
      validset.msia <- phyloseq::subset_samples(ps, country == 'malaysia') %>% phyloseq::sample_data

      trainset.call <- ifelse(valmodel == 'all', 'trainset.all', ifelse(valmodel == 'malaysia', 'trainset.msia', 'none defined'))
      validset.call <- ifelse(valmodel == 'all', 'validset.all', ifelse(valmodel == 'malaysia', 'validset.msia', 'none defined'))

      trainset <- trainset.call %>% get
      validset <- validset.call %>% get

      yt <- trainset[[var]]                 # var for training set
      yv <- validset[[var]]                 # for for validation set

      no.var <- yt %>% unique %>% length

      xt <- subset(x3, rownames(x3) %in% rownames(trainset))    # otu for training set
      xv <- subset(x3, rownames(x3) %in% rownames(validset))    # otu for validation set

      xt.clr <- mixOmics::logratio.transfo(xt, logratio = 'CLR') # transform abundance table with CLR


      plsda.t <- mixOmics::plsda(xt.clr, yt, ncomp = no.comp) # create plsda object
      plsda.t.perf <- mixOmics::perf(plsda.t, validation = 'Mfold', folds = 5, nrepeat = no.repeat, auc = T, progressBar = T) # cross-validation

      ## plot balanced error rate based on plsda model

      er <- plsda.t.perf$error.rate$BER[,1]
      component <- 1:no.comp
      er.df <- data.frame(component, er)

      fig.plsda.er <- ggplot2::ggplot(er.df, aes(x = component, y = er)) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(linetype = 'dashed', color = 'black', size = 0.1, se = F) +
        ggplot2::scale_x_discrete(limits = component) +
        ggplot2::theme_bw() +
        ggplot2::xlab('No. components') +
        ggplot2::ylab('error rate')

      ## plot AUC curve for plsda model

      best.comp.pre <- plsda.t.perf$choice.ncomp[2,1]             # take best no of components based on BER and max distance
      best.comp <- ifelse(best.comp.pre < 2, 2, best.comp.pre)    # make sure 2 is the minimum component no for plotting

      plsda.auc <- mixOmics::auroc(plsda.t, roc.comp = best.comp)

      fig.plsda.auc <- paste0('plsda.auc$graph.Comp', best.comp) %>% rlang::parse_expr %>% eval +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = 'bottom') +
        ggplot2::guides(color = guide_legend(nrow = 6)) #AUC plot

      ## plot plsda ordination for the first 2 axes

      plsda.plot <- mixOmics::plotIndiv(plsda.t, comp = c(1,2),
                                        style = 'ggplot2',
                                        group = yt,
                                        ind.names = F,
                                        ellipse = F,
                                        legend = T,
                                        title = 'PLS-DA')

      fig.plsda.plot <- plsda.plot$graph + ggplot2::theme_bw() + ggplot2::theme(legend.position = 'bottom')

      # validation test via confusion matrix

      xv.clr <- mixOmics::logratio.transfo(xv, logratio = 'CLR') # transform abundance table with CLR

      plsda.test <- mixOmics::predict(plsda.t, newdata = xv.clr, method = 'max.dist')
      plsda.test2 <- plsda.test$class$max.dist %>% as.data.frame
      table(yv, plsda.test2[,best.comp]) %>% prop.table(margin = 1)

      # print error rate to document

      sink('output2.txt', append = F)

      paste('PLSDA and SPLSDA analysis for', var) %>% print

      paste('physeq object used is') %>% print
      ps %>% print

      paste('this test used', no.comp, 'number of components, repeated for', no.repeat, 'times') %>% print

      paste('error rate for plsda model') %>% print
      plsda.t.perf$error.rate %>% print

      paste('best number of component based on t-test') %>% print
      best.comp %>% print

      paste('plsda confusion matrix based on best no of component') %>% print
      table(yv, plsda.test2[,best.comp]) %>% print

      paste('plsda confusion matrix in proportion') %>% print
      table(yv, plsda.test2[,best.comp]) %>% prop.table(margin = 1) %>% print

      sink()

      ### splsda

      splsda.t <- mixOmics::tune.splsda(xt.clr, 						# tune splsda
                                        Y = yt,
                                        ncomp = no.comp,
                                        multilevel = NULL,
                                        test.keepX = c(seq(5,150, 5)),
                                        validation = c('Mfold'),
                                        folds = 5,
                                        dist = 'max.dist',
                                        nrepeat = no.repeat,
                                        cpus = 8,
                                        auc = T,
                                        progressBar = T)

      best.comp2.pre <- splsda.t$choice.ncomp$ncomp				            # see no. optimal component based on t test
      best.comp2 <- ifelse(best.comp2.pre < 2, 2, best.comp2.pre)     # ensure best comp > 1
      choice.keepx <- splsda.t$choice.keepX[1:best.comp2]		          # use the optimal number of components based on best.comp2

      splsda.t2 <- mixOmics::splsda(X = xt.clr, Y = yt, ncomp = best.comp2, keepX = choice.keepx) # run splsda

      splsda.auc <- mixOmics::auroc(splsda.t2, roc.comp = best.comp2)
      fig.splsda.auc <- paste0('splsda.auc$graph.Comp', best.comp2) %>% rlang::parse_expr %>% eval +
        ggplot2::theme_bw() +
        ggplot2::theme(legend.position = 'bottom') +
        ggplot2::guides(color = guide_legend(nrow = 6))

      splsda.plot <- mixOmics::plotIndiv(splsda.t2,
                                         ind.names = F,
                                         col.per.group = color.mixo(1:no.var),
                                         comp = c(1,2),
                                         ellipse = F,
                                         legend = TRUE,
                                         star = F,
                                         centroid = F,
                                         title = 'sPLS-DA')

      fig.splsda.plot <- splsda.plot$graph + ggplot2::theme_bw() + ggplot2::theme(legend.position = 'bottom')

      # splsda cross validation

      splsda.cv <- mixOmics::perf(splsda.t2, validation = 'Mfold', folds = 5, progressBar = T, nrepeat = no.repeat)

      splsda.er <- splsda.cv$error.rate$BER[,1] # splsda error rate based on BER
      splsda.comp <- 1:best.comp2
      splsda.er.df <- data.frame(splsda.comp, splsda.er)

      fig.splsda.er <- ggplot2::ggplot(splsda.er.df, aes(x = splsda.comp, y = splsda.er)) +
        ggplot2::geom_point() +
        ggplot2::geom_smooth(linetype = 'dashed', color = 'black', size = 0.1, se = F) +
        ggplot2::scale_x_discrete(limits = splsda.comp) +
        ggplot2::theme_bw() +
        ggplot2::xlab('No. components') +
        ggplot2::ylab('error rate')

      # validation test via confusion matrix

      splsda.test <- mixOmics::predict(splsda.t2, newdata = xv.clr, method = 'max.dist', dist)
      splsda.test2 <- splsda.test$class$max.dist %>% as.data.frame
      table(yv, splsda.test2[,best.comp2])

      # sink splsda profile

      sink('output2.txt', append = T)

      paste('error rate for splsda model') %>% print
      splsda.cv$error.rate %>% print

      paste('best no of component based on splsda model') %>% print
      best.comp2 %>% print

      paste('confusion matrix based on splsda model') %>% print
      table(yv, splsda.test2[,best.comp2]) %>% print

      paste('splsda confusion matrix in proportion') %>% print()
      table(yv, splsda.test2[,best.comp2]) %>% prop.table(margin = 1) %>% print

      sink()

      # compile plots

      figcomp.plsda1 <- ggpubr::ggarrange(fig.plsda.er, fig.plsda.auc, labels = 'auto')
      figcomp.plsda2 <- ggpubr::ggarrange(figcomp.plsda1, fig.plsda.plot, ncol = 1, labels = c('', 'c'))
      ggsave('fig_plsda.pdf', figcomp.plsda2, units = 'mm', height = 250, width = 250)

      figcomp.splsda1 <- ggpubr::ggarrange(fig.splsda.er, fig.splsda.auc, labels = 'auto')
      figcomp.splsda2 <- ggpubr::ggarrange(figcomp.splsda1, fig.splsda.plot, ncol = 1, labels = c('', 'c'))
      ggsave('fig_splsda.pdf', figcomp.splsda2, units = 'mm', height = 250, width = 250)

      cowplot::plot_grid(fig.plsda.plot, fig.splsda.plot, fig.plsda.er, fig.splsda.er,
                         fig.plsda.auc, fig.splsda.auc, labels = c('auto'), nrow = 3, rel_heights = c(1.5,1,2))

      ggplot2::ggsave('fig_all.pdf', units = 'mm', height = 250, width = 250)


      # contribution plot -- mainly for mainland-immigrant model output to see mutual taxa differentiating both groups

      if (paste0(var, '.', valmodel) == 'country2.malaysia') {

        save.image(file = paste0('splsda_', var, '_', valmodel, '.RData'))

        setwd(path)

      } else {

        splsda.v <- mixOmics::tune.splsda(xv.clr, 						# tune splsda
                                          Y = yv,
                                          ncomp = no.comp,
                                          multilevel = NULL,
                                          test.keepX = c(seq(5,150, 5)),
                                          validation = c('Mfold'),
                                          folds = 5,
                                          dist = 'max.dist',
                                          nrepeat = no.repeat,
                                          cpus = 8,
                                          auc = T,
                                          progressBar = T)

        best.comp3.pre <- splsda.v$choice.ncomp$ncomp				            # see no. optimal component based on t test
        best.comp3 <- ifelse(best.comp3.pre < 2, 2, best.comp3.pre)     # ensure best comp > 1
        choice.keepx2 <- splsda.v$choice.keepX[1:best.comp3]		          # use the optimal number of components based on best.comp2

        splsda.v2 <- mixOmics::splsda(X = xv.clr, Y = yv, ncomp = best.comp3, keepX = choice.keepx2) # run splsda

        contib.t <- lapply(1:best.comp2, function(x){                   # generate contrib dataframe for training set
          mixOmics::plotLoadings(splsda.t2,
                                 comp = x,
                                 method = 'mean',
                                 contrib = 'max')})

        contib.t2 <- lapply(1:best.comp2, function(x) cbind(contib.t[[x]], x, contib.t[[x]] %>% rownames))
        contib.t3 <- do.call(rbind, contib.t2)
        contib.t4 <- contib.t3 %>% dplyr::arrange(dplyr::desc(abs(importance))) %>% dplyr::group_by(x) %>% dplyr::slice(1:7)
        contib.t4$genus <- contib.t4[['contib.t[[x]] %>% rownames']]
        contib.t4$group = 'training'

        contib.v <- lapply(1:best.comp3, function(x){                   # generate contrib dataframe for validation set
          mixOmics::plotLoadings(splsda.v2,
                                 comp = x,
                                 method = 'mean',
                                 contrib = 'max')})

        contib.v2 <- lapply(1:best.comp3, function(x) cbind(contib.v[[x]], x, contib.v[[x]] %>% rownames))
        contib.v3 <- do.call(rbind, contib.v2)
        contib.v4 <- contib.v3 %>% dplyr::arrange(dplyr::desc(abs(importance))) %>% dplyr::group_by(x) %>% dplyr::slice(1:7)
        contib.v4$genus <- contib.v4[['contib.v[[x]] %>% rownames']]
        contib.v4$group = 'validation'

        contib.all <- rbind(contib.t4, contib.v4)                       # combine both training and validation contrib dataframe

        fig.contib <- ggplot2::ggplot(contib.all, aes(x = genus, y = importance %>% abs, fill = GroupContrib %>% as.factor )) +
          ggplot2::geom_col(position = 'dodge') +
          ggplot2::coord_flip() +
          ggplot2::facet_grid(~paste0(group)) +
          ggplot2::guides(color = guide_legend(nrow = 1)) +
          ggplot2::theme_bw() +
          ggplot2::labs(fill = 'component no.') +
          ggplot2::theme(legend.position = 'bottom', legend.title = element_blank()) +
          ggplot2::xlab(element_blank()) +
          ggplot2::ylab('importance') +
          ggsci::scale_fill_d3()

        ggplot2::ggsave('fig_contrib.pdf', units = 'mm', height = 250, width = 250)


        save.image(file = paste0('splsda_', var, '_', valmodel, '.RData'))

        setwd(path)

      }
    }
  }
}
