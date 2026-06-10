wisca_plot <- function(params_plot,cover_plot, pathogen_in, analysis, infection_full){
  
  
  # Epidemiological and WISCA visualisations
  if (analysis %in% names(params_plot)){
    plot_list <- list()
    wisca_ci <- cover_plot %>%
      mutate(CIinterval = paste0(Lower_CI,"-",Upper_CI))
    names(wisca_ci) <- c("Regimen_full","Regimen",str_to_title(analysis),
                         "Coverage","Lower credible bound","Upper credible bound","95% Credible interval")
    wisca_ci$Regimen <- as.factor(wisca_ci$Regimen)
    
    for (pp in 1:length(levels(params_plot[[analysis]]))){


      plot_set <- subset(params_plot, params_plot[[analysis]] == levels(params_plot[[analysis]])[pp])
      cover_set <- subset(wisca_ci, wisca_ci[[str_to_title(analysis)]] == levels(params_plot[[analysis]])[pp])
      
      N <- sum(distinct(plot_set,mo,.keep_all = TRUE)$n.x) 
      plot_set$perc.n <- round(plot_set$prop.n*100,2)
      plot_set$perc.S <- round(plot_set$perc.S, 2)
      plot_set <- arrange(plot_set, by = "perc.n")
      plot_set$fullname <- str_to_sentence(plot_set$fullname)
      
      plot_display <- plot_set %>%
        select(all_of(analysis),mo:n.x, perc.n,regimens,keyantimicrobials:S_n,perc.S) 
      names(plot_display) <- c(str_to_title(analysis),"Microorganism code","Microorganism fullname", "Pathogen incident count", 
                               "Percentage incidence","Antibiotic regimens","Abbrev Antibiotic regimens","Tested count","Pathogen sensitive count",
                               "Percentage pathogen sensitive")
      cover_plot_2 <- subset(cover_plot, cover_plot[[str_to_title(analysis)]] == levels(params_plot[[analysis]])[pp])
      cover_plot_2 <- cover_plot_2 %>%
        rename("Abbrev Antibiotic regimens" = "Regimen")
      
      plot_display <- left_join(x = cover_plot_2[,c("Abbrev Antibiotic regimens","Coverage","Lower credible bound","Upper credible bound","95% Credible interval")], y = plot_display, by = c(str_to_title(analysis), "Abbrev Antibiotic regimens"))
      plot_display <- plot_display %>%
        mutate(`Infection type` = infection_full) %>%
        select(`Antibiotic regimens`, Coverage:`Percentage incidence`,`Tested count`:`Infection type`)
      
      plot_list[[pp]] <- plot_display 
      
      plot_sens_list <- list()
      for (i in 1:length(unique(plot_set$keyantimicrobials))){
                ## Pathogen sensitivity to antimicrobials
        params_sens_plot <- plot_set %>%
          filter(keyantimicrobials == levels(cover_plot$Regimen)[i]) 
        
        if (analysis_level == "ward"){
        plot_sens <- ggplot(params_sens_plot, aes(x = fct_infreq(fullname)))+ geom_col(aes(y = perc.S), fill = "red4") +
          facet_wrap(facets = vars(regimens), nrow = ceiling(length(levels(plot_set$keyantimicrobials))/2), ncol = 2) +
          theme_bw() +
          theme(text = element_text(family = "serif", size = 14),
                axis.text.y = element_text(face = "italic")) +
          labs(title = paste0("Pathogen susceptibility distribution to ",unique(params_sens_plot$regimens)," in ", levels(params_plot[[analysis]])[pp]),
               subtitle = paste0("Facility: ",cover_plot$hospital),
               y = paste0("Pathogen susceptibility "," (N = ",N, ")"),
               x = "Pathogen classification")  + 
          scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
          coord_flip()
        } else {
          plot_sens <- ggplot(params_sens_plot, aes(x = fct_infreq(fullname)))+ geom_col(aes(y = perc.S), fill = "red4") +
            facet_wrap(facets = vars(regimens), nrow = ceiling(length(levels(plot_set$keyantimicrobials))/2), ncol = 2) +
            theme_bw() +
            theme(text = element_text(family = "serif", size = 14),
                  axis.text.y = element_text(face = "italic")) +
            labs(title = paste0("Pathogen susceptibility distribution to ",unique(params_sens_plot$regimens)," in ", levels(params_plot[[analysis]])[pp]),
                 y = paste0("Pathogen susceptibility "," (N = ",N, ")"),
                 x = "Pathogen classification")  + 
            scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
            coord_flip()
        }
        
        plot_sens_list[[i]] <- plot_sens
       
      }
      
      if(analysis_level == "ward"){
      
        
        plot_cover <- ggplot(cover_plot_2, aes(x = Regimen_full, y = Coverage)) + geom_point(size = 3, colour = "red") + 
        geom_text(aes(label = paste0(ceiling(Coverage)," %")),nudge_x = 0.1, size = 5) +
        geom_errorbar(aes(ymin = `Lower credible bound`, ymax = `Upper credible bound`), colour = "purple",width = 0.1, linewidth = 1) +
        theme_bw() +
        theme(text = element_text(family = "serif", size = 14),
              axis.text = element_text(size = 14), axis.text.x = element_text(angle = 0)) + 
        labs(title = paste0("Antibiotic coverage point estimates for selected regimens (with 95% credible intervals)"),
             subtitle = paste0("Facility: ",cover_plot_2$hospital),
             x = "Antibiotic regimens",
             y = "Coverage (95% credible intervals) (%)") +
        scale_y_continuous(expand = c(0,0), limits = c(0,100)) 
      } else {
    
        plot_cover <- ggplot(cover_plot_2, aes(x = Regimen_full, y = Coverage)) + geom_point(size = 3, colour = "red") + 
          geom_text(aes(label = paste0(ceiling(Coverage)," %")),nudge_x = 0.1, size = 5) +
          geom_errorbar(aes(ymin = `Lower credible bound`, ymax = `Upper credible bound`), colour = "purple",width = 0.1, linewidth = 1) +
          theme_bw() +
          theme(text = element_text(family = "serif", size = 14),
                axis.text = element_text(size = 14), axis.text.x = element_text(angle = 0)) + 
          labs(title = paste0("Antibiotic coverage point estimates for selected regimens (with 95% credible intervals)"),
               x = "Antibiotic regimens",
               y = "Coverage (95% credible intervals) (%)") +
          scale_y_continuous(expand = c(0,0), limits = c(0,100)) 
      }
      
      #if want to plot sensitivity and cover together
      plot_epi <- ggpubr::ggarrange(ggpubr::ggarrange(plotlist = plot_sens_list, ncol = 1, nrow = length(plot_sens_list)),
                                    plot_cover, 
                                    ncol = 2, labels = c("A","B"), widths = c(1,1))
      
      print(plot_epi)
      
      
      #if just want to group the sensitivity plots together and print cover separately
      plot_sens <- ggpubr::ggarrange(ggpubr::ggarrange(plotlist = plot_sens_list, ncol = 1, nrow = length(plot_sens_list)))
      print(plot_sens)
      print(plot_cover)
    
    }
    names(plot_list) <- levels(params_plot[[analysis]])
    return(plot_list) #nicely formatted table
    
  } else {
    
    wisca_ci <- cover_plot %>%
      mutate(CIinterval = paste0(Lower_CI,"-",Upper_CI))
    names(wisca_ci) <- c("Regimen_full","Regimen","Coverage","Lower credible bound","Upper credible bound","95% Credible interval")
    wisca_ci$Regimen <- as.factor(wisca_ci$Regimen)
    
    N <- sum(distinct(params_plot,mo,.keep_all = TRUE)$n.x) 
    params_plot$perc.n <- round(params_plot$prop.n*100,2)
    params_plot$perc.S <- round(params_plot$perc.S, 2)
    params_plot <- arrange(params_plot, by = "perc.n")
    params_plot$fullname <- str_to_sentence(params_plot$fullname)
    
    plot_display <- params_plot %>%
      select(mo:n.x, perc.n,regimens,keyantimicrobials:S_n,perc.S) 
    names(plot_display) <- c("Microorganism code","Microorganism fullname", "Pathogen incident count", 
                             "Percentage incidence","Antibiotic regimens","Abbrev Antibiotic regimens","Tested count","Pathogen sensitive count",
                             "Percentage pathogen sensitive")
    
    
    cover_plot_2 <- wisca_ci %>% 
      rename("Abbrev Antibiotic regimens" = "Regimen")
    
    plot_display <- left_join(x = cover_plot_2[,c("Abbrev Antibiotic regimens","Coverage","Lower credible bound","Upper credible bound","95% Credible interval")], y = plot_display, by = c("Abbrev Antibiotic regimens"))
    plot_display <- plot_display %>%
      mutate(`Infection type` = infection_full) %>%
      select(`Antibiotic regimens`, Coverage:`Percentage incidence`,`Tested count`:`Infection type`)
    
    plot_sens_list <- list()
    for (i in 1:length(unique(cover_plot$Regimen))){

      ## Pathogen sensitivity to antimicrobials
      params_sens_plot <- params_plot %>%
        filter(keyantimicrobials == levels(cover_plot$Regimen)[i]) 
      
      plot_sens <- ggplot(params_sens_plot, aes(x = fullname))+ geom_col(aes(y = perc.S), fill = "red4") +
        theme_bw() +
        theme(text = element_text(family = "serif", size = 14),
              axis.text.y = element_text(face = "italic")) +
        labs(title = paste0("Pathogen susceptibility distribution to ",unique(params_sens_plot$regimens)),
             y = paste0("Pathogen susceptibility"," (N = ",N, ")"),
             x = "Pathogen classification")  + 
        scale_y_continuous(expand = c(0,0), limits = c(0,100)) +
        coord_flip()
      
      plot_sens_list[[i]] <- plot_sens
      
      }
    
    
    plot_cover <- ggplot(wisca_ci, aes(x = fct_inorder(Regimen_full), y = Coverage)) + geom_point(size = 3, colour = "red") + 
      geom_text(aes(label = paste0(ceiling(Coverage)," %")),nudge_x = 0.1, size = 5) +
      geom_errorbar(aes(ymin = `Lower credible bound`, ymax = `Upper credible bound`), colour = "purple",width = 0.1, linewidth = 1) +
      theme_bw() +
      theme(text = element_text(family = "serif", size = 14),
            axis.text = element_text(size = 14), axis.text.x = element_text(angle = 0)) + 
      labs(title = paste0("Antibiotic coverage point estimates for selected regimens (with 95% credible intervals)"),
           x = "Antibiotic regimens",
           y = "Coverage (95% credible intervals)(%)") +
      scale_y_continuous(expand = c(0,0), limits = c(min(wisca_ci$`Lower credible bound`)-10,100))
    
#if want to combine sensitivity and cover
    plot_epi <- ggpubr::ggarrange(ggpubr::ggarrange(plotlist = plot_sens_list, ncol = 1, nrow = length(plot_sens_list)),
                                  plot_cover, 
                                  ncol = 2, labels = c("A","B"), widths = c(1,1))
    print(plot_epi)
    
    #if just want to group the sensitivity plots together and print cover separately
    plot_sens <- ggpubr::ggarrange(ggpubr::ggarrange(plotlist = plot_sens_list, ncol = 1, nrow = length(plot_sens_list)))
    print(plot_sens)
    print(plot_cover)
    
    return(plot_display) #nicely formatted table
      }
  
  
}
  
  
  