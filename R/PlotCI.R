PlotCI <- function(output, title="95% SCIs", background="white"){
  
  if(class(output)=="silo"){
    output <- list(Fantasyname=output)
  }
  
  if(length(output)!=1){
    if(is.null(names(output))){
      stop("output must be a NAMED list!")
    }
  }
  
  background <- match.arg(background, choices=c("gray", "grey", "white"))
  
  ### Extract the results
  
  results <- list()
  
  for(l in 1:length(output)){
    results[[l]] <- output[[l]]$Results
  }
  
  names(results) <- names(output)
  
  ###
  
  for(l in 1:length(results)){
    results[[l]]$Comp <- factor(rownames(results[[l]]), levels=rev(rownames(results[[l]])))
    if(length(results)==1){
      results[[l]]$Method <- "Fantasyname"
    }else{
      results[[l]]$Method <- names(results)[[l]]      
    }
  }
  
  if(length(results)==1){
    tab <- rbind(results[[1]])
    colval <- c("black")
  }
  
  if(length(results)==2){
    tab <- rbind(results[[1]], results[[2]])
    colval <- c("black", "gray50")
  }
  
  if(length(results)==3){
    tab <- rbind(results[[1]], results[[2]], results[[3]])
    colval <- c("black", "gray50", "gray80")
  }
  
  if(length(results)==4){
    tab <- rbind(results[[1]], results[[2]], results[[3]], results[[4]])
    colval <- c("black", "gray50", "gray80", "gray20")
  }
  
  if(length(results)==5){
    tab <- rbind(results[[1]], results[[2]], results[[3]], results[[4]], results[[5]])
    colval <- c("black", "gray50", "gray80", "gray20", "gray65")
  }
  
  if(length(results)==6){
    tab <- rbind(results[[1]], results[[2]], results[[3]], results[[4]], results[[5]], results[[6]])
    colval <- c("black", "gray50", "gray80", "gray20", "gray65", "gray35")
  }
  
  if(length(results)==7){
    tab <- rbind(results[[1]], results[[2]], results[[3]], results[[4]], results[[5]], results[[6]], results[[7]])
    colval <- c("black", "gray50", "gray80", "gray20", "gray65", "gray35", "white")
  }
  
  pd <- position_dodge(width=0.5, height=NULL)
  
  if(length(results)==1){
    lp <- "none"
    tit <- paste(title, "\n", sep="")
  }else{
    lp <- "top"
    tit <- title
  }
  
  if(background=="white"){
    backg <- theme_bw()
  }else{
    backg <- theme_grey()
  }
  
  if(tab$Upper[1]=="Inf"){
    niceCIplot <- ggplot(tab, aes(x=Comp, y=Estimate, ymax=Estimate, color=Method)) +
      geom_hline(yintercept=0, colour=gray(0.5), lty=2) +
      geom_errorbar(aes(ymin=Lower, ymax=Upper), lwd=0.8, width=0, position=pd) +
      geom_point(aes(y=Lower), size=7, shape="|", position=pd) +
      geom_point(size=4, position=pd) +
      scale_color_manual(name="Method", values=colval) + 
      coord_flip() +
      xlab("") +
      ggtitle(tit) +
      backg +
      theme(legend.position=lp,
            legend.title=element_blank(),
            legend.text=element_text(size=13),
            legend.key=element_rect(colour="white"),
            plot.title=element_text(size=rel(2)),
            axis.text.x=element_text(size=13),
            axis.text.y=element_text(size=13),
            axis.title.x=element_text(size=15))
  }else{
    if(tab$Lower[1]=="-Inf"){
      niceCIplot <- ggplot(tab, aes(x=Comp, y=Estimate, ymax=Estimate, color=Method)) +
        geom_hline(yintercept=0, colour=gray(0.5), lty=2) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), lwd=0.8, width=0, position=pd) +
        geom_point(aes(y=Upper), size=7, shape="|", position=pd) +
        geom_point(size=4, position=pd) +
        scale_color_manual(name="Method", values=colval) + 
        coord_flip() +
        xlab("") +
        ggtitle(tit) +
        backg +
        theme(legend.position=lp,
              legend.title=element_blank(),
              legend.text=element_text(size=13),
              legend.key=element_rect(colour="white"),
              plot.title=element_text(size=rel(2)),
              axis.text.x=element_text(size=13),
              axis.text.y=element_text(size=13),
              axis.title.x=element_text(size=15))
    }else{
      niceCIplot <- ggplot(tab, aes(x=Comp, y=Estimate, ymax=Estimate, color=Method)) +
        geom_hline(yintercept=0, colour=gray(0.5), lty=2) +
        geom_errorbar(aes(ymin=Lower, ymax=Upper), lwd=0.8, width=0.5, position=pd) +
        geom_point(size=4, position=pd) +
        scale_color_manual(name="Method", values=colval) + 
        coord_flip() +
        xlab("") +
        ggtitle(tit) +
        backg +
        theme(legend.position=lp,
              legend.title=element_blank(),
              legend.text=element_text(size=13),
              legend.key=element_rect(colour="white"),
              plot.title=element_text(size=rel(2)),
              axis.text.x=element_text(size=13),
              axis.text.y=element_text(size=13),
              axis.title.x=element_text(size=15))
    }
  }
  
  print(niceCIplot)
  
}