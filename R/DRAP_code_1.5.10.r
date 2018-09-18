
####################################################################################################
#' @name DataSummary
#' @title Data summary
#' @description Summarize the mean, SD, SE, and sample size of each arm in input data.
#' @usage  DataSummary(data, type, pattern = c("oneAN", "TAN"), measure.var, group.vars)
#'
#' @param data a data frame of measured tumor vloume data or body weight data of mouse.
#' @param type the type of data, "Vomule" or "BodyWeight".
#' @param pattern the pattern of PDX trial design, "oneAN" or "TAN".
#' @param measure.var the measured variable, such as "Volume", "BodyWeight".
#' @param group.vars the group variables.
#'
#' @import plyr
#'
#' @return A data frame including the sample size, mean, SD, SE of group variables.
#' 
#' @examples
#' ## summary the tumor volume data of oneAN pattern
#' data(oneAN.volume.data)
#' oneAN.volume.data[1:10,]
#' oneAN.v.s <- DataSummary(data = oneAN.volume.data,
#'                          type = 'Volume',
#'                          pattern = 'oneAN',
#'                          measure.var = 'Volume',
#'                          group.vars = c('Arms','Times'))
#' head(oneAN.v.s)
#' 
#' ## summary the body weight data of oneAN pattern
#' data(oneAN.bw.data)
#' oneAN.bw.data[1:10,]
#' oneAN.bw.s <-DataSummary(data = oneAN.bw.data,
#'                          type = 'BodyWeight',
#'                          pattern = 'oneAN',
#'                          measure.var = 'BodyWeight',
#'                          group.vars = c('Arms','Times'))
#' head(oneAN.bw.s)
#' 
#' \donttest{
#' ## summary the tumor volume data of TAN pattern
#' data(TAN.volume.data)
#' TAN.volume.data[1:10,]
#' TAN.v.s <- DataSummary(data = TAN.volume.data,
#'                        type = 'Volume',
#'                        pattern = 'TAN',
#'                        measure.var = 'Volume',
#'                        group.vars = c('Tumor','Arms','Times'))
#' head(TAN.v.s)
#' 
#' ## summary the body weight data of TAN pattern
#' data(TAN.bw.data)
#' TAN.bw.data[1:10,]
#' TAN.bw.s <- DataSummary(data = TAN.bw.data,
#'                         type = 'BodyWeight',
#'                         pattern = 'TAN',
#'                         measure.var = 'BodyWeight',
#'                         group.vars = c('Tumor','Arms','Times'))
#' head(TAN.bw.s)
#' }
#' 
#' @export
DataSummary <- function(data, type, pattern = c('oneAN','TAN'), measure.var, group.vars){
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  if(pattern == 'oneAN'){
    
    # make sure the colnames of indata including 'Arms', 'ID', 'Times', 'Volume' or 'BodyWeight', which guarantee normal running for the following code.
    
    if (type == 'Volume'){
      need.cols <- c('Arms', 'ID', 'Times', 'Volume')
      judged.res <- is.element(need.cols,colnames(data))
      
      if ("FALSE" %in% judged.res)
        stop("The input data must include at least four columns, that are 'Arms', 'ID', 'Times' and 'Volume'. Please check the data whether includes these columns and ensure the colnames is consistent with 'Arms', 'ID', 'Times' and 'Volume'")
    }
    
    if (type == 'BodyWeight'){
      need.cols <- c('Arms', 'ID', 'Times','BodyWeight')
      judged.res <- is.element(need.cols,colnames(data))
      
      if ("FALSE" %in% judged.res)
        stop("The input data must include at least four columns, that are 'Arms', 'ID', 'Times' and 'BodyWeight'. Please check the data whether includes these columns and ensure the colnames is consistent with 'Arms', 'ID', 'Times' and 'BodyWeight'")
    }
  }
  
  
  if(pattern == 'TAN'){
    
    # make sure the colnames of indata including 'Arms', 'ID', 'Times', 'Volume' or 'BodyWeight', which guarantee normal running for the following code.
    
    if (type == 'Volume'){
      need.cols <- c('Tumor','Arms', 'ID', 'Times', 'Volume')
      judged.res <- is.element(need.cols,colnames(data))
      
      if ("FALSE" %in% judged.res)
        stop("The input data must include at least four columns, that are 'Tumor', 'Arms', 'ID', 'Times' and 'Volume'. Please check the data whether includes these columns and ensure the colnames is consistent with 'Tumor', 'Arms', 'ID', 'Times' and 'Volume'")
    }
    
    if (type == 'BodyWeight'){
      need.cols <- c('Tumor', 'Arms', 'ID', 'Times','BodyWeight')
      judged.res <- is.element(need.cols,colnames(data))
      
      if ("FALSE" %in% judged.res)
        stop("The input data must include at least four columns, that are 'Tumor', 'Arms', 'ID', 'Times' and 'BodyWeight'. Please check the data whether includes these columns and ensure the colnames is consistent with 'Tumor', 'Arms', 'ID', 'Times' and 'BodyWeight'")
    }
  }
  
  data <- subset(data,measure.var>0)
  
  data.summary <- ddply(.data = data,.variables = group.vars,
                        .fun = function(xx, col) {
                          c(N    = sum(!is.na(xx[[col]])),
                            Mean = mean  (xx[[col]], na.rm=T),
                            SD   = sd    (xx[[col]], na.rm=T)
                          )},
                        measure.var
  )
  
  data.summary <- plyr::rename(data.summary, c("Mean" = measure.var))
  data.summary$SE <- data.summary$SD/sqrt(data.summary$N)  # Calculate standard error
  
  return(data.summary)
}


####################################################################################################
#' @name RelativeChange
#' @title Caculate the relative change for tumor volume data or body weight data.
#' @description Caculate the relative change for tumor volume data or body weight data.
#' @usage RelativeChange(data,type,rm.baseline = TRUE)
#'
#' @param data a data frame of measured tumor vloume data or body weight data of mouse.
#' @param type the type of data, "Vomule" or "BodyWeight".
#' @param rm.baseline a logical parameter, whether remove the baseline.
#'
#' @details If rm.baseline = T, the way of caculating follows like Vt/V0 - 1; if rm.baseline = F, the way of caculating follows like Vt/V0.
#'
#' @return If rm.baseline = T, it returns the relative change of tumor volume data or body weight data; if rm.baseline = F, it returns the fold change of tumor volume data or body weight data.
#' 
#' @examples 
#' data(oneAN.volume.data)
#' ## caculate relative change of tumor volume
#' oneAN.v.rc <- RelativeChange(oneAN.volume.data,type = 'Volume',rm.baseline = TRUE)
#' head(oneAN.v.rc)
#' 
#' ##caculate fold change of tumor volume
#' oneAN.v.fc <- RelativeChange(oneAN.volume.data,type = 'Volume',rm.baseline = FALSE)
#' head(oneAN.v.fc)
#' 
#' ##caculate relative change of body weight data
#' data(oneAN.bw.data)
#' oneAN.bw.rc <- RelativeChange(oneAN.bw.data,type = 'BodyWeight',rm.baseline = TRUE)
#' head(oneAN.bw.rc)
#'
#' @export
RelativeChange <- function(data, type, rm.baseline = TRUE){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  if (type == 'Volume'){
    Volume <- data[ ,'Volume']
    data <- subset(data,Volume > 0)
  }

  if (type == 'BodyWeight'){
    BodyWeight <- data[ ,'BodyWeight']
    data <- subset(data,BodyWeight > 0)
  }

  rc.res <- vector()
  ID <- unique(data$ID)

  for (i in 1:length(ID)){
    s.data <- data[data$ID == ID[i],]

    if (type == 'Volume'){ s.rc <- s.data$Volume/s.data$Volume[1] * 100 }
    if (type == 'BodyWeight'){ s.rc <- s.data$BodyWeight/s.data$BodyWeight[1] * 100 }

    s.rc <- data.frame(ID = s.data$ID,Times = s.data$Times,
                       RelativeChange  = s.rc,stringsAsFactors = F)
    rc.res <- rbind(rc.res,s.rc)
  }

  if (type == 'Volume'){ in.data=data[,setdiff(colnames(data),'Volume')]}
  if (type == 'BodyWeight'){in.data=data[,setdiff(colnames(data),'BodyWeight')] }

  rc.output <- merge(in.data,rc.res,by = c('ID','Times'))
  rc.output <- rc.output[,c(colnames(in.data),'RelativeChange')]

  if(rm.baseline) {
    rc.output$RelativeChange <- rc.output$RelativeChange - 100
  }

  return(rc.output)
}


####################################################################################################
#' @name plotVolumeGC
#' @title Plot growth curves of tumor volume data
#' @description Plot growth curves of tumor volume data.
#' @usage 
#' plotVolumeGC(data, level = c('Animal','Arm'), 
#'              pattern = c("OneAN", "TAN"), orders = NULL, 
#'              position.dodge)
#' 
#' @param data a data frame of measured tumor volume data.
#' @param level the level to present, "Animal" level or "Arm" level.
#' @param pattern the pattern of PDX trial design, "OneAN" or "TAN".
#' @param orders the redefined order for visulization. For oneAN pattern,the order is redifined orders of the 'Arms'; for TAN pattern, the order is redifined orders of the 'Tumor'.
#' @param position.dodge Dodging preserves the vertical position of an geom while adjusting the horizontal position.
#' 
#' @import ggplot2
#' 
#' @examples 
#' 
#' ### oneAN pattern
#' data(oneAN.volume.data)
#' # Animal level
#' plotVolumeGC(data = oneAN.volume.data,
#'              level = 'Animal',pattern = 'oneAN')
#' 
#' # Arm level
#' plotVolumeGC(data = oneAN.volume.data,level = 'Arm',
#'              pattern = 'oneAN',position.dodge = 0.5)
#' 
#' ## reset order of Arms
#' orders <- unique(oneAN.volume.data$Arms)
#' orders <- orders[c(2:6,1)]
#' plotVolumeGC(data = oneAN.volume.data,level = 'Animal',
#'              pattern = 'oneAN',orders = orders)
#' 
#' \donttest{
#' ### TAN pattern
#' data(TAN.volume.data)
#' # Animal level
#' plotVolumeGC(data = TAN.volume.data,level = 'Animal',
#'              pattern = 'TAN')
#' 
#' # Arm level
#' plotVolumeGC(data = TAN.volume.data,level = 'Arm',
#'              pattern = 'TAN',position.dodge = 0.5)
#' 
#' ## reset order of tumors
#' orders <- unique(TAN.volume.data$Tumor)
#' orders <- orders[c(2:4,1)]
#' plotVolumeGC(data = TAN.volume.data,level = 'Animal',
#'              pattern = 'TAN',orders = orders)
#' plotVolumeGC(data = TAN.volume.data,level = 'Arm',
#'              pattern = 'TAN',orders = orders,position.dodge = 0.5)
#' }
#' 
#' @export
plotVolumeGC <- function(data, level = c('Animal','Arm'), pattern = c("OneAN", "TAN"), orders = NULL, position.dodge){

  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  #check the order
  if (!is.null(orders)){

    if(pattern == 'oneAN'){
      arms <- unique(data$Arms)
      judged.order <- is.element(orders,arms)

      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Arms' among data.")

      data$Arms <- factor(data$Arms,levels=orders) # set the order
    }

    if(pattern == 'TAN'){
      tumors <- unique(data$Tumor)
      judged.order <- is.element(orders,tumors)

      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Tumor' among data.")

      data$Tumor <- factor(data$Tumor,levels=orders) # set the order
    }
  }

  Volume <-  data[ , 'Volume']
  data <- subset(data,Volume > 0)
  
  if( level == 'Animal'){
    
    Volume <-  data[ , 'Volume']
    Times<- data[ , 'Times']
    ID <- data[ , 'ID']
    Arms <- data[ , 'Arms']
    
    p <- ggplot(data, aes(x = Times, y = Volume,group = ID,color = Arms)) +
      geom_line(size=0.8)+
      geom_point(cex=1.5,aes(colour = Arms))
  }
  if( level == 'Arm'){
    if(pattern =='oneAN') s.data <- DataSummary(data,type = 'Volume',pattern = pattern, measure.var = "Volume", group.vars = c("Arms","Times"))
    if(pattern == 'TAN') s.data <- DataSummary(data,type = 'Volume',pattern = pattern, measure.var = "Volume", group.vars = c("Tumor","Arms","Times"))

    s.data <- s.data[s.data$N > 1,] #delete group without SD and SE
    
    Volume <-  s.data[ , 'Volume']
    Times<- s.data[ , 'Times']
    Arms <- s.data[ , 'Arms']
    SE <- s.data[ , 'SE']
    
    p <- ggplot(data = s.data, aes(x = Times, y = Volume,color = Arms)) +
      geom_line(position = position_dodge(position.dodge),cex = 1.2) +
      geom_errorbar(aes(ymin = Volume - SE, ymax = Volume + SE), 
                    width = 1.2,
                    position = position_dodge(position.dodge)) +
      geom_point(cex = 2,position = position_dodge(position.dodge))
  }

  p <- p + xlab("Time (d)") + ylab(expression(bold(paste("Tumor Volume (",mm^3,")",sep = " "))))

  p <- p + theme_bw() +
    theme(
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background  = element_rect(fill = "transparent")
    )   #backgroud
  p <- p + theme(
    axis.title.x = element_text(face = "bold",size = 12),
    axis.text.x  = element_text(vjust = 0,size = 12),
    axis.title.y = element_text(face = "bold",size = 12),
    axis.text.y  = element_text(hjust = 1,size = 12)
  )

  if( pattern == 'oneAN' & level == 'Animal'){
    p <- p + theme(legend.position='none')
    p <- p + facet_wrap( ~ Arms,dir = 'h') +
      theme(strip.text = element_text(colour = "black", face = "bold", size = rel(1)),
            strip.background = element_rect(fill = "#E6AB02", size = rel(1.05), linetype = 1)
      )
  }
  if(pattern == 'TAN'){
    
    if(level == 'Arm'){
      
      Tumor <- s.data[ , 'Tumor']
      
      p <- p + facet_wrap( ~ Tumor,dir = 'h') +
        theme(strip.text = element_text(colour = "black", face = "bold", size = rel(1)),
              strip.background = element_rect(fill = "#E6AB02", size = rel(1.05), linetype = 1)
        )
    }
    if(level == 'Animal'){
      
      Tumor <- data[ , 'Tumor']
      
      p <- p + theme(legend.position='none')
      p <- p + facet_grid(Arms ~ Tumor,space = "free_x", scales = "free_x") +
        theme(strip.text = element_text(colour = "black", face = "bold", size = rel(1)),
              strip.background = element_rect(fill = "#E6AB02",size = rel(1.05), linetype = 1)
        )
    }
  }

  plot(p)

}


####################################################################################################
#' @name plotBWeightGC
#' @title Plot growth curves for body weight of mouses
#' @description Plot growth curves for body weight of mouses.
#' @usage 
#' plotBWeightGC(data, level = c('Animal','Arm'), 
#'               pattern = c("OneAN", "TAN"), 
#'               orders = NULL, position.dodge)
#'
#' @param data a data frame of measured body weight data of mouse.
#' @param level the level to present, "Animal" level or "Arm" level.
#' @param pattern the pattern of PDX trial design, "OneAN" or "TAN".
#' @param orders the redefined order for visulization. For oneAN pattern,the order is redifined orders of the 'Arms'; for TAN pattern, the order is redifined orders of the 'Tumor'.
#' @param position.dodge Dodging preserves the vertical position of an geom while adjusting the horizontal position.
#'
#' @import ggplot2
#' 
#' @examples 
#' ### oneAN pattern
#' data(oneAN.bw.data)
#' plotBWeightGC(data = oneAN.bw.data,level = 'Animal',
#'               pattern = 'oneAN')
#' plotBWeightGC(data = oneAN.bw.data,level = 'Arm',
#'               pattern = 'oneAN',position.dodge = 0.5)
#' 
#' \donttest{
#' ### TAN pattern
#' data(TAN.bw.data) 
#' plotBWeightGC(data = TAN.bw.data,level = 'Animal',
#'               pattern = 'TAN')
#' plotBWeightGC(data = TAN.bw.data,level = 'Arm',
#'               pattern = 'TAN',position.dodge = 0.5)
#' }
#' 
#' @export
plotBWeightGC <- function(data, level = c('Animal','Arm'), pattern = c("OneAN", "TAN"), orders = NULL, position.dodge){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  #check the order
  if (!is.null(orders)){

    if(pattern == 'oneAN'){
      arms <- unique(data$Arms)
      judged.order <- is.element(orders,arms)

      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Arms' among data.")

      data$Arms <- factor(data$Arms,levels=orders) # set the order
    }

    if(pattern == 'TAN'){
      tumors <- unique(data$Tumor)
      judged.order <- is.element(orders,tumors)

      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Tumor' among data.")

      data$Tumor <- factor(data$Tumor,levels=orders) # set the order
    }
  }
  
  BodyWeight <- data[ , 'BodyWeight']
  data <- subset(data, BodyWeight > 0)
  
  if( level == 'Animal'){
    
    BodyWeight <- data[ , 'BodyWeight']
    Times <- data[ , 'Times']
    ID <- data[ , 'ID']
    Arms <-data[, 'Arms']
    
    p <- ggplot(data, aes(x = Times, y = BodyWeight,group = ID,color = Arms)) +
      geom_line(size = 0.8) +
      geom_point(cex = 1.5,aes(colour = Arms))
  }
  if( level == 'Arm'){
    if(pattern == 'oneAN') {
      s.data <- DataSummary(data,type = 'BodyWeight',pattern = pattern, measure.var = "BodyWeight", group.vars = c("Arms","Times"))
      
    }
    if(pattern == 'TAN') {
      s.data <- DataSummary(data,type = 'BodyWeight',pattern = pattern, measure.var = "BodyWeight", group.vars = c("Tumor","Arms","Times"))
    }

    s.data <- s.data[s.data$N > 1,] #delete group without SD and SE
    BodyWeight <- s.data[ , 'BodyWeight']
    Times <- s.data[ , 'Times']
    Arms <- s.data[, 'Arms']
    SE <- s.data[ , 'SE']

    p <- ggplot(data = s.data, aes(x = Times, y = BodyWeight,color = Arms))+
      geom_line(position = position_dodge(position.dodge),cex = 1.2) +
      geom_errorbar(aes(ymin = BodyWeight-SE, ymax = BodyWeight+SE), width = 1.2,position = position_dodge(position.dodge)) +
      geom_point(cex = 2,position = position_dodge(position.dodge))
  }

  p <- p + xlab("Time (d)") + ylab("Body Weight (g)")

  p <- p + theme_bw() +
    theme(
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background = element_rect(fill = "transparent")
    )   #backgroud
  p <- p + theme(
    axis.title.x = element_text(face = "bold",size = 12),
    axis.text.x  = element_text(vjust = 0,size = 12),
    axis.title.y = element_text(face = "bold",size = 12),
    axis.text.y  = element_text(hjust = 1,size = 12)
  )

  if( pattern == 'oneAN' & level == 'Animal'){
    p <- p + theme(legend.position='none')
    p <- p + facet_wrap( ~ Arms,dir = 'h') +
      theme(strip.text = element_text(colour = "black", face = "bold", size = rel(1)),
            strip.background = element_rect(fill = "#E6AB02", size = rel(1.05), linetype = 1)
      )
  }
  if(pattern == 'TAN'){
    
    if(level == 'Arm'){
      
      Tumor <- s.data[ , 'Tumor']
      
      p <- p + facet_wrap( ~ Tumor,dir = 'h') +
        theme(strip.text = element_text(colour = "black", face = "bold", size = rel(1)),
              strip.background = element_rect(fill = "#E6AB02", size = rel(1.05), linetype = 1)
        )
    }
    if(level == 'Animal'){
      
      Tumor <- data[ , 'Tumor']
      
      p <- p + theme(legend.position='none')
      p <- p + facet_grid(Arms ~ Tumor,space = "free_x", scales = "free_x") +
        theme(strip.text = element_text(colour = "black", face = "bold", size = rel(1)),
              strip.background = element_rect(fill = "#E6AB02",size = rel(1.05), linetype = 1)
        )
    }
  }

  plot(p)

}


####################################################################################################
#' @name plotRC
#' @title plot the relative change of data
#' @description plot the relative change of tumor volume data or body weight data.
#' @usage plotRC(data, type, pattern, rm.baseline = TRUE, orders = NULL)
#'
#' @param data a data frame of measured tumor volume data or body weight data.
#' @param type the type of data, "Vomule"  or "BodyWeight".
#' @param pattern the pattern of PDX trial design, "oneAN" or "TAN".
#' @param orders the redefined order for visulization. For oneAN pattern,the order is redifined orders of the 'Arms'; for TAN pattern, the order is redifined orders of the 'Tumor'.
#' @param rm.baseline a logical parameter, whether remove the baseline.
#' 
#' @import ggplot2
#' 
#' @examples 
#' ### oneAN pattern
#' data(oneAN.volume.data)
#' plotRC(data = oneAN.volume.data,type = 'Volume',pattern = 'oneAN')
#' 
#' data(oneAN.bw.data)
#' plotRC(data = oneAN.bw.data,type = 'BodyWeight',pattern = 'oneAN')
#' 
#' \donttest{
#' ### TAN pattern
#' data(TAN.volume.data)
#' plotRC(data = TAN.volume.data,type = 'Volume',pattern = 'TAN')
#' 
#' data(TAN.bw.data)
#' plotRC(data = TAN.bw.data,type = 'BodyWeight',pattern = 'TAN')
#' }
#' 
#' @export
plotRC<- function(data, type, pattern, rm.baseline = TRUE, orders = NULL){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  #check the order
  if (!is.null(orders)){

    if(pattern == 'oneAN'){
      arms <- unique(data$Arms)
      judged.order <- is.element(orders,arms)

      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Arms' among data.")

      data$Arms <- factor(data$Arms,levels=orders) # set the order
    }

    if(pattern == 'TAN'){
      tumors <- unique(data$Tumor)
      judged.order <- is.element(orders,tumors)

      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Tumor' among data.")

      data$Tumor <- factor(data$Tumor,levels=orders) # set the order
    }
  }

  rc.res <- RelativeChange(data,type,rm.baseline)
  
  Times <- rc.res[ , 'Times']
  RelativeChange <- rc.res[ , 'RelativeChange']
  ID <- rc.res[ , 'ID']
  Arms <- rc.res[ , 'Arms']
  
  p <- ggplot(rc.res, aes(x = Times, y = RelativeChange,group = ID, color = Arms)) +
    geom_line(cex = 0.8) +geom_point(size = 1.5)+

    xlab('Times (d)')+
    ylab(paste('Relative Change of',type, '(%)'))

  if(type == 'BodyWeight'){
    p <- p + geom_hline(aes(yintercept = -20),col = 'blue')+
      geom_hline(aes(yintercept = 0),col = 'black')+
      geom_hline(aes(yintercept = 20),col = 'blue')
  }

  p <- p + theme_bw()+
    theme(
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background  = element_rect(fill = "transparent")
    )   #backgroud
  p <- p + theme(
    axis.title.x = element_text(face = "bold",size = 12),
    axis.text.x  = element_text(vjust = 0,size = 12),
    axis.title.y = element_text(face = "bold",size = 12),
    axis.text.y  = element_text(hjust = 1,size = 12)
  )

  if(pattern == 'oneAN'){
    p <- p + facet_wrap( ~ Arms,dir = 'h') +
      theme(strip.text = element_text(colour = "black", face = "bold", size = rel(1)),
            strip.background = element_rect(fill = "#E6AB02", size = rel(1.05), linetype = 1)
      )
  }
  if(pattern == 'TAN'){
    
    Tumor <- rc.res[ , 'Tumor']
    
    p <- p + facet_grid(Arms ~ Tumor,space = "free_x", scales = "free_x") +
      theme(strip.text = element_text(colour = "black", face = "bold", size = rel(1)),
            strip.background = element_rect(fill = "#E6AB02",size = rel(1.05), linetype  =  1)
      )
  }

  p <- p + theme(legend.position='none')

  plot(p)
}


####################################################################################################
#' @name plotEndpoint
#' @title visulize the tumor volume data measured at the end of experiment
#' @description  visulize the tumor volume data measured at the end of experiment.
#' @usage plotEndpoint(data, pattern, orders = NULL,type = c('barplot','boxplot'))
#' 
#' @param data a data frame of measured tumor volume data.
#' @param pattern the pattern of PDX trial design, "oneAN" or "TAN".
#' @param orders the redefined order for visulization. For oneAN pattern,the order is redifined orders of the 'Arms'; for TAN pattern, the order is redifined orders of the 'Tumor'.
#' @param type the type to visulize, "barplot" or "boxplot".
#'
#' @import ggplot2
#' 
#' @examples
#' ### oneAN pattern
#' data(oneAN.volume.data)
#' plotEndpoint(oneAN.volume.data,type = 'boxplot',pattern = 'oneAN')
#' plotEndpoint(oneAN.volume.data,type = 'barplot',pattern = 'oneAN')
#' 
#' \donttest{
#' ### TAN pattern
#' data(TAN.volume.data)
#' plotEndpoint(TAN.volume.data,pattern = 'TAN',type = 'barplot')
#' plotEndpoint(TAN.volume.data,pattern = 'TAN',type = 'boxplot')
#' }
#' 
#' @export
plotEndpoint <- function(data, pattern, orders = NULL,type = c('barplot','boxplot')){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  endpoint.data <- data[data$Times == max(data$Times),]
  endpoint.data <- subset(endpoint.data,Volume > 0)
  
  if(pattern == 'oneAN') end.data <- DataSummary(data = endpoint.data,type = 'Volume',pattern = 'oneAN',measure.var = 'Volume',group.vars = 'Arms')
  if(pattern == 'TAN') end.data <- DataSummary(data = endpoint.data,type = 'Volume',pattern = 'TAN',measure.var = 'Volume',group.vars = c('Tumor','Arms'))
  
  #check the order
  if (!is.null(orders)){
    
    if(pattern == 'oneAN'){
      arms <- unique(data$Arms)
      judged.order <- is.element(orders,arms)
      
      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Arms' among data.")
      
      endpoint.data$Arms <- factor(endpoint.data$Arms,levels = orders)
      end.data$Arms <- factor(end.data$Arms,levels = orders) # set the order
    }
    
    if(pattern == 'TAN'){
      tumors <- unique(data$Tumor)
      judged.order <- is.element(orders,tumors)
      
      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Tumor' among data.")
      
      endpoint.data$Tumor <- factor(endpoint.data$Tumor,levels = orders)
      end.data$Tumor <- factor(end.data$Tumor,levels = orders) # set the order
    }
  }
  
  if(pattern == 'oneAN'){
    
    Arms <- end.data[ , 'Arms']
    Volume <- end.data[ , 'Volume']
    SE <- end.data[ , 'SE']
    
    if(type == 'barplot'){
      p <- ggplot(end.data,aes(x = Arms,y = Volume,fill = Arms)) + geom_bar(position = "dodge",stat = "identity")
      p <- p+geom_errorbar(aes(ymin = Volume - SE, ymax = Volume + SE,col = Arms),
                           width = 0.15, position = position_dodge(.9))
    }
    
    if(type == 'boxplot'){
      p <- ggplot(endpoint.data,aes(x = Arms,y = Volume,fill = Arms)) + geom_boxplot()
    }
  }
  
  if(pattern == 'TAN'){
    
    Arms <- end.data[ , 'Arms']
    Volume <- end.data[ , 'Volume']
    SE <- end.data[ , 'SE']
    Tumor <- end.data[ , 'Tumor']
    
    if(type == 'barplot'){
      p <- ggplot(end.data,aes(x = Arms,y = Volume,fill = Tumor)) + geom_bar(position = "dodge",stat = "identity")
      p <- p+geom_errorbar(aes(ymin = Volume-SE, ymax = Volume+SE,col = Tumor),
                           width = 0.15, position = position_dodge(.9))
    }
    
    if(type == 'boxplot'){
      p <- ggplot(endpoint.data,aes(x = Arms,y = Volume,fill = Tumor)) + geom_boxplot()
    }
  }
  
  p <- p+xlab("") + ylab(expression(bold(paste("Tumor Volume (",mm^3,")",sep = " "))))
  
  p <- p + theme_bw()+
    theme(
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background  = element_rect(fill = "transparent")
    )   #backgroud
  p <- p+theme(
    axis.text.x  = element_text(hjust = 1,vjust = 1,size = 12,angle  =  45),
    axis.title.y = element_text(face = "bold",size = 12),
    axis.text.y  = element_text(hjust = 1,size = 12)
  )
  
  plot(p)
  
}


####################################################################################################
#' @name plotGR
#' @title visulize the growth rate of tumor volume
#' @description visulize the growth rate of tumor volume.
#' @usage plotGR(data, pattern, orders = NULL, type = c('barplot','boxplot'))
#'
#' @param data a data frame of measured tumor volume data.
#' @param pattern the pattern of PDX trial design, "oneAN" or "TAN".
#' @param orders the redefined order for visulization. For oneAN pattern,the order is redifined orders of the 'Arms';for TAN pattern, the order is redifined orders of the 'Tumor'.
#' @param type the type to visulize, "barplot" or "boxplot".
#' 
#' @import ggplot2
#' 
#' @examples
#' ### oneAN pattern
#' data(oneAN.volume.data)
#' plotGR(oneAN.volume.data,type = 'boxplot',pattern = 'oneAN')
#' plotGR(oneAN.volume.data,type = 'barplot',pattern = 'oneAN')
#'
#' \donttest{
#' ### TAN pattern
#' data(TAN.volume.data)
#' plotGR(TAN.volume.data,pattern = 'TAN',type = 'barplot')
#' plotGR(TAN.volume.data,pattern = 'TAN',type = 'boxplot')
#' }
#' 
#' @export
plotGR <- function(data, pattern, orders = NULL, type = c('barplot','boxplot')){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  data <- data[data$Volume>0,]
  
  ID <- unique(as.character(data$ID))
  
  GR.out <- vector()
  
  for(i in 1:length(ID)){
    sample.i <- data[data$ID == ID[i],]
    
    point.data <- sample.i[sample.i$Times %in% c(min(sample.i$Times),max(sample.i$Times)),]
    g.rate <- (point.data$Volume[2] - point.data$Volume[1])/(point.data$Times[2] - point.data$Times[1])
    
    if(pattern == 'oneAN') {
      g.rate <- data.frame(Arms = sample.i$Arms[1],ID = sample.i$ID[1],growth.rate = g.rate)
      GR.out <- rbind(GR.out,g.rate)
      
      gr.summary <- ddply(.data = GR.out,.variables = 'Arms',
                          .fun = function(xx, col) {
                            c(N    = sum(!is.na(xx[[col]])),
                              Mean = mean  (xx[[col]], na.rm=T),
                              SD   = sd    (xx[[col]], na.rm=T)
                            )},
                          'growth.rate')
    }
    if(pattern == 'TAN') {
      g.rate <- data.frame(Tumor = sample.i$Tumor[1],Arms = sample.i$Arms[1],ID = sample.i$ID[1],growth.rate = g.rate)
      GR.out <- rbind(GR.out,g.rate)
      gr.summary <- ddply(.data = GR.out,.variables = c('Tumor','Arms'),
                          .fun = function(xx, col) {
                            c(N    = sum(!is.na(xx[[col]])),
                              Mean = mean  (xx[[col]], na.rm=T),
                              SD   = sd    (xx[[col]], na.rm=T)
                            )},
                          'growth.rate')
    }
  }
  
  
  gr.summary <- plyr::rename(gr.summary, c("Mean" = 'growth.rate'))
  
  gr.summary$SE <- gr.summary$SD/sqrt(gr.summary$N)  # Calculate standard error
  
  if (!is.null(orders)){
    
    if(pattern == 'oneAN'){
      arms <- unique(data$Arms)
      judged.order <- is.element(orders,arms)
      
      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Arms' among data.")
      
      GR.out$Arms <- factor(GR.out$Arms,levels = orders)
      gr.summary$Arms <- factor(gr.summary$Arms,levels = orders) # set the order
    }
    
    if(pattern == 'TAN'){
      tumors <- unique(data$Tumor)
      judged.order <- is.element(orders,tumors)
      
      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Tumor' among data.")
      
      GR.out$Tumor <- factor(GR.out$Tumor,levels = orders)
      gr.summary$Tumor <- factor(gr.summary$Tumor,levels = orders) # set the order
    }
  }
  
  if(pattern == 'oneAN'){
    
    Arms <- gr.summary[ , 'Arms']
    growth.rate <- gr.summary[ , 'growth.rate']
    SE <- gr.summary[ , 'SE']
    
    if(type == 'barplot'){
      p <- ggplot(gr.summary,aes(x = Arms,y = growth.rate,fill = Arms)) +
        geom_bar(position = "dodge",stat = "identity")
      p <- p + geom_errorbar(aes(ymin = growth.rate-SE, ymax = growth.rate+SE,col = Arms),
                             width = 0.15, position = position_dodge(.9))
    }
    
    if(type == 'boxplot'){
      p <- ggplot(GR.out,aes(x = Arms,y = growth.rate,fill = Arms)) + geom_boxplot()
    }
  }
  
  if(pattern == 'TAN'){
    
    Arms <- gr.summary[ , 'Arms']
    Tumor <- gr.summary[ , 'Tumor']
    growth.rate <- gr.summary[ , 'growth.rate']
    SE <- gr.summary[ , 'SE']
    
    if(type == 'barplot'){
      p <- ggplot(gr.summary,aes(x = Arms,y = growth.rate,fill = Tumor)) + geom_bar(position = "dodge",stat = "identity")
      p <- p + geom_errorbar(aes(ymin = growth.rate-SE, ymax = growth.rate+SE,col = Tumor),
                             width = 0.15, position = position_dodge(.9))
    }
    
    if(type == 'boxplot'){
      p <- ggplot(GR.out,aes(x = Arms,y = growth.rate,fill = Tumor)) + geom_boxplot()
    }
  }
  
  p <- p + xlab("") + ylab(expression(bold(paste("Growth Rate (",mm^3,"/d)",sep = " "))))
  
  p <- p + theme_bw()+
    theme(
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background  = element_rect(fill = "transparent")
    )   #backgroud
  p <- p+theme(
    axis.title.x = element_text(face = "bold",size = 12),
    axis.text.x  = element_text(hjust = 1,vjust = 1,size = 11,angle  = 45),
    axis.title.y = element_text(face = "bold",size = 12),
    axis.text.y  = element_text(hjust = 1,size = 10)
  )
  
  plot(p)
}


####################################################################################################
#' @name TGI
#' @title  Tumor growth inhibition (TGI) rate caculating
#' @description Caculate tumor growth inhibition (TGI) rate of each Arms based on negative control arm.
#' @usage TGI(data, neg.control, pattern = c('oneAN','TAN'), method = c('deltaV','RCV','AUC'))
#'
#' @param data a data frame of measured tumor volume data.
#' @param neg.control the negative control arm.
#' @param pattern the pattern of PDX trial design, "oneAN" or "TAN".
#' @param method the method used to caculate TGI.
#' 
#' @return The TGI value of each experimental arm. 
#' 
#' @note Tumor growth inhibition (TGI) is one of the most widely used indicators to mesure the effect of tumor growth inhibition of treatment. Even though many forms of caculating TGI have been apllied in publications and the field of drug R&D, the basic form is following: 
#'   TGI = (1 - F(VT)/F(VC)). 
#' F(VT) and F(VC) mean the calculating ways for the treatment arm and control arm respectively. DRAP provides three types of F( ) function to calculate TGI: (1) F=Vt-V0; (2) F=Vt/V0; (3) F = area under the curve of tumor volume. Vt and V0 represent the mean tumor volume at the time t and time 0 respectively. 
#' 
#' @examples 
#' ### oenAN pattern
#' data(oneAN.volume.data)
#' oneAN.tgi <- TGI(data = oneAN.volume.data,
#'                  neg.control = 'Control',
#'                  method = 'AUC',
#'                  pattern = 'oneAN')
#' 
#' \donttest{
#' ###TAN pattern
#' data(TAN.volume.data)
#' TAN.tgi <- TGI(TAN.volume.data,
#'                neg.control = 'Control',
#'                pattern = 'TAN',
#'                method = 'AUC')
#' }
#' 
#' @export
TGI <- function(data, neg.control, pattern = c('oneAN','TAN'), method = c('deltaV','RCV','AUC')){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  if(pattern == 'oneAN'){

    #caculate the mean of tumor volume in each Arms at every time point
    in.data <- DataSummary(data,type = 'Volume',pattern = 'OneAN', measure.var = "Volume", group.vars = c("Arms","Times"))

    tgi.out <- vector()

    arms <- unique(as.character(in.data$Arms))
    arms <- setdiff(arms,neg.control)

    control.v <- in.data[in.data$Arms == neg.control,]
    control.v <- control.v[order(control.v$Times),]

    deltaV.control <- control.v$Volume[-1] - control.v$Volume[1]

    RTV.control <- control.v$Volume[-1]/control.v$Volume[1]

    #caculate AUC
    AUC.control <- vector()
    AUC.control.i <- 0
    for(j in 2:nrow(control.v)){
      AUC.control.j <- (control.v$Volume[j-1] + control.v$Volume[j] - 2*control.v$Volume[1])*(control.v$Times[j] - control.v$Times[j-1])/2
      AUC.control.i <- AUC.control.i + AUC.control.j
      AUC.control <- c(AUC.control,AUC.control.i)
    }

    #caculate TGI for each Arms
    for(i in 1:length(arms)){

      arms.i <- in.data[in.data$Arms == arms[i],]
      arms.i <- arms.i[order(arms.i$Times),]

      if (method == 'deltaV'){
        deltaV.i <- arms.i$Volume[-1]-arms.i$Volume[1]
        tgi.i <- 1-deltaV.i/deltaV.control
      }

      if(method == 'RTV'){
        RTV.i <- arms.i$Volume[-1]/arms.i$Volume[1]
        tgi.i <- 1-RTV.i/RTV.control
      }

      if(method == 'AUC'){
        AUC.arm.i <- vector()
        AUC.i <- 0
        for(j in 2:nrow(arms.i)){
          AUC.j <- (arms.i$Volume[j-1]+arms.i$Volume[j]-2*arms.i$Volume[1])*(arms.i$Times[j]-arms.i$Times[j-1])/2
          AUC.i <- AUC.i + AUC.j
          AUC.arm.i <- c(AUC.arm.i,AUC.i)
        }

        tgi.i <- 1-AUC.arm.i/AUC.control
      }

      tgi.ii <- data.frame(Arms=arms[i],Times=arms.i$Times[-1],TGI=tgi.i)

      tgi.out <- rbind(tgi.out,tgi.ii)
    }
  }

  if(pattern == 'TAN'){

    #caculate the mean of tumor volume in each Arms at every time point
    in.data <- DataSummary(data,type = 'Volume',pattern = 'TAN', measure.var = "Volume", group.vars = c("Tumor","Arms","Times"))

    tgi.out <- vector()

    in.data <- in.data[,c('Tumor','Arms','Times','Volume')]

    #caculate TGI for each Tumor in very time
    tumors <- unique(as.character(in.data$Tumor))

    for(i in 1:length(tumors)){

      #tgi.tumor.i <- vector()

      data.tumor.i <- in.data[in.data$Tumor == tumors[i],]

      control.v <- data.tumor.i[data.tumor.i$Arms == neg.control,]
      control.v <- control.v[order(control.v$Times),]

      deltaV.control <- control.v$Volume[-1] - control.v$Volume[1]

      RTV.control <- control.v$Volume[-1]/control.v$Volume[1]

      #caculate AUC
      AUC.control <- vector()
      AUC.control.i <- 0
      for(j in 2:nrow(control.v)){
        AUC.control.j <- (control.v$Volume[j-1] + control.v$Volume[j]-2*control.v$Volume[1])*(control.v$Times[j] - control.v$Times[j-1])/2
        AUC.control.i <- AUC.control.i + AUC.control.j
        AUC.control <- c(AUC.control,AUC.control.i)
      }


      #caculate TGI for each Arms
      arms.tumor.i <- unique(as.character(data.tumor.i$Arms))
      arms.tumor.i <- setdiff(arms.tumor.i,neg.control)

      for(h in 1:length(arms.tumor.i)){

        data.arms.h <- data.tumor.i[data.tumor.i$Arms == arms.tumor.i[h],]
        data.arms.h <- data.arms.h[order(data.arms.h$Times),]

        if (method == 'deltaV'){
          deltaV.h <- data.arms.h$Volume[-1] - data.arms.h$Volume[1]
          tgi.h <- 1 - deltaV.h/deltaV.control
        }

        if(method == 'RTV'){
          RTV.h <- data.arms.h$Volume[-1]/data.arms.h$Volume[1]
          tgi.h <- 1 - RTV.h/RTV.control
        }

        if(method == 'AUC'){

          AUC.arm.h <- vector()
          AUC.h=0
          for(j in 2:nrow(data.arms.h)){
            AUC.k <- (data.arms.h$Volume[j-1] + data.arms.h$Volume[j] - 2*data.arms.h$Volume[1])*(data.arms.h$Times[j] - data.arms.h$Times[j-1])/2
            AUC.h <- AUC.h + AUC.k
            AUC.arm.h <- c(AUC.arm.h,AUC.h)
          }

          tgi.h <- 1 - AUC.arm.h/AUC.control
        }

        tgi.ii <- data.frame(Tumor = tumors[i],Arms = arms.tumor.i[h],Times = data.arms.h$Times[-1],TGI = tgi.h)
        tgi.out <- rbind(tgi.out,tgi.ii)

      }
    }

  }

  return(tgi.out)
}


####################################################################################################
#' @name plotTGI
#' @title TGI visulizing
#' @description Visulize TGI of each Arms.
#' @usage plotTGI(data, pattern, scope = c('end.point','all.point'),orders=NULL)
#' 
#' @param data the output of function TGI.
#' @param pattern the pattern of PDX trial design, "oneAN" or "TAN".
#' @param scope the scope of time point to visulize, whether "end,point" or "all.point".
#' @param orders the redefined order for visulization. For oneAN pattern,the order is redifined orders of the 'Arms'; for TAN pattern, the order is redifined orders of the 'Tumor'.
#' 
#' @import ggplot2
#' 
#' @examples 
#' ### oneAN pattern
#' data(oneAN.volume.data)
#' oneAN.tgi <- TGI(data = oneAN.volume.data,neg.control = 'Control',method = 'AUC',pattern = 'oneAN')
#' plotTGI(oneAN.tgi,pattern = 'oneAN',scope = 'end.point')
#' plotTGI(oneAN.tgi,pattern = 'oneAN',scope = 'all.point')
#' 
#' \donttest{
#' ### TAN pattern
#' data(TAN.volume.data)
#' TAN.tgi <- TGI(TAN.volume.data,  neg.control = 'Control',pattern = 'TAN',method = 'AUC')
#' plotTGI(TAN.tgi,pattern = 'TAN',scope = 'end.point')
#' plotTGI(TAN.tgi,pattern = 'TAN',scope = 'all.point')
#' }
#' 
#' @export
plotTGI <- function(data, pattern, scope = c('end.point','all.point'),orders=NULL){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  #check the order
  if (!is.null(orders)){

    if(pattern == 'oneAN'){
      arms <- unique(data$Arms)
      judged.order <- is.element(orders,arms)

      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Arms' among data.")

      data$Arms <- factor(data$Arms,levels=orders) # set the order
    }

    if(pattern == 'TAN'){
      tumors <- unique(data$Tumor)
      judged.order <- is.element(orders,tumors)

      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Tumor' among data.")

      data$Tumor <- factor(data$Tumor,levels=orders) # set the order
    }
  }

  if(scope == 'end.point'){

    data <- data[data$Times == max(data$Times),]
    
    Arms <- data[ , 'Arms']
    TGI <- data[ , 'TGI']
    
    if(pattern == 'oneAN') {
      p <- ggplot(data, aes(x = Arms, y = TGI, fill = Arms)) + geom_bar(stat = 'identity')
    }
  
    if(pattern == 'TAN') {
      
      Tumor <- data[ , 'Tumor']
      
      p <- ggplot(data, aes(x = Arms, y = TGI, fill = Tumor)) + geom_bar(position = "dodge",stat = 'identity')
    }
    
    p <- p + xlab('') + ylab('TGI')
    p <- p + theme_bw()+
      theme(
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background  = element_rect(fill = "transparent")
      )   #backgroud
    
    p <- p + theme(
      axis.text.x  = element_text(hjust = 1,vjust = 1,size = 12,angle = 45),
      axis.title.y = element_text(face = "bold",size = 12),
      axis.text.y  = element_text(hjust = 1,size = 12)
    )
  }

  if(scope == 'all.point'){
    
    Times <- data[ , 'Times']
    Arms <- data[ , 'Arms']
    TGI <- data[ , 'TGI']
    
    p <- ggplot(data, aes(x = Times, y = TGI,group = Arms, color = Arms)) + geom_line(cex = 1)+geom_point(size = 2)
    p <- p + xlab('Times (d)') + ylab('TGI')
    
    p <- p + theme_bw()+
      theme(
        panel.background = element_rect(fill = "transparent"),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        plot.background  = element_rect(fill = "transparent")
      )   #backgroud
    
    p <- p+theme(
      axis.title.x = element_text(face = "bold",size = 12),
      axis.text.x  = element_text(hjust = 1,vjust = 0.5,size = 12),
      axis.title.y = element_text(face = "bold",size = 12),
      axis.text.y  = element_text(hjust = 1,size = 12)
    )
  }

  if(pattern == 'TAN' & scope == 'all.point'){
    
    Tumor <- data[ , 'Tumor']
    
    p <- p + facet_wrap( ~ Tumor,dir = 'h') +
      theme(strip.text = element_text(colour = "black", face = "bold", size = rel(1)),
            strip.background = element_rect(fill = "#E6AB02", size = rel(1.05), linetype = 1)
      )
  }

  plot(p)
}


####################################################################################################
#' @name DRAnalysis
#' @title  Assess differences in tumor volume between arms
#' @description Assess differences in tumor volume between arms.
#' @usage DRAnalysis (data, pattern = c('oneAN','TAN'),
#'                    method=c('endpoint.ANOVA','endpoint.KW','endpoint.SRH',
#'                             'GR.ANOVA','GR.KW','GR.SRH','mixed.ANOVA','LMM'))
#' 
#' @param data a data frame of measured tumor volume data.
#' @param pattern the pattern of PDX trial design, "oneAN" or "TAN".
#' @param method the method used to analysis.
#' 
#' @import nlme
#' @import rcompanion
#' 
#' @return The statistical analysis results of the choosen method.
#' 
#' @details DRAnalysis offers both parametic and non-parametic statistic methods to access whether the arms have significant effect on tumor volume. users nedd to choose the proper method based on the feature of data and the suitable conditions.
#' 
#' @examples
#' ###oneAN pattern
#' data(oneAN.volume.data)
#' DRAnalysis(oneAN.volume.data,pattern = 'oneAN',method = 'endpoint.ANOVA')
#' DRAnalysis(oneAN.volume.data,pattern = 'oneAN',method = 'endpoint.KW')
#' DRAnalysis(oneAN.volume.data,pattern = 'oneAN',method = 'GR.ANOVA')
#' DRAnalysis(oneAN.volume.data,pattern = 'oneAN',method = 'mixed.ANOVA')
#' DRAnalysis(oneAN.volume.data,pattern = 'oneAN',method = 'LMM')
#' 
#' ###TAN pattern
#' data(TAN.volume.data)
#' DRAnalysis(TAN.volume.data,pattern = 'TAN',method = 'endpoint.ANOVA')
#' DRAnalysis(TAN.volume.data,pattern = 'TAN',method = 'endpoint.SRH')
#' DRAnalysis(TAN.volume.data,pattern = 'TAN',method = 'LMM')
#' 
#' @export
DRAnalysis <- function(data, pattern = c('oneAN','TAN'),
                       method=c('endpoint.ANOVA','endpoint.KW','endpoint.SRH','GR.ANOVA','GR.KW','GR.SRH','mixed.ANOVA','LMM')){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  Volume <- data[,'Volume']

  data <- subset(data,Volume != 0)
  
  if(pattern == 'oneAN'){
    #get endpoint data
    endpoint.data <- data[data$Times == max(data$Times),]
    endpoint.data <- endpoint.data[,c('Arms','Volume')]
    endpoint.data$Arms=factor(endpoint.data$Arms)

    #get mean growth rate
    ID <- unique(as.character(data$ID))

    GR.out <- vector()

    for(i in 1:length(ID)){
      sample.i <- data[data$ID == ID[i],]

      point.data <- sample.i[sample.i$Times %in% c(min(sample.i$Times),max(sample.i$Times)),]
      g.rate <- (point.data$Volume[2] - point.data$Volume[1])/(point.data$Times[2] - point.data$Times[1])
      g.rate <- data.frame(Arms = sample.i$Arms[1],ID = sample.i$ID[1],growth.rate = g.rate)

      GR.out <- rbind(GR.out,g.rate)
      GR.out$growth.rate <- round(GR.out$growth.rate,2)
    }
    GR.out$Arms <- factor(GR.out$Arms)
    #GR.out$ID <- factor(GR.out$ID)

    dra.res <- switch (method,
                       endpoint.ANOVA = summary(aov(Volume ~ Arms,data = endpoint.data)),
                       endpoint.KW    = kruskal.test(Volume ~ Arms,data = endpoint.data),
                       GR.ANOVA       = summary(aov(growth.rate ~ Arms,data = GR.out)),
                       GR.KW          = kruskal.test(growth.rate ~ Arms,data = GR.out),
                       mixed.ANOVA    = summary(aov(Volume ~ Arms + Error(ID/Times), data = data)),
                       LMM            = summary(lme(Volume ~ Arms, random = ~1|ID/Times,data = data))[[20]]
    )
  }

  if(pattern == 'TAN'){
    #get endpoint data
    endpoint.data <- data[data$Times==max(data$Times),]
    endpoint.data <- endpoint.data[,c('Tumor','Arms','Volume')]
    endpoint.data$Arms <- factor(endpoint.data$Arms)

    #get mean growth rate
    ID <- unique(as.character(data$ID))

    GR.out <- vector()

    for(i in 1:length(ID)){
      sample.i <- data[data$ID == ID[i],]

      point.data <- sample.i[sample.i$Times %in% c(min(sample.i$Times),max(sample.i$Times)),]
      g.rate <- (point.data$Volume[2] - point.data$Volume[1])/(point.data$Times[2] - point.data$Times[1])
      g.rate <- data.frame(Tumor = sample.i$Tumor[1],Arms = sample.i$Arms[1],ID = sample.i$ID[1],growth.rate = g.rate)

      GR.out <- rbind(GR.out,g.rate)
      GR.out$growth.rate <- round(GR.out$growth.rate,2)
    }
    GR.out$Arms <- factor(GR.out$Arms)
    GR.out$ID <- factor(GR.out$ID)

    dra.res <- switch (method,
                       endpoint.ANOVA = summary(aov(Volume ~ Arms + Tumor + Arms*Tumor,data = endpoint.data)),
                       endpoint.SRH   = scheirerRayHare(Volume ~ Arms + Tumor,data = endpoint.data),
                       GR.ANOVA       = summary(aov(growth.rate ~ Arms + Tumor + Arms*Tumor,data = GR.out)),
                       GR.SRH         = scheirerRayHare(growth.rate ~ Arms + Tumor,data = GR.out),
                       mixed.ANOVA    = summary(aov(Volume ~ Arms + Tumor+ Error(ID/Times), data = data)),
                       LMM            = summary(lme(Volume ~ Arms,random = ~1|ID/Times,data = data))[[20]]
    )
  }

  return(dra.res)
}



####################################################################################################
#' @name compareGC
#' @title Compare growth curves of tumor volume between Arms
#' @description Compare the growth curves of tumor volume between pairwise arms with permuation test.
#' 
#' @usage 
#' compareGC(data, compare.to = c('all','neg.control','pos.control'), 
#'           control, n = 1000, fun, adjust)
#' 
#' @param data a data frame of measured tumor volume data.
#' @param compare.to the type of compare. "all" means comparing all pairs of arms, "neg.control" means comparing to negative control arm, "pos.control" means comparing to posative control arm.
#' @param control the control arm. While compare to "neg.control", the negative control arm; While compare to "pos.control", the posative control arm.
#' @param n the times of permuation.
#' @param fun funciton used to caculate the statistics for permuation,MeanT or MeanW.
#' @param adjust the method for adjust p.value.
#' 
#' @details compareGC does pairwise comparisons of growth curves between all pairwise arms, or compared to negative/posative arm . p-values can be obtained with permutation test by setting n to some large value, n = 1000.
#' 
#' @return compareGC returns a list with two components, stat and p.value, containing the observed statistics and the estimated p-value. compareGrowthCurves returns a data frame with components
#' Arm1:	  name of first arm in a comparison;
#' Arm2:	  name of second arm in a comparison;
#' Stat:  	observed value of the statistic;
#' p.value:	estimated p-value;
#' q.value:	p-value adjusted for multiple testing.
#' 
#' @references Elso, C.M., et al. Leishmaniasis host response loci (lmr1-3) modify disease severity through a Th1/Th2-independent pathway. Genes and immunity 2004;5(2):93-100.
#' 
#' @seealso  \code{\link{compareTwoGC}}
#' 
#' @examples
#' \donttest{
#' data(oneAN.volume.data)
#' ## compare each pair of arms
#' compareGC(data = oneAN.volume.data, compare.to = 'all', n = 1000, 
#'           fun = MeanT, adjust = "BH")
#' compareGC(data = oneAN.volume.data, compare.to = 'all', n = 1000, 
#'           fun = MeanW, adjust = "BH")
#' 
#' ## compare to negative control
#' compareGC(data = oneAN.volume.data, 
#'           compare.to = 'neg.control', control = 'Control', 
#'           n = 1000, fun = MeanT, adjust = "BH")
#' 
#' ## compare to positive control
#' compareGC(data = oneAN.volume.data, 
#'           compare.to = 'pos.control', control = 'Treatment_1', 
#'           n = 1000, fun = MeanT, adjust = "BH")
#' }
#'   
#' @export
compareGC <- function (data, compare.to = c('all','neg.control','pos.control'), control = NULL,
                       n = 1000, fun, adjust) {
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  data$Arms=as.character(data$Arms)

  volume.summary <- DataSummary(data = data,type = 'Volume',pattern = 'oneAN',measure.var = 'Volume',group.vars = c('Arms','Times'))
  arms <- unique(volume.summary$Arms)

  # check the sample size of each group
  sampleN5 <- volume.summary[volume.summary$N < 5,]

  if(is.element('FALSE', volume.summary$N >= 5)){
    warning('Exist groups that are less than 5 samples.')

    arms <- setdiff(arms,unique(sampleN5$Arms))
  }

  arms.number <- length(arms)

  if (arms.number < 2)
    stop("Less than 2 groups to compare")

  if (compare.to == 'all'){
    arms.1 <- arms.2  <-  rep("", arms.number * (arms.number - 1)/2)
    statistic <- p.value <- q.value <- rep(0, arms.number * (arms.number - 1)/2)
    pair <- 0
    for (i in 1:(arms.number - 1)) {
      for (j in (i + 1):arms.number) {

        pair <- pair + 1
        two.arms  <-  c(arms[i], arms[j])
        two.arms.data <- data[data$Arms %in% two.arms,]
        two.arms.res <- compareTwoGC(two.arms.data,n = n, fun = fun)

        arms.1[pair] <- arms[i]
        arms.2[pair] <- arms[j]
        statistic[pair] <- two.arms.res$stat
        p.value[pair] <- two.arms.res$p.value
      }
    }
  }

  if (compare.to != 'all'){

    other.arms <- setdiff(arms,control)
    arms.number <- length(other.arms)

    arms.1 <- arms.2 <- rep("", arms.number)
    statistic <- p.value <- q.value <- rep(0, arms.number)
    pair <- 0
    for (i in 1:arms.number) {

      pair <- pair + 1
      two.arms <- c(other.arms[i], control)
      two.arms.data <- data[data$Arms %in% two.arms,]
      two.arms.res <- compareTwoGC(two.arms.data,n = n, fun = fun)

      arms.1[pair] <- other.arms[i]
      arms.2[pair] <- control
      statistic[pair] <- two.arms.res$stat
      p.value[pair] <- two.arms.res$p.value
    }
  }

  res <- data.frame(Arm1 = arms.1, Arm2 = arms.2, Stat = statistic, p.value = p.value)
  res$q.value <- p.adjust(p.value, method = adjust)

  return(res)
}


####################################################################################################
#' @name compareTwoGC
#' @title Compare growth curves of tumor volume between pairwise Arms
#' @description compare the growth curve of tumor volume between pairwise arms with permuation test.
#' @usage compareTwoGC(data, n = n, fun)
#'
#' @param data a data frame of measured tumor volume data.
#' @param n the times of permuation.
#' @param fun funciton used to caculate the statistics for permuation.
#' 
#' @references Elso, C.M., et al. Leishmaniasis host response loci (lmr1-3) modify disease severity through a Th1/Th2-independent pathway. Genes and immunity 2004;5(2):93-100.
#' 
#' @seealso  \code{\link{compareGC}}
#' 
#' @export
compareTwoGC <- function (data, n = n, fun) {
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  arms <- as.character(data$Arms)
  arms <- unique(arms)
  if (length(arms) != 2)
    stop("Must be exactly 2 Arms")

  stat.obs <- fun(data[data$Arms == arms[1],], data[data$Arms == arms[2],])
  asbig  <-  0

  for (i in 1:n) {
    id <- unique(data$ID)
    id.1 <- length(unique(data[data$Arms == arms[1],'ID']))
    id.2 <- length(unique(data[data$Arms == arms[2],'ID']))

    id <- sample(id)
    data.1 <- data[data$ID %in% id[1:id.1],]
    data.2 <- data[data$ID %in% id[(id.1+1):(id.1+id.2)],]
    statistic <- fun(data.1,data.2)

    if (abs(statistic) >= abs(stat.obs))
      asbig <- asbig + 1
  }
  list(stat = stat.obs, p.value = asbig/n)
}

###################################################################################################
#' @name MeanT
#' @title Caculating mean t-statistic between two arms of tumor volume growth curves
#' @description Caculating mean t-statistic between two arms of tumor volume growth curves for all time points.
#' 
#' @param data.1 a data frame of tumor volume for the first arm.
#' @param data.2 a data frame of tumor volume for the second arm.
#'
#' @details This function computes the t-statistic of tumor volumes between two arms at each time point, and returns the mean of t-statistics at all time points. This function is used by compareGC.
#' @return mean t-statistic.
#' 
#' @references Elso, C.M., et al. Leishmaniasis host response loci (lmr1-3) modify disease severity through a Th1/Th2-independent pathway. Genes and immunity 2004;5(2):93-100.
#' 
#' @seealso MeaW
#' 
#' @export
MeanT <- function (data.1, data.2) {
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  if (is.null(dim(data.1)) || is.null(dim(data.2)))
    return(NA)

  data.1.times <- unique(as.vector(data.1$Times))
  data.2.times <- unique(as.vector(data.2$Times))

  if (!setequal(data.1.times,data.2.times))
    stop("Number of time points must match")

  times.point <- data.1.times
  m1 <- m2 <- t.statistic <- rep(0,length(times.point))

  for (t in 1: length(times.point)){
    d.1 <- data.1[data.1$Times == times.point[t],'Volume']
    d.2 <- data.2[data.2$Times == times.point[t],'Volume']
   
    t.statistic[t] <- t.test(d.1,d.2)$statistic
  }
  statistic <- mean(t.statistic,na.rm = TRUE)

  return(statistic)
}


###################################################################################################
#' @name MeanW
#' @title Caculating mean Wilcox-statistic between two arms of tumor volume growth curves
#' @description Caculating mean Wilcox-statistic between two arms of tumor volume growth curves for all time points.
#'
#' @param data.1 a data frame of tumor volume for the first arm.
#' @param data.2 a data frame of tumor volume for the second arm.
#' 
#' @details This function computes the W-statistic of tumor volumes between two arms at each time point, and returns the mean of W-statistics at all time points. This function is used by compareGC.
#' @return mean W-statistic.
#' 
#' @seealso MeaT 
#' @export
MeanW <- function (data.1, data.2) {
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  if (is.null(dim(data.1)) || is.null(dim(data.2)))
    return(NA)
  
  data.1.times <- unique(as.vector(data.1$Times))
  data.2.times <- unique(as.vector(data.2$Times))
  
  if (!setequal(data.1.times,data.2.times))
    stop("Number of time points must match")
  
  times.point <- data.1.times
  m1 <- m2 <- w.statistic <- rep(0,length(times.point))
  
  for (t in 1: length(times.point)){
    d.1 <- data.1[data.1$Times == times.point[t],'Volume']
    d.2 <- data.2[data.2$Times == times.point[t],'Volume']
   
    w.statistic[t] <- wilcox.test(d.1,d.2)$statistic # kruskal.test
  }
  statistic <- mean(w.statistic,na.rm = TRUE)
  
  return(statistic)
}


####################################################################################################
#' @name DRLevel
#' @title Define the drug response level for each animal.
#' @description Define the drug response level of each animal based on tumor volume change.
#' @usage DRLevel(data, method = c('NPDXE.Response','PPTP.Response','RC.Response'), 
#'                criteria, neg.control, rm.neg.control=TRUE)
#' 
#' @param data a data frame of measured volume data.
#' @param method the method used to quantify the drug response level. Currently available methods include NPDXE.Response (default),PPTP.Response,RC.Response.
#' @param criteria the criteria conrresponding to method.
#' @param neg.control the negative control arm.
#' @param rm.neg.control whether remove the negative control arm.
#' 
#' @details Defining drug response level is the general pipeline in clinical trial and precinical animal trial. DRLevel offers three published ways to define response level in function DRLevel. Notably, the criteria of each way for defining response level could be adjusted by users based on the actual experimental data.
#' 
#' @return The response level of each animal.
#' @import dplyr
#' 
#' @references 
#' Gao, H., et al. High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response. Nat Med 2015;21(11):1318-1325.
#' Murphy, B., et al. Evaluation of Alternative In Vivo Drug Screening Methodology: A Single Mouse Analysis. Cancer Res 2016;76(19):5798-5809.
#' Bertotti, A., et al. The genomic landscape of response to EGFR blockade in colorectal cancer. Nature 2015;526(7572):263-267.
#' 
#' @seealso \code{\link{NPDXEResponseLevel}, \link{PPTPResponseLevel}, \link{RCResponseLevel}}
#' 
#' @examples
#' ### build the criteria
#' 
#' # NPDXE.criteria
#' npdxe.criteria <- data.frame(BestResponse.lower = c(-1000,-0.95,-0.5,0.35),
#'                              BestResponse.upper = c(-0.95,-0.5,0.35,1000),
#'                              BestAvgResponse.lower = c(-1000,-0.4,-0.2,0.3), 
#'                              BestAvgResponse.upper = c(-0.4,-0.2,0.3,1000), 
#'                              Level = c( 'CR','PR', 'SD','PD'))
#' npdxe.criteria
#'  
#' # PPTP.criteria 
#' pptp.criteria <- data.frame(min.RC.lower = c(-1,-1,-0.5,-0.5),
#'                             min.RC.upper = c(-0.5,-0.5,1000,1000),
#'                             min.Vol.lower = c(0,100,NA,NA),
#'                             min.Vol.upper = c(100,10000,NA,NA),
#'                             end.RC.lower =c(NA,NA,-1,0.25),
#'                             end.RC.upper = c(NA,NA,0.25,1000),
#'                             Level = c( 'CR','PR', 'SD','PD'))
#' pptp.criteria
#' 
#' # RC.criteria
#' rc.criteria <- data.frame(Response.lower = c(-1000,-0.35,0.35),
#'                           Response.upper = c(-0.35,0.35,1000),
#'                           Level = c( 'CR','SD','PD'))
#' rc.criteria
#' 
#' ### drug response level
#' ##
#' data(oneAN.volume.data)
#' oneAN.drl <- DRLevel(data = oneAN.volume.data, 
#'                      method = 'NPDXE.Response', 
#'                      criteria = npdxe.criteria, 
#'                      neg.control = 'Control')
#' oneAN.drl <- oneAN.drl[order(oneAN.drl$Arms),]
#' head(oneAN.drl)
#' 
#' \donttest{
#' ##
#' data(TAN.volume.data)
#' TAN.drl <- DRLevel(data = TAN.volume.data, 
#'                    method = 'NPDXE.Response', 
#'                    criteria = npdxe.criteria, 
#'                    neg.control = 'Control')
#' head(TAN.drl)
#' 
#' ##
#' data(TAone.volume.data)
#' head(TAone.volume.data)
#' TAone.drl <- DRLevel(data = TAone.volume.data, 
#'                      method = 'NPDXE.Response', 
#'                      criteria = npdxe.criteria, 
#'                      neg.control = 'Control')
#' head(TAone.drl)
#' }
#' 
#' @export
DRLevel <- function(data, method = c('NPDXE.Response','PPTP.Response','RC.Response'), criteria, neg.control, rm.neg.control=TRUE){
  
  if(method == 'NPDXE.Response'){

    warning("The 'Reference' for suitable condition of NPDXE response criteria is: The initial tumor volume is about 200 mm3.")
  
  }
  
  if(method == 'PPTP.Response'){

    warning("The 'Reference' for suitable condition of PPTP response criteria includes: 
            1. The initial tumor volume is between 200 and 500 mm3; 
            2. If the tumor volume readches to 4-fold of its initial volume, the mouse reaches endpoint; else the drug treatment lasts to 42 days.")
  }
  
  if(method == 'RC.Response'){

    warning("The 'Reference' for suitable condition of RC response criteria includes: 
            1. The initial tumor volume is about 400 mm3;
            2. Drug treatment lasts to 21 days.")
  }
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  Volume <- data[,'Volume']
  indata <- subset(data,Volume > 0)
  
  
  drug.response.level <- switch (method,
                                 NPDXE.Response = NPDXEResponseLevel(indata, criteria),
                                 PPTP.Response  = PPTPResponseLevel(indata, criteria),
                                 RC.Response    = RCResponseLevel(indata, criteria)
  )

  indata <- indata[,setdiff(colnames(indata),c('Times','Volume'))]
  indata <- unique(indata)

  drug.response.res <- merge(indata,drug.response.level,by = 'ID')

  if(rm.neg.control){
    drug.response.res <- drug.response.res[as.character(drug.response.res$Arms) %in% setdiff(as.character(drug.response.res$Arms), neg.control), ]
  }

  return(drug.response.res)

}


####################################################################################################
#' @name RCResponseLevel
#' @title Define response level based on relatvie change of tumor volume
#' @description Define the drug response level of each animal based on the relative change of tumor volume.
#'
#' @param data a data frame of measured volume data.
#' @param criteria the criteria used to defined response level.
#' 
#' @details The 'Reference' for suitable condition of RC response criteria includes: 1. The initional tumor volume is about 400 mm3; 2. Drug treatment lasts to 21 days.
#' @return The response level of each animal.
#' @import dplyr
#' 
#' @references 
#' Bertotti, A., et al. The genomic landscape of response to EGFR blockade in colorectal cancer. Nature 2015;526(7572):263-267.
#' 
#' @seealso \code{\link{DRLevel}, \link{NPDXEResponseLevel}, \link{PPTPResponseLevel}}
#' 
#' @export
RCResponseLevel <- function(data, criteria){

  #library(dplyr)
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  data <- data[data$ID %in% setdiff(unique(data$ID),unique(data[data$Volume==0,'ID'])),]

  Response.Level <- vector()
  ID <- unique(data$ID)

  for ( i in 1:length(ID)){
    data.id.i <- data[data$ID == ID[i],]

    data.id.i <- data.id.i[order(data.id.i$Times),]
    volume.change.rate <- (data.id.i$Volume[nrow(data.id.i)] - data.id.i$Volume[1])/data.id.i$Volume[1]

    for (j in 1:nrow(criteria)) {
      if(between(volume.change.rate,criteria[j,'Response.lower'],criteria[j,'Response.upper'])){
        response.level <- as.character(criteria$Level)[j]
      }
    }

    response.level.id.i <- data.frame(ID = ID[i],Volume.CR = volume.change.rate,
                                      Response.Level = response.level)

    Response.Level <- rbind(Response.Level,response.level.id.i)
  }

  return(Response.Level)
}


####################################################################################################
#' @name PPTPResponseLevel
#' @title Define response level based on PPTP response criteria
#' @description Define the drug response level of each animal based on the PPTP response criteria.
#'
#' @param data a data frame of measured volume data.
#' @param criteria the criteria used to defined response level.
#' 
#' @details The 'Reference' for suitable condition of PPTP response criteria includes: 1. The initional tumor volume is between 200 and 500 mm3; 2. If the tumor volume readches to 4-fold of its initial volume, the mouse reaches endpoint; else the drug treatment lasts to 42 days.
#' 
#' @return The response level of each animal.
#' @import dplyr
#' 
#' @references 
#' Murphy, B., et al. Evaluation of Alternative In Vivo Drug Screening Methodology: A Single Mouse Analysis. Cancer Res 2016;76(19):5798-5809.
#' 
#' @seealso \code{\link{DRLevel}, \link{NPDXEResponseLevel}, \link{RCResponseLevel}}
#' 
#' @export
PPTPResponseLevel <- function(data, criteria){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  Response.Level <- vector()
  ID <- unique(data$ID)
  
  for ( i in 1:length(ID)){
    data.id.i <- data[data$ID == ID[i],]
    
    data.id.i <- data.id.i[order(data.id.i$Times),]
    volume.RC <- (data.id.i$Volume - data.id.i$Volume[1])/data.id.i$Volume[1]
    
    volume.RC <- cbind(data.id.i[,c('ID','Times','Volume')],volume.RC)
    volume.RC <- volume.RC[-1,]

    min.RC <- min(volume.RC$volume.RC)
    min.Vol <- min(volume.RC$Volume)
    end.RC <- volume.RC$volume.RC[which.max(volume.RC$Times)]
    
    min.RC.criteria <- unique(c(criteria[criteria$Level %in% c('PD','SD'),'min.RC.lower'],criteria[criteria$Level %in% c('CR','PR'),'min.RC.upper']))
    min.Vol.criteria <- unique(c(criteria[criteria$Level == 'PR','min.Vol.lower'], criteria[criteria$Level == 'CR','min.Vol.upper']))
    end.RC.criteria <- unique(c(criteria[criteria$Level == 'PD','end.RC.lower'],criteria[criteria$Level == 'SD','end.RC.upper']))
    
    if(length(min.RC.criteria) != 1 & length(min.Vol.criteria) != 1 & length(end.RC.criteria) != 1){
      stop("The PPTP response criteria you set exits bug. Please check your inputing criteria")
    }
    
    if(min.RC > min.RC.criteria){
      if(end.RC > end.RC.criteria) {response.level <- 'PD'}
      if(end.RC < end.RC.criteria) {response.level <- 'SD'}
    }
    if(min.RC < min.RC.criteria){
      if(min.Vol > min.Vol.criteria) response.level <- 'PR'
      if(min.Vol < min.Vol.criteria) response.level <- 'PR'
    }
  
    response.level.id.i <- data.frame(ID = ID[i], min.RC = min.RC, min.Vol = min.Vol, end.RC = end.RC, 
                                     Response.Level = response.level)
    
    
    Response.Level <- rbind(Response.Level,response.level.id.i)
  }
  
  return(Response.Level)
}


####################################################################################################
#' @name NPDXEResponseLevel
#' @title Define drug response level based on Novartis PDXE response criteria
#' @description Define the drug response level of each animal based on Novartis PDXE response criteria.
#'
#' @param data data frame of measured volume data.
#' @param criteria the criteria used to defined response level.
#' 
#' @details The 'Reference' for suitable condition of NPDXE response cirteria is: The initional tumor volume is about 200 mm3.
#' 
#' @return The response level of each animal.
#' @import dplyr
#' 
#' @references 
#' Gao, H., et al. High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response. Nat Med 2015;21(11):1318-1325.
#' 
#' @seealso \code{\link{DRLevel}, \link{PPTPResponseLevel}, \link{RCResponseLevel}}
#' @export
NPDXEResponseLevel <- function(data, criteria){
  
  #library(dplyr)
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  level.order <- data.frame(Level  = c('CR','PR', 'SD','PD'),
                            orders = c(1:4),stringsAsFactors = F)
  
  best.response.criteria <- criteria[,c("BestResponse.lower","BestResponse.upper","Level")]
  best.avg.response.criteria <- criteria[,c("BestAvgResponse.lower","BestAvgResponse.upper","Level")]
  
  Response.Level <- vector()
  
  ID <- unique(data$ID)
  
  for ( i in 1:length(ID)){
    data.id.i <- data[data$ID == ID[i],]
    t.max <- max(data.id.i$Times)
    
    if(t.max < 10) next  # next ID
    
    data.id.i <- data.id.i[order(data.id.i$Times),]
    volume.change.rate <- (data.id.i$Volume - data.id.i$Volume[1])/data.id.i$Volume[1]
    
    vc.rate <- cbind(data.id.i[,c('ID','Times')],volume.change.rate)
    
    ##best response
    Times <- vc.rate[,'Times']
    rate_best <- subset(vc.rate,Times >= 10)
    best.response <- min(rate_best$volume.change.rate)
    
    
    avg.response <- vector()
    for(j in 1:nrow(vc.rate)){
      avg.response.j <- mean(volume.change.rate[1:j]) #calulate the mean of volume relative change.
      avg.response <- c(avg.response,avg.response.j)
    }
    avg.response <- cbind(data.id.i[,c('ID','Times')],avg.response)
    
    ##best average response
    best.avg.response <- min(avg.response[avg.response$Times >= 10,'avg.response'])
    
    
    for (j in 1:nrow(best.response.criteria)) {
      if(between(best.response,best.response.criteria[j,'BestResponse.lower'],best.response.criteria[j,'BestResponse.upper'])){
        best.response.level <- as.character(criteria$Level)[j]
      }
    }
    
    for (j in 1:nrow(best.avg.response.criteria)) {
      if(between(best.avg.response,best.avg.response.criteria[j,'BestAvgResponse.lower'],best.avg.response.criteria[j,'BestAvgResponse.upper'])){
        best.avg.response.level <- as.character(criteria$Level)[j]
      }
    }
    
    if(best.response != best.avg.response){
      order.in.level <- level.order[level.order$Level %in% c(best.response.level,best.avg.response.level),]
      response.level.i <- order.in.level[which.max(order.in.level$orders),'Level']  #choose the lower level.
    }else{
      response.level.i <- best.response.level
    }
    
    response.level.id.i <- data.frame(ID = ID[i],
                                      Best.Response = best.response, 
                                      Best.Avg.Response = best.avg.response,
                                      Response.Level = response.level.i,
                                      stringsAsFactors = F)
    
    Response.Level <- rbind(Response.Level,response.level.id.i)
  }
  
  return(Response.Level)
}


####################################################################################################
#' @name DRLevelSummary
#' @title Summary the drug response level
#' @description Summary the frequency and propotion of animals at every drug response level in each group.
#' @usage 
#' DRLevelSummary(data, by)
#' 
#' @param data drug response level, the output of DRLevel.
#' @param by the group summaried by.
#' 
#' @examples
#' # NPDXE.criteria
#' npdxe.criteria <- data.frame(BestResponse.lower = c(-1000,-0.95,-0.5,0.35),
#'                              BestResponse.upper = c(-0.95,-0.5,0.35,1000),
#'                              BestAvgResponse.lower = c(-1000,-0.4,-0.2,0.3), 
#'                              BestAvgResponse.upper = c(-0.4,-0.2,0.3,1000), 
#'                              Level = c( 'CR','PR', 'SD','PD'))
#' npdxe.criteria
#'                             
#' ### oneAN pattern
#' data(oneAN.volume.data)
#' oneAN.drl <- DRLevel(data = oneAN.volume.data, 
#'                      method = 'NPDXE.Response', 
#'                      criteria = npdxe.criteria, 
#'                      neg.control = 'Control')
#' DRLevelSummary(oneAN.drl,by = 'Arms')
#' 
#' \donttest{
#' ### TAN pattern
#' data(TAN.volume.data)
#' TAN.drl <- DRLevel(data = TAN.volume.data, 
#'                    method = 'NPDXE.Response', 
#'                    criteria = npdxe.criteria, 
#'                    neg.control = 'Control')
#' DRLevelSummary(TAN.drl,by = 'Arms')
#' DRLevelSummary(TAN.drl,by = c('Arms','Tumor'))
#' 
#' ### TAone pattern
#' data(TAone.volume.data)
#' TAone.drl <- DRLevel(data = TAone.volume.data, 
#'                      method = 'NPDXE.Response', 
#'                      criteria = npdxe.criteria, 
#'                      neg.control = 'Control')
#' head(TAone.drl)
#' DRLevelSummary(data = TAone.drl,by=c('Arms','Type'))
#' }
#' 
#' @export
DRLevelSummary <- function(data, by){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  colnames(data)[grep('Level',colnames(data),fixed = T)] = 'Level'

  drl.summary <- vector()

  dr.group <- unique(data[,by])
  dr.group <- data.frame(dr.group,stringsAsFactors = F)
  colnames(dr.group) <- by

  for(i in 1:nrow(dr.group)){
    if(length(by) == 1) {
      dr.i <- data[data[,by] == dr.group[i,by],]
    }

    if(length(by) > 1) {
      dr.i <- merge(dr.group[i,],data,by)
    }

    numb <- nrow(dr.i)
    freq <- as.data.frame(table(dr.i$Level),stringsAsFactors = F)
    colnames(freq) <- c('Level','Frequency')
    prop <- prop.table(freq$Frequency)
    
    levels <- unique(as.character(dr.i$Level))
    
    if(length(by) > 1){
      dr.group.i <- vector()
      for(j in 1:length(levels)){
        dr.group.i <- rbind(dr.group.i,dr.group[i,])
      }

      drl.i <- data.frame(dr.group.i,Number = numb,freq, Proportion = prop, stringsAsFactors = F)
    }

    if(length(by) == 1){

      drl.i <- data.frame(dr.group[i,],Number = numb,freq,Proportion = prop,stringsAsFactors = F)
      colnames(drl.i)[1] <- by
    }

    drl.summary <- rbind(drl.summary,drl.i)

  }

  drl.summary <- drl.summary[do.call(order,as.data.frame(drl.summary)),]

  drl.summary$Level <- factor(drl.summary$Level,levels = unique(data$Level))

  rownames(drl.summary) <- NULL

  return(drl.summary)
}


####################################################################################################
#' @name plotDRLevel
#' @title Visulize the drug response level
#' @description Visulize the drug response level for each arm.
#' @usage 
#' plotDRLevel(data, by = c("Arms", "Tumor"), 
#'            pattern = c("oneAN",'TAN'), orders=NULL, 
#'            orders.1=NULL, orders.2=NULL)
#' 
#' @param data drug response level ,the output of DRLevel.
#' @param by the level summaried by.
#' @param pattern the pattern of PDX trial design, "oneAN" or "TAN".
#' @param orders the redefined order while the argument "by" just contain one variable. 
#' @param orders.1 the redefined order of the first varible while the argument "by" contain two variables.
#' @param orders.2 the redefined order of the second varible while the argument "by" contain two variables.
#' 
#' @import ggplot2
#' 
#' @examples
#' # NPDXE.criteria
#' npdxe.criteria <- data.frame(BestResponse.lower = c(-1000,-0.95,-0.5,0.35),
#'                              BestResponse.upper = c(-0.95,-0.5,0.35,1000),
#'                              BestAvgResponse.lower = c(-1000,-0.4,-0.2,0.3), 
#'                              BestAvgResponse.upper = c(-0.4,-0.2,0.3,1000), 
#'                              Level = c( 'CR','PR', 'SD','PD'))
#' npdxe.criteria
#'                             
#' ### oneAN pattern
#' data(oneAN.volume.data)
#' oneAN.drl <- DRLevel(data = oneAN.volume.data, 
#'                      method = 'NPDXE.Response', 
#'                      criteria = npdxe.criteria, 
#'                      neg.control = 'Control')
#' DRLevelSummary(oneAN.drl,by = 'Arms')
#' plotDRLevel(data = oneAN.drl,by = 'Arms',pattern = 'oneAN')
#' 
#' \donttest{
#' ### TAN pattern
#' data(TAN.volume.data)
#' TAN.drl <- DRLevel(data = TAN.volume.data, 
#'                    method = 'NPDXE.Response', 
#'                    criteria = npdxe.criteria, 
#'                    neg.control = 'Control')
#' DRLevelSummary(TAN.drl,by = 'Arms')
#' plotDRLevel(data = TAN.drl, by='Arms', pattern = 'TAN')
#' plotDRLevel(data = TAN.drl, by=c('Arms','Tumor'), pattern = 'TAN')
#' 
#' # rectify orders
#' arms <- unique(TAN.volume.data$Arms)
#' arms <- arms[-1]
#' arms
#' tumors <- unique(TAN.volume.data$Tumor)
#' tumors
#' plotDRLevel(data = TAN.drl, by = c('Arms','Tumor'), 
#'             pattern = 'TAN', orders.1 = arms, orders.2 = tumors)
#' }
#' 
#' @export
plotDRLevel <- function(data, by = c("Arms", "Tumor"), pattern = c("oneAN",'TAN'), orders=NULL, orders.1=NULL, orders.2=NULL){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  drl.summary <- DRLevelSummary(data,by)
  drl.summary <- drl.summary[drl.summary$Frequency > 0, ]

  #check the order
  if (!is.null(orders) | !is.null(orders.1) | !is.null(orders.2)){

    if(is.element('Arms',by) & length(by) == 1){
      arms <- unique(data$Arms)
      judged.order <- is.element(orders,arms)

      if ("FALSE" %in% judged.order)
        stop("The input order is improper. Please ensure the input order is the order of 'Arms' among data.")

      drl.summary$Arms <- factor(drl.summary$Arms,levels=orders) # set the order
    }

    if(is.element('Arms',by) & is.element('Tumor', by) & length(by) == 2){
      arms <- unique(data$Arms)
      judged.order.1 <- is.element(orders.1,arms)
      
      if ("FALSE" %in% judged.order.1)
        stop("The input order is improper. Please ensure the input order is the order of 'Arms' among data.")
      
      drl.summary$Arms <- factor(drl.summary$Arms,levels=orders.1) # set the order
      
      tumors <- unique(data$Tumor)
      judged.order.2 <- is.element(orders.2,tumors)

      if ("FALSE" %in% judged.order.2)
        stop("The input order is improper. Please ensure the input order is the order of 'Tumor' among data.")

      drl.summary$Tumor <- factor(drl.summary$Tumor,levels=orders.2) # set the order
    }
  }

  if(is.element('Arms',by) & length(by) == 1){
    
    Arms <- drl.summary[ , 'Arms']
    Proportion <- drl.summary[ , 'Proportion']
    Level <- drl.summary[ , 'Level']
    
    p <- ggplot(drl.summary,aes(x = Arms,y = Proportion,fill = Level)) +
      geom_bar(stat = 'identity', position = 'stack')
    
  } else{
    
    Arms <- drl.summary[ , 'Arms']
    Tumor <- drl.summary[ , 'Tumor']
    Proportion <- drl.summary[ , 'Proportion']
    Level <- drl.summary[ , 'Level']
    
    p <- ggplot(drl.summary,aes(x = Tumor,y = Proportion,fill = Level)) +
      geom_bar(stat = 'identity', position = 'stack')
  }

  p <- p + xlab("")+ylab('Propotion')

  p <- p + geom_text(mapping = aes(label = drl.summary$Frequency),colour = "black",
                     position = position_stack(0.5),vjust = 0,size = 5)

  p <- p + theme_bw()+
       theme(
         panel.background = element_rect(fill = "transparent"),
         panel.grid.minor = element_blank(),
         panel.grid.major = element_blank(),
         plot.background = element_rect(fill = "transparent")
       )   #backgroud
  p <- p+theme(
    axis.title.x = element_text(face ="bold", size = 12),
    axis.text.x  = element_text(hjust = 1, vjust = 1, size = 12, angle = 45),
    axis.title.y = element_text(face = "bold", size = 12),
    axis.text.y  = element_text(hjust = 1,size = 12)
    )

  if(pattern == 'TAN'){
    p <- p + facet_grid(.~Arms,scales = 'free',space = 'free_x')+
         theme(strip.text =element_text(colour = "black", face ="bold", size =rel(1)),
               strip.background =element_rect(fill ="#E6AB02",size =rel(1.05), linetype = 1)
      )
  }

  plot(p)
}


####################################################################################################
#' @name DRLevelAnalysis
#' @title Caculate response evaluation indexes
#' @description Caculate response evaluation indexes, such as RR and DCR, based on drug response level.
#' @usage DRLevelAnalysis(data, by, measurement = c('RR','DCR','both'))
#' 
#' @param data drug response level ,the output of DRLevel.
#' @param by the level summaried by.
#' @param measurement the measurement of drug response criteria calcuated, whether "RR", "DCR", or "both".
#' 
#' @details RR (response rate) is the percentage of animals that response level is better than partial repsone (PR) after treatment. DCR (decease control rate) is the percentage of animals that response level is better than stable disease (SD) after treatment.
#' 
#' @examples
#' # NPDXE.criteria 
#' npdxe.criteria <- data.frame(BestResponse.lower = c(-1000,-0.95,-0.5,0.35),
#'                              BestResponse.upper = c(-0.95,-0.5,0.35,1000),
#'                              BestAvgResponse.lower = c(-1000,-0.4,-0.2,0.3), 
#'                              BestAvgResponse.upper = c(-0.4,-0.2,0.3,1000), 
#'                              Level = c( 'CR','PR', 'SD','PD'))
#' npdxe.criteria
#' 
#' ### oneAN pattern
#' data(oneAN.volume.data)
#' oneAN.drl <- DRLevel(data = oneAN.volume.data, method = 'NPDXE.Response', 
#'                      criteria = npdxe.criteria, neg.control = 'Control')
#' DRLevelAnalysis(oneAN.drl,by = 'Arms', measurement = 'both')
#' 
#' \donttest{
#' ### TAN pattern
#' data(TAN.volume.data)
#' TAN.drl <- DRLevel(data = TAN.volume.data, method = 'NPDXE.Response', 
#'                    criteria = npdxe.criteria, neg.control = 'Control')
#' DRLevelAnalysis(TAN.drl,by='Arms',measurement = 'both')
#' DRLevelAnalysis(TAN.drl,by=c('Arms','Tumor'),measurement = 'both')
#' 
#' ### TAone pattern
#' data(TAone.volume.data)
#' TAone.drl <- DRLevel(data = TAone.volume.data, method = 'NPDXE.Response', 
#'                      criteria = npdxe.criteria, neg.control = 'Control')
#' head(TAone.drl)
#' TAone.drl.analysis = DRLevelAnalysis(data = TAone.drl,by = c('Arms','Type'),measurement = 'both')
#' head(TAone.drl.analysis)
#' }
#' 
#' @export
DRLevelAnalysis <-function(data, by, measurement = c('RR','DCR','both')){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  drl.summary <-DRLevelSummary(data,by)

  dr.group <-unique(drl.summary[,by])
  dr.group <-data.frame(dr.group,stringsAsFactors = F)
  colnames(dr.group) <-by

  drl.analysis <-vector()

  for(i in 1:nrow(dr.group)){
    if(length(by) == 1) {
      dr.i <-drl.summary[drl.summary[,by] == dr.group[i,by],]
    }

    if(length(by) > 1) {
      dr.i <-merge(dr.group[i,],drl.summary,by)
    }

    upper.1 <-sum(dr.i[dr.i$Level %in% c('CR','PR'),'Frequency'])
    upper.2 <-sum(dr.i[dr.i$Level %in% c('CR','PR','SD'),'Frequency'])
    total <-sum(dr.i$Frequency)

    response.rate <-upper.1/total    #caculate response rate

    desease.control.rate <-upper.2/total  #caculate desease control rate

    if(length(by) > 1){
      dr.group.i <-vector()
      for(j in 1:length(levels)){
        dr.group.i <-rbind(dr.group.i,dr.group[i,])
      }

      drl.i <-data.frame(dr.group.i,Number = total,RR=response.rate,DCR=desease.control.rate)
    }

    if(length(by) == 1){

      drl.i <-data.frame(dr.group[i,],Number = total,RR = response.rate,DCR = desease.control.rate)
      colnames(drl.i)[1] <-by
    }

    drl.analysis <-rbind(drl.analysis,drl.i)

  }

  drl.analysis <-drl.analysis[do.call(order,as.data.frame(drl.analysis)),]
  rownames(drl.analysis) <-NULL

  if(measurement == 'both') drl.analysis.res = drl.analysis
  if(measurement == 'RR')   drl.analysis.res = drl.analysis[,c(by,'Number','RR')]
  if(measurement == 'DCR')  drl.analysis.res = drl.analysis[,c(by,'Number','DCR')]

  return(drl.analysis.res)
}


####################################################################################################
#' @name plotDRResults
#' @title plot the results of drug response level analysis.
#' @description Visulize the result of drug response analysis, such as RR or DCR of each arm.
#' @usage plotDRResults(data, by, min.number, measurement = c('RR','DCR','both'))
#' 
#' @param data drug response level ,the output of DRLevel.
#' @param by the level summaried by.
#' @param min.number the lower cutoff of sample number of every arm
#' @param measurement the response evaluation indexes used,  whether "RR", "DCR", or "both".
#' 
#' @import ggplot2
#' 
#' @examples
#' # NPDXE.criteria 
#' npdxe.criteria <- data.frame(BestResponse.lower = c(-1000,-0.95,-0.5,0.35),
#'                              BestResponse.upper = c(-0.95,-0.5,0.35,1000),
#'                              BestAvgResponse.lower = c(-1000,-0.4,-0.2,0.3), 
#'                              BestAvgResponse.upper = c(-0.4,-0.2,0.3,1000), 
#'                              Level = c( 'CR','PR', 'SD','PD'))
#' npdxe.criteria
#' 
#' data(TAN.volume.data)
#' TAN.drl <- DRLevel(data = TAN.volume.data, 
#'                    method = 'NPDXE.Response', 
#'                    criteria = npdxe.criteria, 
#'                    neg.control = 'Control')
#' plotDRResults(TAN.drl,by = 'Arms',measurement = 'DCR',min.number = 5)
#' 
#' @export
plotDRResults <-function(data, by, min.number, measurement = c('RR','DCR','both')){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  drl.res <- DRLevelAnalysis(data, by, measurement)
  drl.res <- drl.res[drl.res$Number >= min.number, ]

  if(measurement == 'both') {
    drl.res <- drl.res[order(drl.res$RR,decreasing = T),]

    rr.data <- drl.res[,colnames(drl.res) %in% setdiff(colnames(drl.res),'DCR')]
    rr.data$measurement <-'RR'
    colnames(rr.data)[which(colnames(rr.data) == 'RR')] <- 'Rate'
    dcr.data <- drl.res[,colnames(drl.res) %in% setdiff(colnames(drl.res),'RR')]
    dcr.data$measurement <-'DCR'
    colnames(dcr.data)[which(colnames(dcr.data) == 'DCR')] <- 'Rate'
    drl.res <- rbind(rr.data,dcr.data)
    drl.res$measurement <- factor(drl.res$measurement,levels = c('RR','DCR'))

    x.order <- factor(drl.res[,by],levels = rr.data[,by])
    Rate <- drl.res[ , 'Rate']

    p <- ggplot(drl.res,aes(x = x.order,y = Rate,fill = measurement))+
      geom_bar(position = "dodge",stat = "identity")

  }
  else {
    drl.res <- drl.res[order(drl.res[,measurement],decreasing = T),]

    colnames(drl.res)[which(colnames(drl.res) == measurement)] <-'Rate'
    
    x.order <- factor(drl.res[,by],levels = drl.res[,by])
    Rate <- drl.res[ , 'Rate']

    p <-ggplot(drl.res,aes(x = x.order,y = Rate))+
      geom_bar(position = "dodge",stat = "identity")
  }

  p <- p + xlab("") + ylab("Rate (%)")

  p <- p + theme_bw()+
    theme(
      panel.background = element_rect(fill = "transparent"),
      panel.grid.minor = element_blank(),
      panel.grid.major = element_blank(),
      plot.background  = element_rect(fill = "transparent")
    )   #backgroud

  p <- p + theme(
    plot.margin  = unit(c(1,1,1,3), 'lines'),
    axis.title.x = element_text(face = "bold",size = 12),
    axis.text.x  = element_text(hjust = 1,vjust = 1,size = 10,angle = 45),
    axis.title.y = element_text(face = "bold",size = 12),
    axis.text.y  = element_text(hjust = 0,size = 12),
  )

  plot(p)

}


####################################################################################################
#' @name plotWaterfall
#' @title  plot waterfall 
#' @description Visulize the drug response for each animal or each Tumor in the form of waterfall.
#'
#' @param data drug response level, the output of DRLevel.
#' @param max.threshold, the max threshold of volume change meseaure to visulize.
#' @param order.by the meseasure ordered by.
#' 
#' @import ggplot2
#' 
#' @examples
#' \donttest{
#' ## NPDXE.criteria
#' npdxe.criteria <- data.frame(BestResponse.lower = c(-1000,-0.95,-0.5,0.35),
#'                              BestResponse.upper = c(-0.95,-0.5,0.35,1000),
#'                              BestAvgResponse.lower = c(-1000,-0.4,-0.2,0.3), 
#'                              BestAvgResponse.upper = c(-0.4,-0.2,0.3,1000), 
#'                              Level = c( 'CR','PR', 'SD','PD'))
#' npdxe.criteria
#' 
#' data(TAone.volume.data)
#' TAone.drl <- DRLevel(data = TAone.volume.data, 
#'                      method = 'NPDXE.Response', 
#'                      criteria = npdxe.criteria, 
#'                      neg.control = 'Control')
#' 
#' # Choose the data of GC responsed to BYL719 + LJM716
#' bl.gc.drl <- TAone.drl[TAone.drl$Arms == "BYL719 + LJM716" & TAone.drl$Type == 'GC',]
#' plotWaterfall(data = bl.gc.drl, max.threshold = 100,order.by = 'Best.Response')
#' }
#' 
#' @export
plotWaterfall <- function(data, max.threshold = 100,
                          order.by = c('Best.Avg.Response','Best.Response','RC.Response')){
  
  if(class(data) != 'data.frame'){
    data<-as.data.frame(data)
  }
  
  colnames(data)[which(colnames(data)==order.by)] <- 'order.by'
  data <- data[order(data$order.by,decreasing = T),]
  data$order.by <- data$order.by*100

  aa <- max.threshold
  bb <- -max.threshold

  data$order.by[data$order.by > aa] <- aa
  data$order.by[data$order.by < bb] <- bb

  data$orders <- as.factor(1:nrow(data))

  #control the order of legend
  colnames(data)[grep('Level',colnames(data),fixed = T)] <- 'Level'

  levels <- unique(as.character(data$Level))
  
  orders <- data[ , 'orders']
  Level <- data[ , 'Level']

  p <- ggplot(data, aes(x = orders, y = order.by,fill = Level)) +
    scale_fill_discrete(name = "Response Level",limits = levels) +
    scale_color_discrete(guide = "none") +
    labs(list(x = NULL,y = "Tumor Volume Change (%)")) +
    theme_classic() %+replace%
    theme(axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.ticks.x = element_blank(),
          axis.text.y  = element_text(hjust = 1,vjust = 0.5,size = 12),
          axis.title.y = element_text(face = "bold",angle = 90,size = 12)) +
    coord_cartesian(ylim = c(min(data$order.by),max(data$order.by)))

  p <- p + geom_bar(stat = "identity", width = 0.8, position = position_dodge(width = 0.8))

  p <- p+geom_abline(intercept = 0, slope = 0,size = 0.5,colour = 'black')

  p <- p + theme(legend.position = c(0.95,0.95),legend.justification = c(0.95,0.95)) +
    theme(legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.title = element_text( size = 12, face = "bold"),
          legend.text = element_text( size = 12))
  plot(p)

}


####################################################################################################
#' @name oneAN.volume.data
#' @docType data
#' @title Example of volume data for oneAN pattern.
#' @description An example of volume data for oneAN pattern, which includes the required imformation for DRAP package.
#' @usage data(oneAN.volume.data)
#' 
#' @details Tumor volume data for oenAN pattern must at least include the columns "Arms" "ID" "Times" and "Volume".
#' 
#' @examples
#' data(oneAN.volume.data)
#' oneAN.volume.data[1:10,]
#' 
#' @keywords datasets
#' @export


####################################################################################################
#' @name oneAN.bw.data
#' @docType data
#' @title Example of body weight data for oneAN pattern.
#' @description An example of body weight data for oneAN pattern, which includes the required imformation for DRAP package.
#' @usage data(oneAN.bw.data)
#' 
#' @details Body weigth data for oenAN pattern must at least include the columns "Arms" "ID" "Times" and "BodyWeight".
#' 
#' @examples
#' data(oneAN.bw.data)
#' oneAN.bw.data[1:10,]
#' 
#' @keywords datasets
#' @export


####################################################################################################
#' @name TAone.volume.data
#' @docType data
#' @title Example of volume data for TAone pattern.
#' @description An example of volume data for TAone pattern, which includes the required imformation for DRAP package.
#' @usage data(TAone.volume.data)
#' 
#' @details The dataset derived from Novartis Institutes for BioMedical Research PDX encyclopedia (NIBR PDXE), which includes tumor volume data for 6 tumor types, 277 tumors, and total 4771 animals responded to 61 treatments. The dataset used in DRAP includes information "Tumor" "ID" "Type" "Arms" "Times" "Volume".
#' @references 
#' Gao, H., et al. High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response. Nat Med 2015;21(11):1318-1325.
#' 
#' @examples
#' data(TAone.volume.data)
#' TAone.volume.data[1:10,]
#' 
#' @keywords datasets
#' @export


####################################################################################################
#' @name TAone.bw.data
#' @docType data
#' @title Example of body weight data for TAone pattern.
#' @description An example of body weight data for TAone pattern, which includes the required imformation for DRAP package.
#' @usage  data(TAone.bw.data)
#' 
#' @details The dataset derived from Novartis Institutes for BioMedical Research PDX encyclopedia (NIBR PDXE), which includes body weight data for 6 tumor types, 277 tumors, and total 4771 animals responded to 61 treatments. The dataset used in DRAP includes information "Tumor" "ID" "Type" "Arms" "Times" "BodyWeight".
#' @references 
#' Gao, H., et al. High-throughput screening using patient-derived tumor xenografts to predict clinical trial drug response. Nat Med 2015;21(11):1318-1325.
#' 
#' @examples
#' data(TAone.bw.data)
#' TAone.bw.data[1:10,]
#' 
#' @keywords datasets
#' @export


####################################################################################################
#' @name TAN.volume.data
#' @docType data
#' @title Example of volume data for TAN pattern.
#' @description An example of volume data for TAN pattern, which includes the required imformation for DRAP package.
#' @usage  data(TAN.volume.data)
#' 
#' @details Tumor volume data for TAN pattern must at least include the columns "Arms" "Tumor" "ID" "Times" and "Volume".
#' 
#' @examples
#' data(TAN.volume.data)
#' TAN.volume.data[1:10,]
#' 
#' @keywords datasets
#' @export


####################################################################################################
#' @name TAN.bw.data
#' @docType data
#' @title Example of body weight data for TAN pattern.
#' @description An example of body weight data for TAN pattern, which includes the required imformation for DRAP package.
#' @usage  data(TAN.bw.data)
#' 
#' @details Body weigth data for TAN pattern must at least include the columns "Arms" "Tumor" "ID" "Times" and "BodyWeight".
#' 
#' @examples
#' data(TAN.bw.data)
#' TAN.bw.data[1:10,]
#' 
#' @keywords datasets
#' @export

