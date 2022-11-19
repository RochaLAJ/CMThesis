library(shiny)
library(rvest)
library(dplyr)
library(ggplot2)
library(stringi)
library(car)
library(RColorBrewer)

LCR <- read.csv('./data/FullList_LCR_filtered.csv', header=T, row.names=1)
DM <- read.csv('./data/FullList_DM_filtered.csv', header=T, row.names=1)

biotype_meaningLCR <- read.csv('./data/biotypesLCR.csv', header=T, row.names=1)
biotype_meaningLCR$Biotypes <- factor(biotype_meaningLCR$Biotypes, 
						       levels = biotype_meaningLCR$Biotypes[order(biotype_meaningLCR$Freq, decreasing = T)]) 

biotype_meaningDM <- read.csv('./data/biotypesDM.csv', header=T, row.names=1)
biotype_meaningDM$Biotypes <- factor(biotype_meaningDM$Biotypes, 
							   levels = biotype_meaningDM$Biotypes[order(biotype_meaningDM$Freq, decreasing = T)]) 

top10LCR <- read.csv('./data/top10_table_LCR.csv', header=T, row.names=1)
top10DM <- read.csv('./data/top10_table_DM.csv', header=T, row.names=1)


server <- function(input, output, session) {
       output$tab1 = renderTable({
        			results_top10 <- switch(input$radio1,
		  			            	        A = top10LCR,
		    					            B = top10DM)
                    results_top10 %>%
                    	 filter(Biotypes==input$biotypeTable) %>%
                    	 arrange(desc(Average)) %>%
                    	 top_n(10)
		})

      	output$selected_rows = renderTable({
      	 	if (is.null(input$plot_click$y)) return()
        	 else {
       				 keeprows <- round(input$plot_click$y) == as.numeric(biotype_meaningLCR$Biotypes)
       				 results = biotype_meaningLCR[keeprows, ]
       				 results[4]
       			  }
    	})

    	output$x_value = renderText({
      		if (is.null(input$plot_click$y)) return()
      		 else {
      				lvls <- levels(biotype_meaningLCR$Biotypes)
      				name <- lvls[round(input$plot_click$y)]
      				HTML("You've selected <code>", name, "[Source: https://www.gencodegenes.org/pages/biotypes.html]")
      			  }
      	})
 
		output$freqplot = renderPlot({
		    	     biotype_meaning <- switch(input$radio1,
		    							       A = biotype_meaningLCR,
		    							       B = biotype_meaningDM)
	     			 ggplot(biotype_meaning,
	 	        	 aes(x =  Biotypes,Percentage, fill=brewer.pal(12,'Set3'), alpha=0.97)) + 
	 	 			 geom_bar(stat = "identity",colour='darkgreen') + 
	 	 			 theme(legend.position='none',panel.background = element_rect(fill='transparent'),
	 	 									 plot.background  = element_rect(fill='transparent')) +
	 	 			 coord_flip() +
	 	 			 ggtitle("Percentage Distribution By Biotype") +
	 	 			 xlab("Biotypes") + ylab("Percentage Value")
	 	  }, bg='transparent', res=100)
	 
     output$pubmed =  renderTable({ 
     		if (is.null(input$search)) return()
     	  	  else {
	  				 link = 'https://www.ncbi.nlm.nih.gov/gene/'
	  				 pub_search <- LCR %>% 
	  						  select(Symbol, EntrezID) %>% 
	  						  filter(Symbol==input$search)
	  				 stri_sub(link, nchar(link)+1) <- pub_search[2]
	  				 read 		<- read_html(link)
	  				 description <-  read %>%
	        		  		  html_nodes("dd") %>% 
	  		          		  html_text()
	  		    	 description[9:11]
	  		       }
	    })


	  output$box = renderPlot({
	  	     if (is.null(input$biotype)) return()
	  		   else {
	  		         data <- switch(input$radio1,
	  		   		  	     	        A = LCR,
	  		   		    	            B = DM)
	  	             box <- data %>% filter(Biotypes==input$biotype)  %>%
	  	     		 select(NED1, NED2, NED3, NED4, 
	  	     		        DS1, DS2, DS3, DS4, Symbol)
	  	             Boxplot(log2(box[1:8]), col=rep(c("red", "blue"), 
	  	        		times = c(4, 4)), pch=19, id=list(labels=box$Symbol,n=5),
	  	        		xlab='Samples', ylab='Counts Adjusted by Log2',
	  	        		main='Outliers')
	  	        	}
	  }, bg="transparent", res=80)
	}
