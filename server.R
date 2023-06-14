library(shiny)

function(input, output, session) {

################################################################################################################################
### plot of RE versus number of time periods
################################################################################################################################  
  output$plot1 <- renderPlotly({
    
    ICC=input$ICC1                    # intraclass correlation coefficient
    CAC=input$CAC1                    # cluster autocorrelation

    nr.subj=input$nr.subj1             # number of subjects per cluster-period
    min.nr.per=input$nr.per1[1]       # minimum number of periods to be evaluated
    max.nr.per=input$nr.per1[2]       # maximum number of periods to be evaluated

    costs.c=input$costs.c1            # costs per cluster
    costs.s=input$costs.s1            # costs per subject
    costs.x=input$costs.x1            # costs per crossover
    budget=input$budget1              # budget
    
    validate(
    need(input$costs.c1>0, 'Input error: costs per cluster should be larger than zero'),
    need(input$costs.s1>0, 'Input error: costs per subject should be larger than zero'),
    need(input$costs.x1>0, 'Input error: costs per treatment switch should be larger than zero'),
    need(input$budget1>0,  'Input error: budget should be larger than zero'),
    )
        
    results=matrix(NA,12,6)                                                                                   # matrix to store results
    
    for(ii in min.nr.per:max.nr.per)
    {
      nr.per=ii                                                                                               # number of periods
      
      ### calculation of covariance matrix
      exponent <- abs(matrix(1:nr.per - 1, nrow = nr.per, ncol = nr.per, byrow = TRUE) - (1:nr.per - 1))
      VV <- ICC*CAC^exponent + diag(1-ICC,nrow=nr.per)/nr.subj                                                # correlation matrix
      VVinv=solve(VV)                                                                                         # inverse correlation matrix
      
      ### calculation of design matrix
      sequence.A=rep(0:1,times=6)
      sequence.B=1-sequence.A
      XX.A=cbind(rep(1,nr.per),diag(1,nrow=nr.per)[,2:nr.per],sequence.A[1:nr.per])
      XX.B=cbind(rep(1,nr.per),diag(1,nrow=nr.per)[,2:nr.per],sequence.B[1:nr.per])
      var.gamma1=solve((t(XX.A)%*%VVinv%*%XX.A+t(XX.B)%*%VVinv%*%XX.B))                                       # covariance matrix
      costs=costs.c+costs.s*nr.subj*nr.per+costs.x*(nr.per-1)                                                 # costs in each cluster for current sequence    
      N=floor(0.5*budget/costs)
      var.gamma2=solve((N*t(XX.A)%*%VVinv%*%XX.A+N*t(XX.B)%*%VVinv%*%XX.B))                                   # covariance matrix
      results[ii,]=(c(nr.per,var.gamma1[nr.per+1,nr.per+1],costs*2*N,var.gamma2[nr.per+1,nr.per+1],NA,NA))
    }
    
    results[min.nr.per:max.nr.per,5]=round(min(results[min.nr.per:max.nr.per,2])/results[min.nr.per:max.nr.per,2],3)   # RE given fixed number of clusters
    results[min.nr.per:max.nr.per,6]=round(min(results[min.nr.per:max.nr.per,4])/results[min.nr.per:max.nr.per,4],3)   # RE given budgetary constraint
    colnames(results)=c("nr.per","variance 1","costs","variance 2","RE 1","RE 2")                             # column names in results matrix 

    results=as.data.frame(results)
    key <- row.names(results)
    plot=plot_ly(results, x = ~results[,1], y = ~results[,5], mode='lines+markers',type="scatter",marker=list(size=10,color='red'),line=list(width=0.1,color='black'),name="fixed number of clusters") 
      plot=plot %>% add_lines(y=~results[,6],key = ~key,mode='lines',type="scatter",marker=list(size=10,color='blue'),line=list(color='black'),name="budgetary constraint")
      plot=plot %>% layout(dragmode = "select",
             xaxis = list(title = "Number of time periods", range = c(0, 12.5), showgrid = F,zeroline=TRUE,showline = TRUE),
             yaxis = list(title = "Relative Efficiency", range=c(0,1.05),showgrid = T,zeroline=TRUE,showline = TRUE))
  })
 
################################################################################################################################
### table of RE versus number of time periods
################################################################################################################################  
  output$ResultsTable <- DT::renderDataTable({
    
    ICC=input$ICC1                    # intraclass correlation coefficient
    CAC=input$CAC1                    # cluster autocorrelation
    
    nr.subj=input$nr.subj1             # number of subjects per cluster-period
    min.nr.per=input$nr.per1[1]       # minimum number of periods to be evaluated
    max.nr.per=input$nr.per1[2]       # maximum number of periods to be evaluated
    
    costs.c=input$costs.c1            # costs per cluster
    costs.s=input$costs.s1            # costs per subject
    costs.x=input$costs.x1            # costs per crossover
    budget=input$budget1              # budget
    
    validate(
      need(input$costs.c1>0, 'Input error: costs per cluster should be larger than zero'),
      need(input$costs.s1>0, 'Input error: costs per subject should be larger than zero'),
      need(input$costs.x1>0, 'Input error: costs per treatment switch should be larger than zero'),
      need(input$budget1>0,  'Input error: budget should be larger than zero'),
    )
    
    results=matrix(NA,12,6)                                                                                   # matrix to store results
    
    for(ii in min.nr.per:max.nr.per)
    {
      nr.per=ii                                                                                               # number of periods
      
      ### calculation of covariance matrix
      exponent <- abs(matrix(1:nr.per - 1, nrow = nr.per, ncol = nr.per, byrow = TRUE) - (1:nr.per - 1))
      VV <- ICC*CAC^exponent + diag(1-ICC,nrow=nr.per)/nr.subj                                                # correlation matrix
      VVinv=solve(VV)                                                                                         # inverse correlation matrix
      
      ### calculation of design matrix
      sequence.A=rep(0:1,times=6)
      sequence.B=1-sequence.A
      XX.A=cbind(rep(1,nr.per),diag(1,nrow=nr.per)[,2:nr.per],sequence.A[1:nr.per])
      XX.B=cbind(rep(1,nr.per),diag(1,nrow=nr.per)[,2:nr.per],sequence.B[1:nr.per])
      var.gamma1=solve((t(XX.A)%*%VVinv%*%XX.A+t(XX.B)%*%VVinv%*%XX.B))                                       # covariance matrix
      costs=costs.c+costs.s*nr.subj*nr.per+costs.x*(nr.per-1)                                                 # costs in each cluster for current sequence    
      N=floor(0.5*budget/costs)
      var.gamma2=solve((N*t(XX.A)%*%VVinv%*%XX.A+N*t(XX.B)%*%VVinv%*%XX.B))                                   # covariance matrix
      results[ii,]=(c(nr.per,var.gamma1[nr.per+1,nr.per+1],costs*2*N,var.gamma2[nr.per+1,nr.per+1],NA,NA))
    }
    
    results[min.nr.per:max.nr.per,5]=min(results[min.nr.per:max.nr.per,2])/results[min.nr.per:max.nr.per,2]   # RE given fixed number of clusters
    results[min.nr.per:max.nr.per,6]=min(results[min.nr.per:max.nr.per,4])/results[min.nr.per:max.nr.per,4]   # RE given budgetary constraint
    results=as.data.frame(results)
    results=round(results[min.nr.per:max.nr.per,c(1,3,5,6)],3)
    colnames(results)=c("number of periods","total costs","RE fixed number of clusters","RE budgetary constraint")    # column names in results matrix 
    datatable(results,rownames=FALSE, options = list(pageLength = 25,dom = 'tl'))
  })
  
################################################################################################################################
### table of RE for each treatment switch
################################################################################################################################  
  output$ResultsTable2 <- DT::renderDataTable({

    ICC=input$ICC2                    # intraclass correlation coefficient
    CAC=input$CAC2                    # cluster autocorrelation
    
    nr.subj=input$nr.subj2            # number of subjects per cluster-period
    nr.per=input$nr.per2              # minimum number of periods to be evaluated
    
    costs.c=input$costs.c2            # costs per cluster
    costs.s=input$costs.s2            # costs per subject
    costs.x=input$costs.x2            # costs per crossover
    budget=input$budget2              # budget
    
    validate(
      need(input$costs.c1>0, 'Input error: costs per cluster should be larger than zero'),
      need(input$costs.s1>0, 'Input error: costs per subject should be larger than zero'),
      need(input$costs.x1>0, 'Input error: costs per treatment switch should be larger than zero'),
      need(input$budget1>0,  'Input error: budget should be larger than zero'),
    )

    nr.seq=2^nr.per                                                                                                     # number of sequences for chosen number of periods
    
    ### calculation of covariance matrix
    exponent <- abs(matrix(1:nr.per - 1, nrow = nr.per, ncol = nr.per, byrow = TRUE) - (1:nr.per - 1))
    VV <- ICC*CAC^exponent + diag(1-ICC,nrow=nr.per)/nr.subj                                                            # correlation matrix
    VVinv=solve(VV)                                                                                                     # inverse correlation matrix
    
    ### calculation of design matrix
    XX=cbind(rep(1,nr.per),diag(1,nrow=nr.per)[,2:nr.per])
    list_treat=list(0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1, 0:1)
    sequences=expand.grid(list_treat[1:nr.per])                                                                         # all possible treatment sequences
    
    results=matrix(NA,nr.seq/2,7)                                                                                       # matrix to store results
    for(ii in 1:(nr.seq/2))
    {
      nr.switches= sum(abs(diff(as.numeric(sequences[ii,]))))                                                           # number of switches in current sequence
      XX.A=as.matrix(cbind(XX,t(sequences[ii,])))                                                                       # design matrix clusters in group A
      XX.B=as.matrix(cbind(XX,t(sequences[nr.seq+1-ii,])))                                                              # design matrix clusters in group B
      var.gamma1=solve((t(XX.A)%*%VVinv%*%XX.A+t(XX.B)%*%VVinv%*%XX.B))                                                 # covariance matrix
      
      costs=costs.c+costs.s*nr.subj*nr.per+costs.x*nr.switches                                                          # costs in each cluster of current sequence    
      N=floor(0.5*budget/costs)
      var.gamma2=solve(N*(t(XX.A)%*%VVinv%*%XX.A+N*t(XX.B)%*%VVinv%*%XX.B))                                             # covariance matrix
      
      results[ii,]=(c(ii,nr.switches,var.gamma1[nr.per+1,nr.per+1],costs*2*N,var.gamma2[nr.per+1,nr.per+1],NA,NA))
    }
    
    results[,6]=round(min(results[,3])/results[,3],3)                                                                   # RE given fixed number of clusters
    results[,7]=round(min(results[,5])/results[,5],3)                                                                   # RE given budgetary constraint
    
    sequences[sequences==0]="A"
    sequences[sequences==1]="B"
    order=do.call(paste,sequences)
    results=data.frame(order[1:(length(order)/2)],results[,c(2,4,6,7)])
    colnames(results)=c("order of treatments","nr.switches","costs","RE fixed number of clusters","RE budgetary constraint")  # column names in results matrix 
    datatable(results,rownames=FALSE, options = list(pageLength = 2^nr.per,dom = 'tl'))
  })
  
}

