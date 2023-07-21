library(shiny)
source('functions.R')
function(input, output, session) {
    
  
    dat.out <- reactive({
      
       w <- input$I.e + input$I.l
       
       la <- c(input$I.e/input$P.E, input$I.l/input$P.E, input$I.l/input$P.L)
       
       fu <- input$fu
       
       
       se.e.exp.e <- input$se.e.exp * (1-input$exp.ctam) + input$se.e.ctr * input$exp.ctam
       se.e.ctr.e <- input$se.e.exp * input$ct.ctam + input$se.e.ctr * (1 - input$ct.ctam)
       
       se.l.exp.e <- input$se.l.exp * (1-input$exp.ctam) + input$se.l.ctr * input$exp.ctam
       se.l.ctr.e <- input$se.l.exp * input$ct.ctam + input$se.l.ctr * (1 - input$ct.ctam)
       
       
       test.fu.tt = c(seq(from=0,by=input$t.interval,length.out=input$n.test),fu+2)
       
       prob.exp <- get.prob.s(se.e = se.e.exp.e, se.l = se.l.exp.e, test.fu.tt = test.fu.tt,
                              w = w,la = la)
       
       prob.ctr <- get.prob.s(se.e = se.e.ctr.e, se.l = se.l.ctr.e, test.fu.tt = test.fu.tt,
                              w = w,la = la)
      
       prob.exp.l <- prob.exp$all.late * (1 - input$non.scr) + 
         input$non.scr * input$I.l * prob.exp$time
       
       prob.ctr.l <- prob.ctr$all.late * (1 - input$non.scr) + 
         input$non.scr * input$I.l * prob.exp$time
       
       if(input$a.yr>0){
         prob.exp.l <- acurral.ave(prob=prob.exp.l, a.yr=input$a.yr)
         prob.ctr.l <- acurral.ave(prob=prob.ctr.l, a.yr=input$a.yr)
       }
       
       
       rel.red <- (1 - prob.exp.l / prob.ctr.l) * 100
       power <- power.appr(p1=prob.exp.l,p0=prob.ctr.l,n1=input$n.exp,n0=input$n.ctr,
                           alpha=input$alpha)
       
       dat <- data.frame(time=prob.exp$time,prob.exp.l=prob.exp.l,
                             prob.ctr.l=prob.ctr.l,rel.red=rel.red,
                             power=power)
       
      dat
    })

    
    output$cum.prob.plot <- renderPlot({
      
      dat <- dat.out()
      
      plot(dat$time, dat$prob.exp.l,type='l',col='red',
           ylab='Cumulative late-stage cancer probability',
           xlab='Years since 1st enrollment',lwd=3,lty=2,
           ylim=c(0,max(dat$prob.exp.l,dat$prob.ctr.l)),
           xlim=c(input$a.yr,input$fu+2))
      lines(dat$time, dat$prob.ctr.l,col='blue',lwd=3)
      legend('topleft',c('experimental arm', 'control arm'),lwd=3, lty=c(2,1),
             col = c('red','blue'))
      abline(v=input$fu,lty=3)
      title('absolute effect',adj=0)
    })
    
    output$rel.red.plot <- renderPlot({
      
      dat <- dat.out()
      
      plot(dat$time, dat$rel.red,type='l',col='black',
           ylab='relative reduction (%)',
           xlab='Years since 1st enrollment',lwd=3,lty=1,
           ylim=c(0,max(dat$rel.red)),
           xlim=c(input$a.yr,input$fu+2))
      abline(v=input$fu,lty=3)
      title('relative effect',adj=0)
    })
    
    
    output$power.plot <- renderPlot({
      
      dat <- dat.out()
      
      plot(dat$time, dat$power,type='l',col='black',
           ylab='power',
           xlab='Years since 1st enrollment',lwd=3,lty=1,
           ylim=c(0,1),
           xlim=c(input$a.yr,input$fu+2))
      abline(v=input$fu,lty=3)
    })
    
    
    output$power.rst<- renderUI({
      dat <- dat.out()
      dat.t <- dat[round(dat$time,2)==round(input$fu,2),]
      power <- round(dat.t['power'] * 100)
      rel.red <- round(dat.t['rel.red'])
      exp.l <- dat.t['prob.exp.l']
      ctr.l <- dat.t['prob.ctr.l']
      out <- paste0('Given the inputs, the study is expected to have a power of ',
                    power,'% at Year ', input$fu, ' to detect a relative reduction of ', rel.red, 
                    '% (=1-',signif(exp.l,2),'/',signif(ctr.l,2), ').')
      out
    })
    
    output$natural.history <- renderUI({
      w <- signif(input$I.e + input$I.l,2)
      
      la <- signif(c(input$I.e/input$P.E, input$I.l/input$P.E, input$I.l/input$P.L),2)
      withMathJax(paste0("$$w =", w, ", \\lambda_1 =", la[1],", \\lambda_2 =", la[2],
                         ", \\lambda_3 =", la[3],"$$"))
      
    })
    
}
