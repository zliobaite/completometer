# 2021 03 26 I.Zliobaite

Nnow <- 6363
Nfos <- 4501+742
d <- 1
t <- 23-5.333

C <- (Nfos/Nnow)*(d/(d+t))

dd <- seq(0,4,0.01)

CC <- (Nfos/Nnow)*(dd/(dd+t))

pdf('fig_sensitivity1.pdf',height = 3.5,width=4)
plot(dd,CC*100,type='l',lwd=3,xlab = 'Species duration, Myr',ylab = 'Completeness, %',ylim = c(0,16))
dev.off()

alf <- 2
bet <- ((alf + (1-alf)*(23-5.333)/23) + alf)/2

C2 <- Nfos/(Nnow*alf + Nnow*bet*t/d)

aa <- seq(0.5,2,0.01)
bb <-  ((aa + (1-aa)*(23-5.333)/23) + aa)/2

CC2 <- Nfos/(Nnow*aa + Nnow*bb*t/d)

pdf('fig_sensitivity2.pdf',height = 3.5,width=4)
plot(aa,CC2*100,type='l',lwd=3,xlab = '(increasing) -- Diversity -- (decreasing)',ylab = 'Completeness, %',ylim = c(0,16),log="x")
dev.off()
