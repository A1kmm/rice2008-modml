png('Rice2008-modmlversion-isosarcometricfca.png', width = 178, height = 200, res=150, units='mm', antialias='subpixel', bg='white')
d <- read.csv('Isosarcometric-2.3.csv')
plot(d$Amount.of.Calcium..2...concentration.in.Cardiac.Muscle.Site * 1000, d$Active.Force, type='l', log='x', main='Isosarcometric F-Ca graph for Rice2008 model\n(from ModML-Solver)', ylab='Normalised active force (unitless)', xlab='[Ca2+](microM)', ylim=c(0.0, 1.0), xlim=c(0.2, 30.0))
text(22, 1, 'SL(micron)')
for (i in 14:23/10) {
  d <- read.csv(sprintf('Isosarcometric-%0.1f.csv', i)); lines(d$Amount.of.Calcium..2...concentration.in.Cardiac.Muscle.Site * 1000, d$Active.Force, type='l')
  text(22, d$Active.Force[length(d$Active.Force)], i)
}
dev.off()

png('Rice2008-modmlversion-isometricfca.png', width = 178, height = 200, res=150, units='mm', antialias='subpixel', bg='white')
d <- read.csv('Isometric-50.csv')
plot(d$Amount.of.Calcium..2...concentration.in.Cardiac.Muscle.Site * 1000, d$Active.Force, type='l', log='x', main='Isometric F-Ca graph for Rice2008 model\n(from ModML-Solver)', ylab='Normalised active force (unitless)', xlab='[Ca2+](microM)', ylim=c(0.0, 0.85), xlim=c(0.2, 30.0))
text(22, 0.85, 'KSE')

for (kse in c("1.0", "1.4", "2.0", "3.0", "5.0", "10", "50")) {
  df <- read.csv(sprintf('Isometric-%s.csv', kse)); lines(df$Amount.of.Calcium..2...concentration.in.Cardiac.Muscle.Site * 1000, df$Active.Force);
  text(22, df$Active.Force[length(d$Active.Force)], kse)
}
dev.off()
