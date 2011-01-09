png('Rice2008-modmlversion-isosarcometricfca.png', width = 480, height = 480, bg='white')
d <- read.csv('Isosarcometric-2.3.csv')
plot(d$Amount.of.Calcium..2...concentration.in.Cardiac.Muscle.Site * 1000, d$Active.Force, type='l', log='x', main='Isosarcometric F-Ca graph for Rice2008 model\n(from ModML-Solver)', ylab='Normalised active force (unitless)', xlab='[Ca](microM)')
for (i in 14:23/10) { d <- read.csv(sprintf('Isosarcometric-%0.1f.csv', i)); lines(d$Amount.of.Calcium..2...concentration.in.Cardiac.Muscle.Site * 1000, d$Active.Force, type='l') }
dev.off()

png('Rice2008-modmlversion-isometricfca.png', width = 480, height = 480, bg='white')
d <- read.csv('Isometric-50.csv')
plot(d$Amount.of.Calcium..2...concentration.in.Cardiac.Muscle.Site * 1000, d$Active.Force, type='l', log='x', main='Isometric F-Ca graph for Rice2008 model\n(from ModML-Solver)', ylab='Normalised active force (unitless)', xlab='[Ca](microM)')
for (kse in c("1.0", "1.4", "2.0", "3.0", "5.0", "10", "50")) { df <- read.csv(sprintf('Isometric-%s.csv', kse)); lines(df$Amount.of.Calcium..2...concentration.in.Cardiac.Muscle.Site * 1000, df$Active.Force); }
dev.off()
