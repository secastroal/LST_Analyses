#63 mean age of 44.6 years and SD= 9.7; 79%
#female) and a control group (66 participants assigned to a waiting
#                             list with a mean age of 43.2 years and SD= 9.5; 73% female

s1 <- 63
s2 <- 66
m1 <- 44.6
m2 <- 43.2
sd1 <- 9.7
sd2 <- 9.5

gm <- (m1*s1 + m2*s2)/(s1 + s2)

var1 <- sd1 ^ 2
var2 <- sd2 ^ 2

ess1 <- var1 * (s1 - 1)
ess2 <- var2 * (s2 - 1)

ess <- ess1 + ess2

gss1 <- ((m1 -gm) ^ 2)* s1
gss2 <- ((m2 -gm) ^ 2)* s2

gss <- gss1 + gss2

gvar <- (ess + gss)/ (s1 + s2 -1)

gsd <- sqrt(gvar)

gm; gsd

s1*0.79
s2 * 0.73

(48+49) / 129

(0.79 * s1 + 0.73 * s2)/ 129
