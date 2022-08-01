TM.growth(1,1,0:220,elas=1,sens=2)

TM.survival(-2,0:220,elas=1,sens=2)

plot(0:220,
TM.flower(9.6,0:220,elas=1,sens=2)
)

plot(0:220,
     TM.seeds(14,0:220,elas=1,sens=2)
)

plot(0:220,
TM.recruitment(d=0:220,elas=1,sens=2)
)


hold.lambda<-c()
for(i in 1:10){
  # "06_BootRes"
  # Run resampling subroutine for wind speeds, terminal velocities, and demography
  set.seed(seeds[i,2]);print(seeds[i,2])
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/06_BootRes.R")
  print(mean(LATR_full$max.ht_t,na.rm=T))
  # "05_CDataAnalysis_NS.R"
  # Create demography models for use in SIPM, using best models from initial run of this file
  # This avoids model uncertainty; only model coefficients change, not the overall structure
  source("https://raw.githubusercontent.com/TrevorHD/LTEncroachment/master/05_CDataAnalysis.R")
  hold.lambda[i] <- lambda(TransMatrix(dens = 0, mat.size = 100)$IPMmat)
}
