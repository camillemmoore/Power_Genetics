
power_envir.calc(P_e = 0.2, MAF = 0.1, N = 200, Case.Rate = 0.5, Alpha = 0.05, OR_G = 1.5, OR_E = 2, OR_GE = 1.8, Test.Model = "Dominant",True.Model="Dominant")
# N=200;Case.Rate = 0.5;MAF = 0.1;P_e = 0.2;OR_G = 1.5;OR_E = 2;OR_GE = 1.8;Alpha = 0.05;True.Model = "All";Test.Model = "All"
# this one hangs
power_envir.calc(N=234, Case.Rate = 0.4, MAF = 0.2, OR_G = c(1.2,4), OR_E = c(1.5,2,2.6), P_e = 0.2, OR_GE = c(1.4,2,4), Alpha = 0.05, True.Model = "All", Test.Model = "All")
# fixed!
#testing
# N=234;Case.Rate = 0.4;MAF = 0.2;P_e = 0.2;OR_G = c(1.2,4);OR_E = c(1.5,2,2.6);OR_GE = c(1.4,2,4);Alpha = 0.05;True.Model = "All";Test.Model = "All"

# error:
power_envir.calc(N=200, Case.Rate = 0.5, MAF = 0.2, OR_G = c(1.2,4), OR_E = c(1.5,2,2.6), P_e = 0.2, OR_GE = c(1.4,2,4), Alpha = 0.05, True.Model = "All", Test.Model = c("Dominant", "Recessive", "Additive"))
# N=200;Case.Rate = 0.5;MAF = 0.2;P_e = 0.2;OR_G = c(1.2,4);OR_E = c(1.5,2,2.6);OR_GE = c(1.4,2,4);Alpha = 0.05;True.Model = "All";Test.Model = "All"
power_envir.calc(N = c(100, 200, 1000), Case.Rate= c(0.1, 0.2, 0.5), MAF = c(0.05, 0.1, 0.2, 0.3), OR_G = c(1.2, 2), OR_E = c(1.5, 2), OR_GE = c(1.5, 2), P_e = c(0.1, 0.2, 0.3), Alpha = 0.05, Test.Model="All", True.Model="All")

power_envir.calc(N = c(100, 200, 1000), Case.Rate= c(0.1, 0.2, 0.5), MAF = c(0.05, 0.1, 0.2, 0.3), OR_G = c(1.2, 2), OR_E = c(1.5, 2), OR_GE = c(1.5, 2), P_e = c(0.1, 0.2, 0.3), Alpha = 0.05, Test.Model="All", True.Model="All")

foo <- power_envir.calc(N = 100, Case.Rate= c(0.05), MAF = c(0.05, 0.1, 0.2), OR_G = c(0.1, 0.4), OR_E = c(0.05), OR_GE = c(0.1,0.8), P_e = c(0.1, 0.2, 0.3), Alpha = 0.05, Test.Model="All", True.Model="All");
#N=100;Case.Rate= 0.2;MAF = 0.05;OR_G = 0.2;OR_E = c(0.5,0.2);OR_GE = c(0.5,0.8);P_e = 0.2;Alpha = 0.05;Test.Model="All";True.Model="All"


for(ii in 1:3){
	# print(paste0("ii:", ii))
	for(jj in 1:4){
		# print(paste0("jj:", jj))
		for(kk in 1:2){
			# print(paste0("kk:", kk))
			for(ll in 1:2){
				# print(paste0("ll:", ll))
				for(nn in 1:2){
					# print(paste0("nn:", nn))
					for(mm in 1:3){
						print(paste0("ii:", ii, "; jj:", jj, "; kk:", kk, "; ll:", ll, "; nn:", nn, "; mm:", mm))
						foo <- power_envir.calc(N = c(100, 200, 1000)[ii], Case.Rate= c(0.5), MAF = c(0.05, 0.1, 0.2, 0.3)[jj], 
							OR_G = c(1.2, 2)[kk], OR_E = c(1.5, 2)[ll], OR_GE = c(1.5, 2)[nn], P_e = c(0.1, 0.2, 0.3)[mm], 
							Alpha = 0.05, Test.Model="All", True.Model="All", compareQuanto = 1)
					}
				}
			}
		}
	}
}
foo <- power_envir.calc(N = c(100, 200, 1000)[ii], Case.Rate= c(0.5), MAF = c(0.05, 0.1, 0.2, 0.3)[jj], 
							OR_G = c(1.2, 2)[kk], OR_E = 2, OR_GE = c(1.5, 2)[nn], P_e = c(0.1, 0.2, 0.3)[mm], 
							Alpha = 0.05, Test.Model="All", True.Model="All")

foo <- power_envir.calc(N = 100, Case.Rate= 0.2, MAF = 0.1, OR_G = 0.3, OR_E = c(0.05), OR_GE = c(0.1), P_e = c(0.2), Alpha = 0.05, Test.Model="All", True.Model="All");

# not working
ss_envir.calc(power= 0.08185038, Case.Rate = 0.3, P_e = 0.1, OR_G = 1.4, OR_E = 0.3, OR_GE = 2, MAF = 0.1, Alpha = 0.05, True.Model = "All", Test.Model = "All")
power= 0.08185038;Case.Rate = 0.3;P_e = 0.1;OR_G = 1.4;OR_E = 0.3;OR_GE = 2;MAF = 0.1;Alpha = 0.05;True.Model = "All";Test.Model = "All"
