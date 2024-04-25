#install.packages("rEDM")
#library(rEDM)
data(sardine_anchovy_sst)

#anchovy vs sst
anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, lib_column = "anchovy",target_column = "np_sst", lib_sizes = seq(10, 60, by = 10), num_samples = 1000, random_libs = TRUE, replace = TRUE, silent = TRUE)
sst_xmap_anchovy <- ccm(sardine_anchovy_sst, E = 3, lib_column = "np_sst", target_column = "anchovy",lib_sizes = seq(10, 60, by = 10), num_samples = 1000, random_libs = TRUE,replace = TRUE, silent = TRUE)
a_xmap_t_means <- ccm_means(anchovy_xmap_sst)
t_xmap_a_means <- ccm_means(sst_xmap_anchovy)

plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red",xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.25))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topleft", legend = c("anchovy xmap SST", "SST xmap anchovy"), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)


#sardine vs sst
sardine_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, lib_column = "sardine",target_column = "sio_sst", lib_sizes = seq(10, 60, by = 10), num_samples = 1000, random_libs = TRUE, replace = TRUE, silent = TRUE)
sst_xmap_sardine <- ccm(sardine_anchovy_sst, E = 3, lib_column = "sio_sst", target_column = "sardine",lib_sizes = seq(10, 60, by = 10), num_samples = 1000, random_libs = TRUE,replace = TRUE, silent = TRUE)
s_xmap_t_means <- ccm_means(sardine_xmap_sst)
t_xmap_s_means <- ccm_means(sst_xmap_sardine)

plot(s_xmap_t_means$lib_size, pmax(0, s_xmap_t_means$rho), type = "l", col = "red",xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.4))
lines(t_xmap_s_means$lib_size, pmax(0, t_xmap_s_means$rho), col = "blue")
legend(x = "topleft", legend = c("sardine xmap SST", "SST xmap sardine"), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)


#sardine vs anchovy
sardine_xmap_anchovy <- ccm(sardine_anchovy_sst, E = 3, lib_column = "sardine",target_column = "anchovy", lib_sizes = seq(10, 60, by = 10), num_samples = 1000, random_libs = TRUE, replace = TRUE, silent = TRUE)
anchovy_xmap_sardine <- ccm(sardine_anchovy_sst, E = 3, lib_column = "anchovy", target_column = "sardine",lib_sizes = seq(10, 60, by = 10), num_samples = 1000, random_libs = TRUE,replace = TRUE, silent = TRUE)
s_xmap_a_means <- ccm_means(sardine_xmap_anchovy)
a_xmap_s_means <- ccm_means(anchovy_xmap_sardine)

plot(s_xmap_a_means$lib_size, pmax(0, s_xmap_a_means$rho), type = "l", col = "red",xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.4))
lines(a_xmap_s_means$lib_size, pmax(0, a_xmap_s_means$rho), col = "blue")
legend(x = "topleft", legend = c("sardine xmap anchovy", "anchovy xmap sardine"), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)


#sardine and  anchovy vs year
normalized_sardine=(sardine_anchovy_sst$sardine+abs(min(sardine_anchovy_sst$sardine)))/(max(sardine_anchovy_sst$sardine)-min(sardine_anchovy_sst$sardine))*4
normalized_anchovy=(sardine_anchovy_sst$anchovy+abs(min(sardine_anchovy_sst$anchovy)))/(max(sardine_anchovy_sst$anchovy)-min(sardine_anchovy_sst$anchovy))*4
plot(sardine_anchovy_sst$year, normalized_sardine, type = "l", col = "red",xlab = "Year", ylab = "fisheries landings")
lines(sardine_anchovy_sst$year, normalized_anchovy, col = "blue")
legend(x = "topleft", legend = c("sardine", "anchovy"), col = c("red","blue"), lwd = 1, bty = "n", inset = 0.02, cex = 0.8)


#attractor
x=sardine_anchovy_sst$anchovy
plot(head(x,-1),tail(x,-1),type="l")
