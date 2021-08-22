#' Runs PMMH procedure for the toy SV model
#'
#' @param data measurements or y
#' @param initVals starting values as a vector in the follwoing order
#'   phiX, sigmaX, and betaY
#' @param particles number of particles
#' @param resampleFreq frequency at which resampling is performed; if negative,
#'   then resampling is never performed; if between [0,1), then resampling is
#'   performed when  the ESS falls below that proportion of the number of
#'   particles and when it is greater than or equal to 1, resampling is carried
#'   out when the ESS falls below that value (note: if this parameter is larger
#'   than the total number of particles, then resampling will always be
#'   performed!)
#'
#' @return returns PMMH output
#' @export
svModelALTracking <- function(data,
                              initVals,
                              particles = 8,
                              resampleFreq = 0.5,
                              printALinfo = FALSE) {
  if (missing(data)) {
    msgErr <- paste0("'data' argument is missing: ",
                      "either real data or simulate e.g. from ",
                      "'SVmodelPMMH::generateDataSimulSv()'")
    stop(msgErr)
  }
  if (particles > 8 ) {
    warning("Number of particles exceeds number of  different
            colours used in the output plot; bad visibility of
            traced ancestral lines.")
  }
  
  res <- SVmodelRcppSMC::svModelALtrackingImp(as.matrix(data),
                                              initVals,
                                              particles,
                                              resampleFreq,
                                              printALinfo) 
  plotAL <- data.frame(iteration = rep(1:length(data), times = particles),
                       particle_no = paste0("particle ", rep(1:particles, each = length(data))),
                       particle_value = as.vector(t(res[["ancestorLines"]])))
  # plotAL <- t(res[["ancestorLines"]])
  # matplot(plotAL, type = c("b"),pch = 1, col = 1:particles) #plot
  # legend("topleft", legend = 1:particles, col=1:particles, pch=1)
  # title(main = "Particle ancestral lines (full genealogy)",
  #       xlab = "Iteration",
  #       ylab = "Particle values")
  cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73",
                  "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
  cbbPalette <- rep(cbbPalette, length.out = particles)

  plotOut <- ggplot2::ggplot(plotAL, ggplot2::aes(x = iteration, y = particle_value)) + 
    ggplot2::geom_line(ggplot2::aes(color = particle_no, linetype = particle_no)) + 
    ggplot2::scale_color_manual(values = cbbPalette) +
    ggplot2::labs(x = "Iteration",
                  y = "Particle values",
                  title = "Particle ancestral lines (full genealogy)") + 
    ggplot2::theme_bw() + ggplot2::theme(
        plot.title = ggplot2::element_text(face = "bold", size = 12),
        axis.ticks = ggplot2::element_line(colour = "grey70", size = 0.2),
        panel.grid.major = ggplot2::element_line(colour = "grey70", size = 0.2),
        panel.grid.minor = ggplot2::element_blank(),
        legend.background = ggplot2::element_rect(fill = "white", size = 4, colour = "white"))
  
  plot(plotOut)
  
  return(invisible(res))
}