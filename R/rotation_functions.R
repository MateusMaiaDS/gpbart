# Rotation function
A <- function(theta) {
        S <- sin(theta)
        C <- cos(theta)
        matrix(c(
                C, -S,
                S, C),
        ncol = 2, nrow = 2, byrow = TRUE)
}
