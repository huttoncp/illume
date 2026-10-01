## POST HOC (chosen after the registered run, not a verdict): the floor
## variant. A candidate whose estimate is already more than e^6 below its line
## is held by value; the others are judged by the pushes, as ruled.
d <- utils::read.csv(commandArgs(TRUE)[1]); gs <- utils::read.csv(commandArgs(TRUE)[2])[, c("seed", "rel")]
d <- d[d$candidate, ]
## distance below the line, on the log scale the pushes use
d$below <- with(d, ifelse(study == "disp" | (study == "cases" & cell %in% 1:3),
                          (-value / 2) - log(1e-2),            # log(1/sqrt(disp)) - log(line)
                   ifelse(study == "re", NA, NA)))
g <- merge(d[d$study == "gauss", c("study", "cell", "rep", "seed", "path")], gs, by = "seed")
d$below[d$study == "gauss"] <- log(g$rel[match(paste(d$seed[d$study == "gauss"], d$path[d$study == "gauss"]),
                                               paste(g$seed, g$path))] / 0.2)
## the random effects: the SD on main's scale against 0.1 (of sd(y) for gaussian), from the value
d$below[d$study == "re"] <- NA      # none of the flips or unconverged fits is a random effect
far <- !is.na(d$below) & d$below < -6
d$floor <- ifelse(far, "held by value", d$new)
k <- function(x) paste(x$study, x$cell, x$rep)
p1 <- d[d$path == "P1", ]; p2 <- d[d$path == "P2", ]
m <- merge(p1[, c("study", "cell", "rep", "new", "floor")], p2[, c("study", "cell", "rep", "new", "floor")],
           by = c("study", "cell", "rep"), suffixes = c(".1", ".2"))
cat("flips across paths: new", sum(m$new.1 != m$new.2), " floor variant", sum(m$floor.1 != m$floor.2), "\n")
cat("unconverged: new", sum(d$new == "unconverged"), " floor variant", sum(d$floor == "unconverged"), "\n")
cat("held by value:", sum(d$floor == "held by value"), "by study/path:\n"); print(table(d$study[far], d$path[far]))
print(d[d$floor == "unconverged", c("path", "study", "cell", "rep", "value", "push_6", "conv")], row.names = FALSE, digits = 4)
