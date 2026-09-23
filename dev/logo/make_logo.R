## man/figures/logo.png from the sticker artwork in this folder.
##
## The artwork is a hex on a black square. Everything outside the hex's
## light border ring is made transparent, so the sticker sits cleanly on
## GitHub's light and dark themes and on the pkgdown site, then it is
## cropped to the hex and scaled for the README and the pkgdown home page.
##
## The ring's outer edge is anti-aliased against the black, so its outermost
## pixels are part ring, part black. Kept opaque they leave a dark fringe on
## a white page; dropped, the edge is jagged. They are "un-matted" instead:
## a pixel that is fraction a of the ring colour and 1 - a of black gets
## alpha a and the ring colour.
##
## Run from the package root:  Rscript dev/logo/make_logo.R
## Needs jpeg, png and magick; none of them is a dependency of the package.
##
## pkgdown/favicon/ was made from the result by pkgdown::build_favicons(),
## which sends the logo to realfavicongenerator.net -- and which
## pkgdown::init_site() runs by itself whenever there is a logo and no
## favicons, so keep them committed or every CI build calls the service. Two
## edits followed: favicon.svg embedded the whole 480 px logo (573 KB) for an
## icon drawn at 16-32 px, so it now embeds a 96 px copy (34 KB); and
## site.webmanifest's icon paths began with "/", which on a project site
## served under /illume/ points outside it, so they are relative. If the logo
## changes, delete pkgdown/favicon/, rebuild them, and repeat both.

src   <- "dev/logo/illume-hex-source.jpg"
out   <- "man/figures/logo.png"
width <- 480L       # px; displayed at 139 px tall in the README, so ~4x for sharp screens
thr   <- 50         # 0-255: brighter than this, a pixel belongs to the sticker

A <- jpeg::readJPEG(src)                    # rows x cols x 3, in [0, 1]
nr <- dim(A)[1]; nc <- dim(A)[2]
mx <- 255 * pmax(A[, , 1], A[, , 2], A[, , 3])
lum <- 255 * (0.2126 * A[, , 1] + 0.7152 * A[, , 2] + 0.0722 * A[, , 3])

## OUTSIDE. The hex is convex, so along any row the sticker is one run of
## pixels: what lies before its first bright pixel or after its last is
## outside. (The black inside the ring never counts: it lies within the run.)
bright <- mx > thr
outside <- matrix(TRUE, nr, nc)
for (r in seq_len(nr)) {
  on <- which(bright[r, ])
  if (length(on)) outside[r, min(on):max(on)] <- FALSE
}

## chessboard distance from the edge, both ways, out to 3 px
shift <- function(M, dr, dc, fill) {
  out <- matrix(fill, nr, nc)
  rs <- max(1, 1 + dr):min(nr, nr + dr); cs <- max(1, 1 + dc):min(nc, nc + dc)
  out[rs, cs] <- M[rs - dr, cs - dc]
  out
}
grow <- function(M) {
  G <- M
  for (dr in -1:1) for (dc in -1:1) G <- G | shift(M, dr, dc, FALSE)
  G
}
d_in <- matrix(Inf, nr, nc)    # inside pixels: how far from the outside
d_out <- matrix(Inf, nr, nc)   # outside pixels: how far from the inside
Go <- outside; Gi <- !outside
for (k in 1:3) {
  Go <- grow(Go); d_in[Go & !outside & !is.finite(d_in)] <- k
  Gi <- grow(Gi); d_out[Gi & outside & !is.finite(d_out)] <- k
}

## the ring's own colour near each edge pixel: the mean of ring pixels 3 px
## in, within a small window. The ring is shaded -- about 210 at the top and
## sides, 135 along the bottom -- so a single colour for all of it would not do.
ring <- which(d_in == 3 & mx > 100)
band <- which((d_in <= 2) | (d_out <= 2))
rr <- ((band - 1) %% nr) + 1; cc <- ((band - 1) %/% nr) + 1
acc <- matrix(0, length(band), 3); cnt <- numeric(length(band))
is_ring <- logical(nr * nc); is_ring[ring] <- TRUE
for (dr in -4:4) for (dc in -4:4) {
  r2 <- rr + dr; c2 <- cc + dc
  ok <- r2 >= 1 & r2 <= nr & c2 >= 1 & c2 <= nc
  idx <- (c2[ok] - 1) * nr + r2[ok]
  hit <- is_ring[idx]
  w <- which(ok)[hit]
  for (ch in 1:3) acc[w, ch] <- acc[w, ch] + A[, , ch][idx[hit]]
  cnt[w] <- cnt[w] + 1
}
has <- cnt > 0
ringcol <- acc[has, , drop = FALSE] / cnt[has]
ringlum <- 255 * (0.2126 * ringcol[, 1] + 0.7152 * ringcol[, 2] + 0.0722 * ringcol[, 3])

## assemble RGBA: opaque inside, clear outside, un-matted along the edge
alpha <- matrix(0, nr, nc); alpha[!outside] <- 1
R <- A[, , 1]; G <- A[, , 2]; B <- A[, , 3]
bi <- band[has]
a <- pmin(pmax(lum[bi] / ringlum, 0), 1)
alpha[bi] <- a
R[bi] <- ringcol[, 1]; G[bi] <- ringcol[, 2]; B[bi] <- ringcol[, 3]
cat(sprintf("edge pixels un-matted: %d (%d without a ring colour nearby, left as they were)\n",
            length(bi), sum(!has)))

## crop to the sticker, with a 2 px margin
keep <- which(alpha > 0.01, arr.ind = TRUE)
r0 <- max(1, min(keep[, 1]) - 2); r1 <- min(nr, max(keep[, 1]) + 2)
c0 <- max(1, min(keep[, 2]) - 2); c1 <- min(nc, max(keep[, 2]) + 2)
rgba <- array(0, c(r1 - r0 + 1, c1 - c0 + 1, 4))
rgba[, , 1] <- R[r0:r1, c0:c1]; rgba[, , 2] <- G[r0:r1, c0:c1]
rgba[, , 3] <- B[r0:r1, c0:c1]; rgba[, , 4] <- alpha[r0:r1, c0:c1]
cat(sprintf("cropped to %d x %d px (rows %d-%d, cols %d-%d)\n",
            c1 - c0 + 1, r1 - r0 + 1, r0, r1, c0, c1))

full <- tempfile(fileext = ".png")
png::writePNG(rgba, full)
img <- magick::image_read(full)
img <- magick::image_resize(img, sprintf("%dx", width), filter = "Lanczos")
dir.create(dirname(out), showWarnings = FALSE, recursive = TRUE)
magick::image_write(img, out, format = "png", quality = 95)
info <- magick::image_info(magick::image_read(out))
cat(sprintf("wrote %s: %d x %d px, %.0f KB\n", out, info$width, info$height,
            file.size(out) / 1024))
