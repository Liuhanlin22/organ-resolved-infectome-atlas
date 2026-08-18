#!/usr/bin/env Rscript

# Figure 6A: deterministic ring layout with order blocks and straight edges.
required <- c("ggplot2", "ggrepel", "scales", "jsonlite")
missing <- required[!vapply(required, requireNamespace, logical(1), quietly = TRUE)]
if (length(missing) > 0) stop("Missing R package(s): ", paste(missing, collapse = ", "))
suppressPackageStartupMessages({library(ggplot2); library(ggrepel); library(scales)})

arg <- grep("^--file=", commandArgs(FALSE), value = TRUE)
script_dir <- if (length(arg) == 1) dirname(normalizePath(sub("^--file=", "", arg), winslash = "/")) else normalizePath(getwd(), winslash = "/")
read_tsv <- function(path) read.delim(path, check.names = FALSE, stringsAsFactors = FALSE,
                                      fileEncoding = "UTF-8-BOM")
nodes <- read_tsv(file.path(script_dir, "Figure6A_node_table.tsv"))
edges <- read_tsv(file.path(script_dir, "Figure6A_edge_table.tsv"))
edges$RPM <- as.numeric(edges$RPM); edges$log10_RPM_plus1 <- as.numeric(edges$log10_RPM_plus1)
original_path <- file.path(script_dir, "Figure6A_node_positions_original.tsv")
if (!file.exists(original_path)) original_path <- file.path(script_dir, "Figure6A_node_positions.tsv")
original_positions <- read_tsv(original_path)
stopifnot(sum(nodes$node_type == "Host species") == 42, sum(nodes$node_type != "Host species") == 75)
stopifnot(nrow(edges) == 118, all(edges$RPM > 0))

node_colors <- c("Host species" = "#c4e4f6", "Virus" = "#ef186c", "Bacteria" = "#d4b9da", "Eukaryota" = "#a1d99b")
order_colour_table <- read_tsv(file.path(script_dir, "Figure6A_order_colors.tsv"))
order_colors <- setNames(order_colour_table$colour, order_colour_table$host_order)
hosts <- nodes[nodes$node_type == "Host species", , drop = FALSE]
pathogens <- nodes[nodes$node_type != "Host species", , drop = FALSE]
outer <- pathogens[pathogens$node_zone == "outer", , drop = FALSE]
inner <- pathogens[pathogens$node_zone == "inner", , drop = FALSE]
all_orders <- unique(hosts$host_order)
if (!setequal(all_orders, names(order_colors))) stop("Host-order colour mapping is incomplete")

# User-specified host-order background colours. Named assignment keeps the
# mapping stable when the optimiser changes the circular order.
order_colors <- c(
  "Anguilliformes" = "#e0f3db",
  "Carangiformes" = "#ccebc5",
  "Mugiliformes" = "#fff7ec",
  "Aulopiformes" = "#fee8c8",
  "Lutjaniformes" = "#fdd49e",
  "Pleuronectiformes" = "#fff7fb",
  "Gadiformes" = "#ece7f2",
  "Clupeiformes" = "#d0d1e6",
  "Spariformes" = "#fde0dd",
  "Siluriformes" = "#fcc5c0",
  "Chaetodontiformes" = "#f9f9ed",
  "Perciformes" = "#fcfbde",
  "Carcharhiniformes" = "#fffacc",
  "Gobiiformes" = "#d9f0a3",
  "Tetraodontiformes" = "#f2e0df",
  "Eupercaria incertae sedis" = "#e5f5f9",
  "Scombriformes" = "#ccece6",
  "Myctophiformes" = "#e0ecf4",
  "Carangaria incertae sedis" = "#bfd3e6",
  "Myliobatiformes" = "#f7fcf0"
)
stopifnot(setequal(names(order_colors), all_orders), !anyDuplicated(names(order_colors)))
order_to_hosts <- lapply(all_orders, function(ord) hosts$node_id[hosts$host_order == ord]); names(order_to_hosts) <- all_orders
pathogen_edges <- split(seq_len(nrow(edges)), edges$target)
node_zone <- setNames(nodes$node_zone, nodes$node_id); node_label <- setNames(nodes$label, nodes$node_id)

# Independent edge pairs only: shared endpoints are legal junctions, not crossings.
edge_i <- integer(0); edge_j <- integer(0)
for (i in seq_len(nrow(edges) - 1L)) for (j in (i + 1L):nrow(edges)) {
  if (!length(intersect(c(edges$source[i], edges$target[i]), c(edges$source[j], edges$target[j])))) {
    edge_i <- c(edge_i, i); edge_j <- c(edge_j, j)
  }
}
crossing_breakdown <- function(pos) {
  orient_vec <- function(ax, ay, bx, by, cx, cy) (bx - ax) * (cy - ay) - (by - ay) * (cx - ax)
  a <- pos[edges$source[edge_i], , drop = FALSE]; b <- pos[edges$target[edge_i], , drop = FALSE]
  c <- pos[edges$source[edge_j], , drop = FALSE]; d <- pos[edges$target[edge_j], , drop = FALSE]
  ab_c <- orient_vec(a[, 1], a[, 2], b[, 1], b[, 2], c[, 1], c[, 2])
  ab_d <- orient_vec(a[, 1], a[, 2], b[, 1], b[, 2], d[, 1], d[, 2])
  cd_a <- orient_vec(c[, 1], c[, 2], d[, 1], d[, 2], a[, 1], a[, 2])
  cd_b <- orient_vec(c[, 1], c[, 2], d[, 1], d[, 2], b[, 1], b[, 2])
  proper <- abs(ab_c) >= 1e-10 & abs(ab_d) >= 1e-10 & abs(cd_a) >= 1e-10 & abs(cd_b) >= 1e-10 &
    ((ab_c > 0) != (ab_d > 0)) & ((cd_a > 0) != (cd_b > 0))
  z1 <- node_zone[edges$target[edge_i]]; z2 <- node_zone[edges$target[edge_j]]
  c(total = sum(proper), outer_outer = sum(proper & z1 == "outer" & z2 == "outer"),
    inner_inner = sum(proper & z1 == "inner" & z2 == "inner"), inner_outer = sum(proper & z1 != z2))
}
edge_length_sum <- function(pos) sum(sqrt((pos[edges$source, 1] - pos[edges$target, 1])^2 + (pos[edges$source, 2] - pos[edges$target, 2])^2))
circular_mean <- function(theta, weights = rep(1, length(theta))) {
  if (!length(theta)) return(0); weights <- weights / sum(weights)
  atan2(sum(weights * sin(theta)), sum(weights * cos(theta))) %% (2 * pi)
}
clamp_radius <- function(x, y, maximum) {
  rr <- sqrt(x^2 + y^2); if (rr > maximum && rr > 0) c(x * maximum / rr, y * maximum / rr) else c(x, y)
}

# With 42 hosts and 20 order boundaries, d=2*pi/62. Adding one extra d
# after every block makes same-order gaps d and between-order gaps exactly 2d.
host_positions_for_sequence <- function(order_sequence, order_to_hosts) {
  step <- 2 * pi / (sum(lengths(order_to_hosts[order_sequence])) + length(order_sequence))
  cursor <- pi / 2; out <- data.frame(); host_sequence <- character(0)
  for (ord in order_sequence) {
    for (id in order_to_hosts[[ord]]) {
      theta_unwrapped <- cursor
      theta <- theta_unwrapped %% (2 * pi)
      out <- rbind(out, data.frame(node_id = id, label = node_label[id], node_type = "Host species",
                                   node_zone = "host_ring", host_order = ord, theta = theta,
                                   theta_unwrapped = theta_unwrapped,
                                   x = 0.98 * cos(theta), y = 0.98 * sin(theta)))
      host_sequence <- c(host_sequence, id); cursor <- cursor + step
    }
    cursor <- cursor + step
  }
  rownames(out) <- NULL; list(data = out, sequence = host_sequence)
}

make_positions <- function(order_sequence, order_to_hosts) {
  hp <- host_positions_for_sequence(order_sequence, order_to_hosts)
  pos <- matrix(NA_real_, nrow = nrow(nodes), ncol = 2, dimnames = list(nodes$node_id, c("x", "y")))
  pos[hp$data$node_id, ] <- as.matrix(hp$data[, c("x", "y")])
  host_rank <- setNames(seq_along(hp$sequence), hp$sequence); host_theta <- setNames(hp$data$theta, hp$data$node_id)
  outer_by_order <- lapply(order_sequence, function(ord) {
    ids <- outer$node_id[vapply(outer$node_id, function(pid) unique(edges$host_order[pathogen_edges[[pid]]])[1] == ord, logical(1))]
    if (!length(ids)) return(character(0))
    score <- vapply(ids, function(pid) mean(host_rank[edges$source[pathogen_edges[[pid]]]]), numeric(1))
    ids[order(score, node_label[ids])]
  })
  outer_sequence <- unlist(outer_by_order, use.names = FALSE)
  base <- 2 * pi * (seq_along(outer_sequence) - 1) / length(outer_sequence)
  desired <- vapply(outer_sequence, function(pid) circular_mean(host_theta[edges$source[pathogen_edges[[pid]]]]), numeric(1))
  theta <- (base + circular_mean((desired - base) %% (2 * pi))) %% (2 * pi)
  stopifnot(is.matrix(pos), is.character(outer_sequence), length(theta) == length(outer_sequence))
  pos[outer_sequence, ] <- cbind(1.125 * cos(theta), 1.125 * sin(theta))
  for (pid in inner$node_id) {
    ii <- pathogen_edges[[pid]]; w <- 0.5 + edges$log10_RPM_plus1[ii]
    xy0 <- colSums(pos[edges$source[ii], , drop = FALSE] * w) / sum(w) * 0.72
    if (sqrt(sum(xy0^2)) < 0.07) {a <- circular_mean(host_theta[edges$source[ii]], w); xy0 <- 0.09 * c(cos(a), sin(a))}
    pos[pid, ] <- clamp_radius(xy0[1], xy0[2], 0.66)
  }
  list(pos = pos, hosts = hp$data, host_sequence = hp$sequence, outer_sequence = outer_sequence)
}
layout_score <- function(pos) crossing_breakdown(pos)["total"] * 100 + edge_length_sum(pos)

affinity_seed_sequence <- function() {
  affinity <- matrix(0, length(all_orders), length(all_orders), dimnames = list(all_orders, all_orders))
  totals <- setNames(numeric(length(all_orders)), all_orders)
  for (pid in inner$node_id) {
    linked <- unique(edges$host_order[pathogen_edges[[pid]]]); if (length(linked) < 2) next
    contribution <- 1 / max(1, length(linked) - 1)
    for (aa in seq_len(length(linked) - 1L)) for (bb in (aa + 1L):length(linked)) {
      affinity[linked[aa], linked[bb]] <- affinity[linked[aa], linked[bb]] + contribution
      affinity[linked[bb], linked[aa]] <- affinity[linked[bb], linked[aa]] + contribution
      totals[linked[aa]] <- totals[linked[aa]] + contribution; totals[linked[bb]] <- totals[linked[bb]] + contribution
    }
  }
  sequence <- all_orders[which.max(totals)]; unvisited <- setdiff(all_orders, sequence)
  while (length(unvisited)) {
    left_w <- affinity[sequence[1], unvisited]; right_w <- affinity[tail(sequence, 1), unvisited]
    use_left <- left_w >= right_w; best <- order(pmax(left_w, right_w), totals[unvisited], unvisited, decreasing = TRUE)[1]
    chosen <- unvisited[best]; if (use_left[best]) sequence <- c(chosen, sequence) else sequence <- c(sequence, chosen)
    unvisited <- setdiff(unvisited, chosen)
  }
  sequence
}

optimise_order_blocks <- function() {
  set.seed(20260815); starts <- list(all_orders, affinity_seed_sequence()); best_sequence <- all_orders; best_score <- Inf
  for (start in starts) {
    current <- start; current_score <- layout_score(make_positions(current, order_to_hosts)$pos)
    local_best <- current; local_score <- current_score
    for (iter in seq_len(2400)) {
      candidate <- current; ix <- sort(sample(2:length(candidate), 2))
      if (runif(1) < 0.55) candidate[ix] <- rev(candidate[ix]) else candidate[ix[1]:ix[2]] <- rev(candidate[ix[1]:ix[2]])
      score <- layout_score(make_positions(candidate, order_to_hosts)$pos); temp <- max(0.08, 8 * (1 - iter / 2400))
      if (score < current_score || runif(1) < exp(min(0, (current_score - score) / temp))) {current <- candidate; current_score <- score}
      if (score < local_score) {local_best <- candidate; local_score <- score}
    }
    if (local_score < best_score) {best_sequence <- local_best; best_score <- local_score}
  }
  best_sequence
}

optimise_hosts_within_orders <- function(order_sequence) {
  set.seed(20260816); current_map <- order_to_hosts; current_score <- layout_score(make_positions(order_sequence, current_map)$pos)
  for (ord in order_sequence) {
    ids <- current_map[[ord]]; if (length(ids) < 2) next; attempts <- 120 * length(ids)
    for (iter in seq_len(attempts)) {
      candidate_ids <- ids; ix <- sample(seq_along(ids), 2); candidate_ids[ix] <- candidate_ids[rev(ix)]
      candidate_map <- current_map; candidate_map[[ord]] <- candidate_ids
      score <- layout_score(make_positions(order_sequence, candidate_map)$pos); temp <- max(0.05, 3 * (1 - iter / attempts))
      if (score < current_score || runif(1) < exp(min(0, (current_score - score) / temp))) {
        ids <- candidate_ids; current_map <- candidate_map; current_score <- score
      }
    }
    current_map[[ord]] <- ids
  }
  current_map
}

inner_collision_penalty <- function(pos) {
  xy <- pos[inner$node_id, , drop = FALSE]; d <- as.matrix(dist(xy))
  sum((pmax(0, 0.13 - d)[upper.tri(d)])^2) * 500
}
final_score <- function(pos) layout_score(pos) + inner_collision_penalty(pos)

optimise_inner_nodes <- function(pos) {
  set.seed(20260817); anchors <- pos[inner$node_id, , drop = FALSE]; relaxed <- anchors
  for (iter in seq_len(260)) {
    step <- 0.025 * (1 - iter / 300); updates <- relaxed
    for (pid in inner$node_id) {
      xy <- relaxed[pid, ]; delta <- (anchors[pid, ] - xy) * 0.12
      for (other_id in setdiff(inner$node_id, pid)) {
        v <- xy - relaxed[other_id, ]; ds <- sum(v^2) + 1e-5
        if (ds < 0.12) delta <- delta + v * (0.0025 / ds)
      }
      updates[pid, ] <- clamp_radius(xy[1] + step * delta[1], xy[2] + step * delta[2], 0.66)
    }
    relaxed <- updates
  }
  pos[inner$node_id, ] <- relaxed
  current <- pos; current_score <- final_score(current); best <- current; best_score <- current_score
  for (iter in seq_len(5200)) {
    pid <- sample(inner$node_id, 1); candidate <- current; scale <- 0.085 * (1 - iter / 6200) + 0.012
    candidate[pid, ] <- clamp_radius(current[pid, 1] + rnorm(1, 0, scale), current[pid, 2] + rnorm(1, 0, scale), 0.66)
    score <- final_score(candidate); temp <- max(0.08, 7 * (1 - iter / 5200))
    if (score < current_score || runif(1) < exp(min(0, (current_score - score) / temp))) {current <- candidate; current_score <- score}
    if (score < best_score) {best <- candidate; best_score <- score}
  }
  best
}

order_sequence <- optimise_order_blocks()
order_to_hosts <- optimise_hosts_within_orders(order_sequence)
layout <- make_positions(order_sequence, order_to_hosts)
positions <- optimise_inner_nodes(layout$pos)

host_pos <- layout$hosts[match(unlist(order_to_hosts[order_sequence], use.names = FALSE), layout$hosts$node_id), ]
host_pos$theta <- atan2(positions[host_pos$node_id, 2], positions[host_pos$node_id, 1]) %% (2 * pi)
host_pos$x <- positions[host_pos$node_id, 1]; host_pos$y <- positions[host_pos$node_id, 2]
outer_pos <- outer[match(layout$outer_sequence, outer$node_id), , drop = FALSE]
outer_pos$x <- positions[outer_pos$node_id, 1]; outer_pos$y <- positions[outer_pos$node_id, 2]
inner_pos <- inner
inner_pos$x <- positions[inner_pos$node_id, 1]; inner_pos$y <- positions[inner_pos$node_id, 2]

# Smooth capsule-shaped order bands only cover their host nodes. Using
# unwrapped angles prevents the 357-20 degree Eupercaria block from filling
# nearly the whole ring. Rounded end caps replace the previous sharp corners.
host_radius <- 0.98
band_corner_radius <- 0.07
band_inner <- host_radius - band_corner_radius
band_outer <- host_radius + band_corner_radius
band_cap_angle <- asin(band_corner_radius / host_radius)
band_ranges <- do.call(rbind, lapply(seq_along(order_sequence), function(i) {
  ord <- order_sequence[i]
  h <- host_pos[host_pos$host_order == ord, , drop = FALSE]
  start <- min(h$theta_unwrapped) - band_cap_angle
  end <- max(h$theta_unwrapped) + band_cap_angle
  data.frame(
    order_index = i, host_order = ord, host_count = nrow(h),
    start_degrees_unwrapped = start * 180 / pi,
    end_degrees_unwrapped = end * 180 / pi,
    start_degrees_normalized = (start %% (2 * pi)) * 180 / pi,
    end_degrees_normalized = (end %% (2 * pi)) * 180 / pi,
    crosses_zero_degrees = floor(start / (2 * pi)) != floor(end / (2 * pi)),
    fill_colour = unname(order_colors[ord]),
    stringsAsFactors = FALSE
  )
}))
make_rounded_band <- function(ord) {
  h <- host_pos[host_pos$host_order == ord, , drop = FALSE]
  start <- min(h$theta_unwrapped)
  end <- max(h$theta_unwrapped)
  if (nrow(h) == 1) {
    phi <- seq(0, 2 * pi, length.out = 121)
    return(data.frame(
      host_order = ord,
      x = host_radius * cos(start) + band_corner_radius * cos(phi),
      y = host_radius * sin(start) + band_corner_radius * sin(phi)
    ))
  }
  arc_points <- max(40, ceiling((end - start) * 180 / pi * 2))
  outer_theta <- seq(start, end, length.out = arc_points)
  inner_theta <- seq(end, start, length.out = arc_points)
  end_alpha <- seq(0, pi, length.out = 41)
  start_alpha <- seq(pi, 2 * pi, length.out = 41)
  cap_xy <- function(theta, alpha) {
    cbind(
      host_radius * cos(theta) + band_corner_radius *
        (cos(alpha) * cos(theta) - sin(alpha) * sin(theta)),
      host_radius * sin(theta) + band_corner_radius *
        (cos(alpha) * sin(theta) + sin(alpha) * cos(theta))
    )
  }
  end_cap <- cap_xy(end, end_alpha)
  start_cap <- cap_xy(start, start_alpha)
  data.frame(
    host_order = ord,
    x = c(band_outer * cos(outer_theta), end_cap[, 1],
          band_inner * cos(inner_theta), start_cap[, 1]),
    y = c(band_outer * sin(outer_theta), end_cap[, 2],
          band_inner * sin(inner_theta), start_cap[, 2])
  )
}
band_data <- do.call(rbind, lapply(order_sequence, make_rounded_band))

radial_label <- function(theta, radius) {
  deg <- theta * 180 / pi; left <- cos(theta) < 0
  data.frame(label_x = radius * cos(theta), label_y = radius * sin(theta),
             label_angle = ifelse(left, deg + 180, deg), label_hjust = ifelse(left, 1, 0))
}
outer_lab <- cbind(outer_pos, radial_label(atan2(outer_pos$y, outer_pos$x), 1.125))
pt_to_mm <- 25.4 / 72.27
host_label_size <- 1.55 + 2 * pt_to_mm
inner_label_size <- 2.58 + 2 * pt_to_mm
outer_label_size <- 2.24 + 1 * pt_to_mm
node_pos <- nodes
node_pos$x <- positions[node_pos$node_id, 1]; node_pos$y <- positions[node_pos$node_id, 2]
node_pos$angle_degrees <- (atan2(node_pos$y, node_pos$x) %% (2 * pi)) * 180 / pi
node_pos$radius <- sqrt(node_pos$x^2 + node_pos$y^2)
edges_plot <- data.frame(
  x_host = positions[edges$source, 1], y_host = positions[edges$source, 2],
  x_pathogen = positions[edges$target, 1], y_pathogen = positions[edges$target, 2],
  log10_RPM_plus1 = edges$log10_RPM_plus1
)

edge_breaks <- seq(min(edges$log10_RPM_plus1), max(edges$log10_RPM_plus1), length.out = 4)
write.table(node_pos, file.path(script_dir, "Figure6A_node_positions.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(order_index = seq_along(order_sequence), host_order = order_sequence,
                       colour = unname(order_colors[order_sequence])),
            file.path(script_dir, "Figure6A_order_colors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(band_ranges, file.path(script_dir, "Figure6A_order_band_ranges.tsv"),
            sep = "\t", quote = FALSE, row.names = FALSE)
write.table(data.frame(log10_RPM_plus1 = edge_breaks,
                       plotted_linewidth = scales::rescale(edge_breaks, to = c(0.18, 2.25))),
            file.path(script_dir, "Figure6A_edge_legend.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

p <- ggplot() +
  geom_polygon(data = band_data, aes(x, y, group = host_order, fill = host_order), colour = NA, alpha = 1) +
  geom_segment(data = edges_plot,
               aes(x = x_host, y = y_host, xend = x_pathogen, yend = y_pathogen,
                   linewidth = log10_RPM_plus1), colour = "#737373", alpha = 0.43, lineend = "round") +
  geom_point(data = node_pos, aes(x, y, colour = node_type, size = node_type), stroke = 0) +
  geom_text(data = host_pos, aes(x, y, label = label), family = "Arial", fontface = "italic",
            size = host_label_size, colour = "#202020", lineheight = 0.82) +
  geom_text(data = outer_lab,
            aes(label_x, label_y, label = label, angle = label_angle, hjust = label_hjust),
            family = "Arial", fontface = "italic", size = outer_label_size, colour = "#111111") +
  geom_label_repel(data = inner_pos, aes(x, y, label = label), family = "Arial", fontface = "italic",
                   size = inner_label_size, colour = "#111111", fill = alpha("white", 0.88),
                   label.size = 0.35, label.r = grid::unit(0.08, "lines"),
                   label.padding = grid::unit(0.10, "lines"), seed = 20260815, box.padding = 0.30,
                   point.padding = 0.12, min.segment.length = 0, segment.colour = "#8f8f8f",
                   segment.size = 0.26, force = 8, force_pull = 0.30, max.iter = 50000,
                   max.time = 20, max.overlaps = Inf, xlim = c(-0.86, 0.86), ylim = c(-0.86, 0.86)) +
  annotate("text", x = -1.52, y = 1.50, label = "a", family = "Arial", fontface = "bold", size = 8.7) +
  scale_fill_manual(values = order_colors, breaks = order_sequence, name = "Host order") +
  scale_colour_manual(values = node_colors, breaks = names(node_colors), name = "Node") +
  scale_size_manual(values = c("Host species" = 12.0, "Virus" = 3.15, "Bacteria" = 3.15, "Eukaryota" = 3.15),
                    breaks = names(node_colors), name = "Node") +
  scale_linewidth_continuous(range = c(0.18, 2.25), breaks = edge_breaks,
                             labels = number_format(accuracy = 0.1), name = expression(log[10] * "(RPM+1)")) +
  guides(colour = guide_legend(order = 1, override.aes = list(size = c(12.0, 3.15, 3.15, 3.15))),
         size = "none", fill = guide_legend(order = 2, ncol = 1, byrow = TRUE, override.aes = list(alpha = 1)),
         linewidth = guide_legend(order = 3)) +
  coord_fixed(xlim = c(-1.56, 1.56), ylim = c(-1.56, 1.56), clip = "off") +
  theme_void(base_family = "Arial") +
  theme(legend.position = "right", legend.box = "vertical",
        legend.title = element_text(size = 14, face = "bold"), legend.text = element_text(size = 11.2),
        legend.key.height = grid::unit(0.44, "cm"), legend.key.width = grid::unit(0.78, "cm"),
        plot.margin = margin(18, 18, 18, 18))

ggsave(file.path(script_dir, "Figure6A.pdf"), p, width = 18.5, height = 15.2, units = "in", device = cairo_pdf, bg = "white")
ggsave(file.path(script_dir, "Figure6A.png"), p, width = 18.5, height = 15.2, units = "in", dpi = 600, bg = "white")

before_pos <- matrix(NA_real_, nrow = nrow(nodes), ncol = 2, dimnames = list(nodes$node_id, c("x", "y")))
before_pos[original_positions$node_id, ] <- as.matrix(original_positions[, c("x", "y")])
before <- crossing_breakdown(before_pos); after <- crossing_breakdown(positions)
host_r <- node_pos$radius[node_pos$node_type == "Host species"]
outer_r <- node_pos$radius[node_pos$node_zone == "outer"]
inner_dist <- as.matrix(dist(inner_pos[, c("x", "y")]))
metrics <- list(
  random_seed = 20260815L, nodes = nrow(nodes), host_nodes = nrow(hosts), pathogen_nodes = nrow(pathogens),
  edges = nrow(edges), host_orders = length(order_sequence), host_radius_mean = mean(host_r),
  host_radius_cv = sd(host_r) / mean(host_r), outer_radius_mean = mean(outer_r),
  outer_radius_cv = sd(outer_r) / mean(outer_r), outer_radius_target = 1.125,
  within_order_gap_radians = 2 * pi / 62, between_order_gap_radians = 4 * pi / 62,
  host_marker_size = 12.0, outer_marker_size = 3.15,
  host_label_size_mm = host_label_size, inner_label_size_mm = inner_label_size,
  outer_label_size_mm = outer_label_size, band_inner = band_inner,
  band_outer = band_outer, band_corner_radius = band_corner_radius,
  host_order_legend_columns = 1L, legend_title_size_pt = 14,
  legend_text_size_pt = 11.2,
  order_band_white_gap_radians = 4 * pi / 62 - 2 * band_cap_angle,
  order_band_white_gap_degrees = (4 * pi / 62 - 2 * band_cap_angle) * 180 / pi,
  order_band_count = nrow(band_ranges), cross_zero_order_count = sum(band_ranges$crosses_zero_degrees),
  order_sequence = order_sequence,
  crossings_before = as.list(before), crossings_after = as.list(after),
  crossing_reduction_percent = round(100 * (before["total"] - after["total"]) / max(1, before["total"]), 2),
  min_inner_node_distance = min(inner_dist[upper.tri(inner_dist)])
)
jsonlite::write_json(metrics, file.path(script_dir, "Figure6A_layout_metrics.json"), pretty = TRUE, auto_unbox = TRUE)
scalar_metrics <- metrics[!vapply(metrics, is.list, logical(1))]
write.table(data.frame(metric = names(scalar_metrics),
                       value = vapply(scalar_metrics, function(x) paste(x, collapse = ";"), character(1))),
            file.path(script_dir, "Figure6A_layout_metrics.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
