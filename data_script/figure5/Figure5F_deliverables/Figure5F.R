#!/usr/bin/env Rscript

# Figure 5F — descriptive observed pathogen co-infection network
# Analysis unit: 112 fish with complete eight-organ sampling.
# Network edges show observed co-detection counts and are descriptive.

required_packages <- c("openxlsx", "igraph", "ggplot2", "ggraph", "ggrepel", "scales")
missing_packages <- required_packages[
  !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0) {
  stop("Missing R package(s): ", paste(missing_packages, collapse = ", "))
}

suppressPackageStartupMessages({
  library(openxlsx)
  library(igraph)
  library(ggplot2)
  library(ggraph)
  library(ggrepel)
  library(scales)
})

script_arg <- grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)
script_dir <- if (length(script_arg) == 1) {
  dirname(normalizePath(sub("^--file=", "", script_arg), winslash = "/"))
} else {
  normalizePath(getwd(), winslash = "/")
}

input_file <- file.path(script_dir, "Figure5F_input.xlsx")
if (!file.exists(input_file)) stop("Input file not found: ", input_file)

rpm_df <- read.xlsx(input_file, sheet = "host_by_pathogen_RPM", check.names = FALSE)
annotations <- read.xlsx(input_file, sheet = "pathogen_annotations", check.names = FALSE)

if (!"host_id" %in% names(rpm_df)) stop("Missing host_id column")
matrix_columns <- setdiff(names(rpm_df), "host_id")
pathogen_names <- annotations$Species.name
if (nrow(rpm_df) != 112 || length(matrix_columns) != 75 || length(pathogen_names) != 75) {
  stop("Expected 112 fish x 75 pathogens; observed ", nrow(rpm_df), " x ", length(pathogen_names))
}
if (anyDuplicated(rpm_df$host_id) || anyDuplicated(pathogen_names)) {
  stop("Duplicated host or pathogen identifiers detected")
}
rpm <- as.matrix(rpm_df[, matrix_columns, drop = FALSE])
colnames(rpm) <- pathogen_names
storage.mode(rpm) <- "numeric"
if (anyNA(rpm) || any(rpm < 0)) stop("RPM matrix contains missing or negative values")

present <- rpm > 0
host_richness <- rowSums(present)
infected_hosts <- colSums(present)
total_detections <- vapply(
  seq_along(pathogen_names),
  function(j) sum(host_richness[present[, j]]),
  numeric(1)
)
average_pathogens <- ifelse(infected_hosts > 0, total_detections / infected_hosts, 0)

ann_idx <- match(pathogen_names, annotations$Species.name)
node_all <- data.frame(
  name = pathogen_names,
  classification = annotations$Classification[ann_idx],
  novelty = annotations$Novelty[ann_idx],
  pathogen_characteristics = annotations$Pathogen.characteristics[ann_idx],
  kingdom = annotations$Kingdom[ann_idx],
  family = annotations$Family[ann_idx],
  infected_hosts = as.integer(infected_hosts),
  total_pathogen_detections = as.integer(total_detections),
  average_pathogens_per_infected_host = as.numeric(average_pathogens),
  stringsAsFactors = FALSE,
  check.names = FALSE
)

detected <- pathogen_names[infected_hosts > 0]
pair_index <- combn(detected, 2, simplify = FALSE)
pairwise <- do.call(
  rbind,
  lapply(pair_index, function(pair) {
    x <- present[, pair[1]]
    y <- present[, pair[2]]
    both <- sum(x & y)
    x_only <- sum(x & !y)
    y_only <- sum(!x & y)
    neither <- sum(!x & !y)
    tab <- matrix(c(both, x_only, y_only, neither), nrow = 2, byrow = TRUE)
    ft <- fisher.test(tab, alternative = "greater")
    denom <- sqrt(sum(x) * (length(x) - sum(x)) * sum(y) * (length(y) - sum(y)))
    phi <- if (denom > 0) (both * neither - x_only * y_only) / denom else NA_real_
    data.frame(
      source = pair[1],
      target = pair[2],
      source_infected_hosts = sum(x),
      target_infected_hosts = sum(y),
      co_infected_hosts = both,
      expected_co_infected_hosts = sum(x) * sum(y) / length(x),
      odds_ratio = unname(ft$estimate),
      phi_coefficient = phi,
      fisher_p = ft$p.value,
      stringsAsFactors = FALSE
    )
  })
)
pairwise$fisher_BH_q <- p.adjust(pairwise$fisher_p, method = "BH")
pairwise$significant_positive_BH <- with(
  pairwise,
  fisher_BH_q < 0.05 & co_infected_hosts > expected_co_infected_hosts & odds_ratio > 1
)
pairwise <- pairwise[order(pairwise$fisher_BH_q, pairwise$fisher_p), ]

edges <- pairwise[pairwise$co_infected_hosts >= 1, ]
edges <- edges[order(-edges$co_infected_hosts, edges$source, edges$target), ]
network_names <- unique(c(edges$source, edges$target))
nodes <- node_all[node_all$name %in% network_names, ]

g0 <- graph_from_data_frame(
  edges[, c("source", "target", "co_infected_hosts", "phi_coefficient", "fisher_p", "fisher_BH_q")],
  directed = FALSE,
  vertices = nodes
)

community_membership <- if (ecount(g0) > 0) {
  membership(cluster_fast_greedy(g0, weights = E(g0)$co_infected_hosts))
} else {
  setNames(rep(1L, vcount(g0)), V(g0)$name)
}
nodes$community <- unname(community_membership[nodes$name])

# Circular order: highest infected-host count first, then proceed counter-clockwise.
# Toxocara sp. CTSZ has the highest occurrence and is fixed at the bottom.
nodes <- nodes[
  order(
    -nodes$infected_hosts,
    -nodes$average_pathogens_per_infected_host,
    nodes$name
  ),
]
nodes$plot_order <- seq_len(nrow(nodes))
plot_angle <- -pi / 2 + 2 * pi * (nodes$plot_order - 1) / nrow(nodes)
nodes$plot_x <- cos(plot_angle)
nodes$plot_y <- sin(plot_angle)

g <- graph_from_data_frame(
  edges[, c("source", "target", "co_infected_hosts", "phi_coefficient", "fisher_p", "fisher_BH_q")],
  directed = FALSE,
  vertices = nodes
)

nodes$degree <- degree(g)[nodes$name]
nodes$weighted_degree <- strength(g, weights = E(g)$co_infected_hosts)[nodes$name]
nodes$betweenness <- betweenness(
  g,
  directed = FALSE,
  weights = 1 / E(g)$co_infected_hosts,
  normalized = TRUE
)[nodes$name]
nodes$local_clustering <- transitivity(g, type = "local", isolates = "zero")[nodes$name]
nodes$component <- components(g)$membership[nodes$name]

for (column in names(nodes)) {
  set_vertex_attr(g, column, value = nodes[[column]][match(V(g)$name, nodes$name)])
}

write.csv(node_all, file.path(script_dir, "Figure5F_nodes_all.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(nodes, file.path(script_dir, "Figure5F_nodes_network.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(edges, file.path(script_dir, "Figure5F_edges_observed.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write.csv(pairwise, file.path(script_dir, "Figure5F_pairwise_Fisher_BH.csv"), row.names = FALSE, fileEncoding = "UTF-8")
write_graph(g, file.path(script_dir, "Figure5F_network.graphml"), format = "graphml")

summary_table <- data.frame(
  metric = c(
    "Total fish", "Pathogen-positive fish", "Pathogen-negative fish",
    "Total pathogens", "Detected pathogens", "Network nodes",
    "Observed co-infection edges", "BH-significant positive pairs"
  ),
  value = c(
    nrow(rpm), sum(host_richness > 0), sum(host_richness == 0),
    ncol(rpm), sum(infected_hosts > 0), vcount(g), ecount(g),
    sum(pairwise$significant_positive_BH)
  )
)

wb <- createWorkbook()
header_style <- createStyle(
  fontName = "Arial", fontSize = 10, textDecoration = "bold",
  fgFill = "#D9EAF7", halign = "center", valign = "center", wrapText = TRUE
)
body_style <- createStyle(fontName = "Arial", fontSize = 9, valign = "center")
add_table <- function(sheet_name, data) {
  addWorksheet(wb, sheet_name)
  writeData(wb, sheet_name, data, headerStyle = header_style, withFilter = TRUE)
  addStyle(wb, sheet_name, body_style, rows = 2:(nrow(data) + 1), cols = 1:ncol(data), gridExpand = TRUE)
  freezePane(wb, sheet_name, firstRow = TRUE)
  setColWidths(wb, sheet_name, cols = 1:ncol(data), widths = "auto")
}
add_table("summary", summary_table)
add_table("nodes_all", node_all)
add_table("nodes_network", nodes)
add_table("edges_observed", edges)
add_table("pairwise_Fisher_BH", pairwise)
saveWorkbook(wb, file.path(script_dir, "Figure5F_network_results.xlsx"), overwrite = TRUE)

size_breaks <- sort(unique(c(1, 8, max(nodes$infected_hosts))))
edge_breaks <- sort(unique(edges$co_infected_hosts))
color_max <- max(7, ceiling(max(nodes$average_pathogens_per_infected_host)))

set.seed(20260815)
p <- ggraph(
  g,
  layout = "manual",
  x = nodes$plot_x[match(V(g)$name, nodes$name)],
  y = nodes$plot_y[match(V(g)$name, nodes$name)]
) +
  geom_edge_link(
    aes(width = co_infected_hosts),
    edge_colour = "#777777", alpha = 0.28, lineend = "round",
    show.legend = TRUE
  ) +
  scale_edge_width_continuous(
    name = "Strength of\nco-infection",
    range = c(0.18, 1.35), breaks = edge_breaks
  ) +
  geom_node_point(
    aes(size = infected_hosts, colour = average_pathogens_per_infected_host),
    shape = 16
  ) +
  geom_node_text(
    aes(label = name),
    repel = TRUE, max.overlaps = Inf,
    family = "Arial", fontface = "italic", size = 3.0,
    box.padding = 0.28, point.padding = 0.18,
    segment.colour = "#A6A6A6", segment.size = 0.22,
    min.segment.length = 0
  ) +
  scale_size_continuous(
    name = "Number of\ninfected hosts",
    range = c(2.5, 10), breaks = size_breaks
  ) +
  scale_colour_gradientn(
    name = "Avg. pathogens\nper host",
    colours = c("#FDD49E", "#F768A1", "#AE017E", "#49006A"),
    limits = c(0, color_max),
    breaks = seq(0, color_max, length.out = 5),
    oob = squish
  ) +
  coord_equal(clip = "off", xlim = c(-1.50, 1.50), ylim = c(-1.35, 1.35)) +
  labs(tag = "f") +
  guides(
    size = guide_legend(order = 1, override.aes = list(colour = "black")),
    colour = guide_colourbar(order = 2, barheight = unit(36, "mm")),
    edge_width = guide_legend(order = 3)
  ) +
  theme_void(base_family = "Arial", base_size = 10) +
  theme(
    legend.position = "left",
    legend.box = "vertical",
    legend.title = element_text(face = "bold", size = 10),
    legend.text = element_text(size = 8.5),
    plot.tag = element_text(face = "bold", size = 22),
    plot.tag.position = c(0.01, 0.99),
    plot.margin = margin(10, 18, 10, 10)
  )

ggsave(
  file.path(script_dir, "Figure5F.pdf"), p,
  width = 10.5, height = 9.0, units = "in",
  device = cairo_pdf, bg = "white"
)
ggsave(
  file.path(script_dir, "Figure5F.png"), p,
  width = 10.5, height = 9.0, units = "in",
  dpi = 600, bg = "white"
)
