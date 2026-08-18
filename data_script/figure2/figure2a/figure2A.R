library(readxl)
library(dplyr)
library(networkD3)
library(htmlwidgets)
library(webshot)

# Read the data
data0 <- readxl::read_xlsx("figure2A_data.xlsx")

# Clean and prepare the data
data <- data0 %>%
  group_by(Timeline, `Lineages Distribution`, `Geographic Location`, `P.M.A.s`) %>%
  summarise(Count = n(), .groups = 'drop') %>%
  mutate(`P.M.A.s` = ifelse(`P.M.A.s` == "", NA, `P.M.A.s`))

# Define node orders
timeline_order <- c("DNA virus", "RNA virus", "Bacteria", "Eukaryota")
lineages_order <- c("Known species", "Novel species")
pathogen_order <- c("Pathogen", "Nonpathogen")

# Get unique PMAs (excluding NA)
pmas_order <- sort(unique(na.omit(data$`P.M.A.s`)))

# Create nodes data frame
nodes <- data.frame(
  name = c(
    timeline_order,
    lineages_order,
    pathogen_order,
    pmas_order
  )
)

# Create links data
links <- data %>%
  # First layer: Timeline to Lineages Distribution
  select(source = Timeline, target = `Lineages Distribution`, value = Count) %>%
  group_by(source, target) %>%
  summarise(value = sum(value), .groups = 'drop') %>%
  # Second layer: Lineages Distribution to Pathogen/Nonpathogen
  bind_rows(
    data %>%
      select(source = `Lineages Distribution`, target = `Geographic Location`, value = Count) %>%
      group_by(source, target) %>%
      summarise(value = sum(value), .groups = 'drop')
  ) %>%
  # Third layer: Pathogen/Nonpathogen to P.M.A.s (only for Pathogens)
  bind_rows(
    data %>%
      filter(`Geographic Location` == "Pathogen") %>%
      select(source = `Geographic Location`, target = `P.M.A.s`, value = Count) %>%
      filter(!is.na(target)) %>%
      group_by(source, target) %>%
      summarise(value = sum(value), .groups = 'drop')
  )

# Convert source/target to node indices (0-based for networkD3)
links$IDsource <- match(links$source, nodes$name) - 1
links$IDtarget <- match(links$target, nodes$name) - 1

# Custom color palette
customColors <- c(
  "#abdaec", "#fa9fb5", "#f49e39", "#c7e9c0",  # Timeline colors
  "#F0C2A2", "#B5A1E3",                        # Lineages Distribution colors
  "#CB5148", "#6A9BCB",                        # Pathogen/Nonpathogen colors
  "#67001f", "#ce1256", "#e7298a", "#df65b0",  # PMA colors
  "#c994c7", "#d4b9da", "#c5b0d5", "#08306b",
  "#08519c", "#4292c6", "#6baed6", "#9ecae1",
  "#00441b", "#006d2c", "#238b45", "#74c476",
  "#a1d99b", "#67000d", "#cb181d", "#ef3b2c",
  "#fb6a4a", "#fcbba1", "#3f007d", "#6a51a3",
  "#807dba", "#9e9ac8", "#cc4c02", "#fe9929"
)

# Create color scale
colourScale <- sprintf('d3.scaleOrdinal().range(["%s"])', 
                       paste(customColors, collapse = '", "'))

# Create Sankey plot
p <- sankeyNetwork(
  Links = links, 
  Nodes = nodes, 
  Source = "IDsource", 
  Target = "IDtarget", 
  Value = "value", 
  NodeID = "name", 
  units = "Count", 
  fontSize = 18, 
  nodeWidth = 20,
  height = 600,
  width = 2600,
  nodePadding = 10,
  iterations = 32,
  sinksRight = FALSE,
  colourScale = colourScale
)

# Add custom JavaScript for gradient links
p <- htmlwidgets::onRender(
  p,
  '
  function(el, x) {
    // Helper function to create valid ID
    function createValidID(name) {
      return name.replace(/[\\s&]+/g, "_").replace(/[^\\w-]/g, "");
    }
    
    // Select the SVG
    var svg = d3.select(el).select("svg");
    
    // Wait for the layout to complete
    setTimeout(function() {
      // Apply gradients to all links
      svg.selectAll(".link").each(function(d) {
        // Create unique gradient ID
        var gradientID = "gradient-" + createValidID(d.source.name) + "-" + createValidID(d.target.name);
        
        // Get the source and target colors
        var sourceColor = d3.select(el).selectAll(".node").filter(function(node) { 
          return node.name === d.source.name; 
        }).select("rect").style("fill");
        
        var targetColor = d3.select(el).selectAll(".node").filter(function(node) { 
          return node.name === d.target.name; 
        }).select("rect").style("fill");
        
        // Create gradient definition
        var defs = svg.append("defs");
        
        var gradient = defs.append("linearGradient")
          .attr("id", gradientID)
          .attr("gradientUnits", "userSpaceOnUse")
          .attr("x1", d.source.x + d.source.dx / 2)
          .attr("y1", d.source.y + d.source.dy / 2)
          .attr("x2", d.target.x + d.target.dx / 2)
          .attr("y2", d.target.y + d.target.dy / 2);
        
        gradient.append("stop")
          .attr("offset", "0%")
          .attr("stop-color", sourceColor);
        
        gradient.append("stop")
          .attr("offset", "100%")
          .attr("stop-color", targetColor);
        
        // Apply gradient to the link
        d3.select(this).style("stroke", "url(#" + gradientID + ")");
      });
    }, 500); // Delay to ensure layout is complete
  }
  '
)

# Save the visualization
saveNetwork(p, "sankey1.html")  # Save as HTML
# webshot("sankey1.html", "sankey1.pdf", vwidth = 1200, vheight = 700)

# Take a screenshot
# tryCatch({
#   webshot("sankey1.html", "sankey1.png", vwidth = 1200, vheight = 700)
# }, error = function(e) {
# })
