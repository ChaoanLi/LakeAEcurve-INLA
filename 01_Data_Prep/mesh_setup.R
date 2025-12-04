################################################################################
# SPDE Mesh Construction Module
# Authors: Chaoan Li, Yinuo Zhu | STAT 647, Texas A&M University
################################################################################

#' Build INLA SPDE mesh over lake region
#' @param data_list List from load_and_prep_data()
#' @param max_edge Maximum triangle edge length (inner, outer) vector
#' @param cutoff Minimum node spacing
#' @param offset Boundary extension distance (inner, outer) vector
#' @return mesh object
build_mesh <- function(data_list,
                       max_edge = c(100, 500),
                       cutoff = 50,
                       offset = c(100, 500)) {
  
  library(INLA)
  library(sf)
  
  cat("====================================\n")
  cat("Building SPDE mesh...\n")
  cat("====================================\n")
  
  lake_boundary <- data_list$lake_boundary
  
  # 1. Extract boundary coordinates
  cat("\n1. Extracting boundary coordinates...\n")
  
  boundary_coords <- st_coordinates(lake_boundary)[, 1:2]
  
  # Remove duplicate first/last point if present
  if (all(boundary_coords[1, ] == boundary_coords[nrow(boundary_coords), ])) {
    boundary_coords <- boundary_coords[-nrow(boundary_coords), ]
  }
  
  cat(sprintf("   Boundary points: %d\n", nrow(boundary_coords)))
  
  # 2. Build mesh
  cat("\n2. Building mesh...\n")
  cat(sprintf("   max.edge: c(%.1f, %.1f)\n", max_edge[1], max_edge[2]))
  cat(sprintf("   cutoff: %.1f\n", cutoff))
  cat(sprintf("   offset: c(%.1f, %.1f)\n", offset[1], offset[2]))
  
  mesh <- inla.mesh.2d(
    loc.domain = boundary_coords,
    max.edge = max_edge,
    cutoff = cutoff,
    offset = offset
  )
  
  cat(sprintf("\nMesh built successfully.\n"))
  cat(sprintf("   Number of vertices: %d\n", mesh$n))
  cat(sprintf("   Number of triangles: %d\n", nrow(mesh$graph$tv)))
  
  return(mesh)
}


#' Visualize the mesh
#' @param mesh INLA mesh object
#' @param data_list List from load_and_prep_data()
#' @param save_path Path to save plot (optional)
#' @return ggplot object
plot_mesh <- function(mesh, data_list, save_path = NULL) {
  
  library(ggplot2)
  library(sf)
  
  cat("Plotting mesh...\n")
  
  mesh_coords <- mesh$loc[, 1:2]
  mesh_df <- data.frame(x = mesh_coords[, 1], y = mesh_coords[, 2])
  
  triangles <- mesh$graph$tv
  triangle_list <- lapply(1:nrow(triangles), function(i) {
    idx <- triangles[i, ]
    data.frame(
      x = c(mesh_coords[idx, 1], mesh_coords[idx[1], 1]),
      y = c(mesh_coords[idx, 2], mesh_coords[idx[1], 2]),
      triangle = i
    )
  })
  triangle_df <- do.call(rbind, triangle_list)
  
  p <- ggplot() +
    geom_path(data = triangle_df, 
              aes(x = x, y = y, group = triangle),
              color = "gray70", linewidth = 0.3, alpha = 0.5) +
    geom_point(data = mesh_df, aes(x = x, y = y),
               size = 0.5, color = "blue", alpha = 0.5) +
    geom_sf(data = data_list$lake_boundary, 
            fill = NA, color = "red", linewidth = 1) +
    geom_sf(data = data_list$perm_water_polygon, 
            fill = "lightblue", color = "darkblue", alpha = 0.3, linewidth = 0.5) +
    coord_sf() +
    labs(title = sprintf("SPDE Mesh for %s", data_list$lake_name),
         subtitle = sprintf("%d vertices, %d triangles", 
                            mesh$n, nrow(triangles))) +
    theme_minimal() +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5)
    )
  
  print(p)
  
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    cat(sprintf("   Mesh plot saved to: %s\n", save_path))
  }
  
  return(p)
}
