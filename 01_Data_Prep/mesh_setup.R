################################################################################
# SPDE Mesh Construction Module
# Authors: Chaoan Li, Yinuo Zhu | STAT 647, Texas A&M University
################################################################################

#' 构建 INLA SPDE mesh
#' Build INLA SPDE mesh over lake region
#' 
#' @param data_list load_and_prep_data() 返回的列表
#' @param max_edge 三角形最大边长（内部，外部）向量，单位与坐标系一致
#' @param cutoff 最小节点间距
#' @param offset 边界外扩展距离（内部，外部）向量
#' @return mesh 对象
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
  
  # ---- 1. 提取边界坐标 ----
  cat("\n1. Extracting boundary coordinates...\n")
  
  # 获取边界的外环坐标
  boundary_coords <- st_coordinates(lake_boundary)[, 1:2]
  
  # 确保是唯一的坐标（去除首尾重复点）
  if (all(boundary_coords[1, ] == boundary_coords[nrow(boundary_coords), ])) {
    boundary_coords <- boundary_coords[-nrow(boundary_coords), ]
  }
  
  cat(sprintf("   Boundary points: %d\n", nrow(boundary_coords)))
  
  # ---- 2. 构建 mesh ----
  cat("\n2. Building mesh...\n")
  cat(sprintf("   max.edge: c(%.1f, %.1f)\n", max_edge[1], max_edge[2]))
  cat(sprintf("   cutoff: %.1f\n", cutoff))
  cat(sprintf("   offset: c(%.1f, %.1f)\n", offset[1], offset[2]))
  
  # 使用 inla.mesh.2d 构建 mesh
  mesh <- inla.mesh.2d(
    loc.domain = boundary_coords,
    max.edge = max_edge,
    cutoff = cutoff,
    offset = offset
  )
  
  cat(sprintf("\n✓ Mesh built successfully!\n"))
  cat(sprintf("   Number of vertices: %d\n", mesh$n))
  cat(sprintf("   Number of triangles: %d\n", nrow(mesh$graph$tv)))
  
  return(mesh)
}


#' 可视化 mesh
#' Visualize the mesh
#' 
#' @param mesh INLA mesh 对象
#' @param data_list load_and_prep_data() 返回的列表
#' @param save_path 保存图片的路径（可选）
plot_mesh <- function(mesh, data_list, save_path = NULL) {
  
  library(ggplot2)
  library(sf)
  
  cat("Plotting mesh...\n")
  
  # 提取 mesh 三角形顶点
  mesh_coords <- mesh$loc[, 1:2]
  mesh_df <- data.frame(x = mesh_coords[, 1], y = mesh_coords[, 2])
  
  # 提取三角形
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
  
  # 绘图
  p <- ggplot() +
    # 绘制三角形边
    geom_path(data = triangle_df, 
              aes(x = x, y = y, group = triangle),
              color = "gray70", linewidth = 0.3, alpha = 0.5) +
    # 绘制 mesh 节点
    geom_point(data = mesh_df, aes(x = x, y = y),
               size = 0.5, color = "blue", alpha = 0.5) +
    # 叠加湖区边界
    geom_sf(data = data_list$lake_boundary, 
            fill = NA, color = "red", linewidth = 1) +
    # 叠加永久水域
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
  
  # 保存
  if (!is.null(save_path)) {
    ggsave(save_path, p, width = 10, height = 8, dpi = 300)
    cat(sprintf("   Mesh plot saved to: %s\n", save_path))
  }
  
  return(p)
}

