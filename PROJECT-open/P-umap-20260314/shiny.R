# app.R
library(shiny)
library(shinythemes)
library(plotly)
library(DT)
library(tidyverse)
library(Seurat)
library(viridis)

# 加载预处理的数据
load("./DATA/shiny_data.RData")

# UI界面
ui <- navbarPage(
  "单细胞交互式质控平台",
  theme = shinytheme("flatly"),
  
  # 主标签页：UMAP可视化
  tabPanel("UMAP质控",
    fluidPage(
      fluidRow(
        column(3,
          wellPanel(
            h4("质控参数设置"),
            
            # 基因选择
            selectizeInput("gene_select", 
                          "选择基因:",
                          choices = rownames(counts_matrix)[1:100],  # 可调整为所有基因
                          selected = rownames(counts_matrix)[1],
                          multiple = FALSE),
            
            # 分组变量
            selectInput("group_var", 
                       "分组变量:",
                       choices = colnames(metadata_with_umap),
                       selected = "orig.ident"),
            
            # 颜色方案
            selectInput("color_scheme", 
                       "颜色方案:",
                       choices = c("Viridis" = "viridis",
                                 "Plasma" = "plasma",
                                 "Magma" = "magma",
                                 "Inferno" = "inferno")),
            
            # 点大小控制
            sliderInput("point_size", 
                       "点大小:",
                       min = 0.1, max = 3, value = 0.5, step = 0.1),
            
            # 透明度控制
            sliderInput("point_alpha", 
                       "透明度:",
                       min = 0.1, max = 1, value = 0.7, step = 0.1),
            
            hr(),
            h4("小提琴图质控"),
            
            # 表达量阈值滑块
            uiOutput("violin_slider")
          )
        ),
        
        column(9,
          fluidRow(
            column(6, plotlyOutput("umap_plot", height = "500px")),
            column(6, plotlyOutput("violin_plot", height = "500px"))
          ),
          fluidRow(
            column(12,
              wellPanel(
                h4("统计信息"),
                verbatimTextOutput("stats_info")
              )
            )
          )
        )
      )
    )
  ),
  
  # 数据表格标签页
  tabPanel("数据表格",
    fluidPage(
      fluidRow(
        column(12,
          DTOutput("metadata_table")
        )
      )
    )
  ),
  
  # 高级分析标签页
  tabPanel("高级分析",
    fluidPage(
      fluidRow(
        column(6,
          h4("特征图对比"),
          plotlyOutput("feature_plot1", height = "400px"),
          plotlyOutput("feature_plot2", height = "400px")
        ),
        column(6,
          h4("聚类统计"),
          plotOutput("cluster_stats", height = "800px")
        )
      )
    )
  )
)

# Server逻辑
server <- function(input, output, session) {
  
  # 动态创建小提琴图滑块
  output$violin_slider <- renderUI({
    req(input$gene_select)
    
    # 获取选定基因的表达量范围
    gene_counts <- counts_matrix[input$gene_select, ]
    min_val <- min(gene_counts)
    max_val <- max(gene_counts)
    
    tagList(
      sliderInput("expr_threshold",
                 label = paste(input$gene_select, "表达量阈值:"),
                 min = min_val,
                 max = max_val,
                 value = c(min_val, max_val),
                 step = (max_val - min_val)/100),
      
      checkboxInput("invert_selection",
                   "反向选择",
                   value = FALSE)
    )
  })
  
  # 根据阈值筛选细胞
  filtered_cells <- reactive({
    req(input$gene_select, input$expr_threshold)
    
    gene_expr <- counts_matrix[input$gene_select, ]
    
    # 应用阈值筛选
    if(input$invert_selection) {
      keep_cells <- gene_expr < input$expr_threshold[1] | 
                   gene_expr > input$expr_threshold[2]
    } else {
      keep_cells <- gene_expr >= input$expr_threshold[1] & 
                   gene_expr <= input$expr_threshold[2]
    }
    
    return(keep_cells)
  })
  
  # 生成UMAP图
  output$umap_plot <- renderPlotly({
    req(metadata_with_umap, filtered_cells())
    
    # 准备数据
    plot_data <- metadata_with_umap
    plot_data$selected <- filtered_cells()
    plot_data$gene_expr <- counts_matrix[input$gene_select, ]
    
    # 根据分组变量设置颜色
    if(input$group_var %in% colnames(metadata_with_umap)) {
      color_col <- metadata_with_umap[[input$group_var]]
      
      p <- plot_data %>%
        plot_ly(
          x = ~umap1,
          y = ~umap2,
          color = ~eval(parse(text = input$group_var)),
          colors = viridis_pal(option = input$color_scheme)(10),
          type = 'scatter',
          mode = 'markers',
          marker = list(
            size = input$point_size,
            opacity = input$point_alpha
          ),
          text = ~paste(
            "Cell: ", rownames(metadata_with_umap),
            "<br>Group: ", eval(parse(text = input$group_var)),
            "<br>", input$gene_select, ": ", round(gene_expr, 2),
            "<br>Selected: ", selected
          ),
          hoverinfo = 'text'
        )
    } else {
      # 连续变量（表达量）着色
      p <- plot_data %>%
        plot_ly(
          x = ~umap1,
          y = ~umap2,
          color = ~gene_expr,
          colors = viridis_pal(option = input$color_scheme)(100),
          type = 'scatter',
          mode = 'markers',
          marker = list(
            size = input$point_size,
            opacity = input$point_alpha,
            showscale = TRUE,
            colorbar = list(title = input$gene_select)
          ),
          text = ~paste(
            "Cell: ", rownames(metadata_with_umap),
            "<br>", input$gene_select, ": ", round(gene_expr, 2),
            "<br>Selected: ", selected
          ),
          hoverinfo = 'text'
        )
    }
    
    # 高亮筛选的细胞
    p <- p %>%
      add_trace(
        data = plot_data[plot_data$selected, ],
        x = ~umap1,
        y = ~umap2,
        type = 'scatter',
        mode = 'markers',
        marker = list(
          size = input$point_size + 0.5,
          color = 'black',
          opacity = 0.3,
          line = list(color = 'red', width = 2)
        ),
        showlegend = FALSE,
        inherit = FALSE
      )
    
    p %>%
      layout(
        title = paste("UMAP降维图 - 高亮筛选细胞"),
        xaxis = list(title = "UMAP1"),
        yaxis = list(title = "UMAP2")
      )
  })
  
  # 生成小提琴图
  output$violin_plot <- renderPlotly({
    req(input$gene_select)
    
    plot_data <- data.frame(
      expression = counts_matrix[input$gene_select, ],
      group = metadata_with_umap[[input$group_var]],
      selected = filtered_cells()
    )
    
    # 创建小提琴图
    p <- plot_data %>%
      plot_ly(
        x = ~group,
        y = ~expression,
        type = 'violin',
        box = list(visible = TRUE),
        meanline = list(visible = TRUE),
        color = ~group,
        colors = viridis_pal(option = input$color_scheme)(length(unique(plot_data$group))),
        split = ~group,
        points = FALSE,
        hoverinfo = 'y+name'
      ) %>%
      layout(
        title = paste(input$gene_select, "表达分布"),
        xaxis = list(title = input$group_var),
        yaxis = list(title = "表达量"),
        shapes = list(
          list(
            type = "line",
            x0 = -1,
            x1 = length(unique(plot_data$group)),
            y0 = input$expr_threshold[1],
            y1 = input$expr_threshold[1],
            line = list(color = "red", dash = "dash")
          ),
          list(
            type = "line",
            x0 = -1,
            x1 = length(unique(plot_data$group)),
            y0 = input$expr_threshold[2],
            y1 = input$expr_threshold[2],
            line = list(color = "red", dash = "dash")
          )
        )
      )
    
    return(p)
  })
  
  # 统计信息
  output$stats_info <- renderPrint({
    req(filtered_cells())
    
    total_cells <- nrow(metadata_with_umap)
    selected_cells <- sum(filtered_cells())
    percentage <- round(selected_cells / total_cells * 100, 2)
    
    cat("=== 质控统计 ===\n")
    cat("总细胞数:", total_cells, "\n")
    cat("筛选后细胞数:", selected_cells, "\n")
    cat("保留比例:", percentage, "%\n")
    cat("\n分组统计:\n")
    
    # 按分组统计
    group_counts <- table(metadata_with_umap[[input$group_var]])
    selected_group_counts <- table(metadata_with_umap[[input$group_var]][filtered_cells()])
    
    for(group in names(group_counts)) {
      total <- group_counts[group]
      selected <- ifelse(group %in% names(selected_group_counts), 
                        selected_group_counts[group], 0)
      perc <- round(selected / total * 100, 2)
      cat(sprintf("  %s: %d/%d (%.2f%%)\n", group, selected, total, perc))
    }
  })
  
  # 数据表格
  output$metadata_table <- renderDT({
    filtered_data <- metadata_with_umap
    filtered_data$selected <- filtered_cells()
    filtered_data$gene_expr <- counts_matrix[input$gene_select, ]
    
    datatable(filtered_data,
              options = list(
                pageLength = 20,
                scrollX = TRUE,
                dom = 'Bfrtip'
              )) %>%
      formatStyle(
        'selected',
        backgroundColor = styleEqual(c(TRUE, FALSE), c('#90EE90', '#FFB6C1'))
      )
  })
  
  # 高级分析：特征图对比
  output$feature_plot1 <- renderPlotly({
    req(input$gene_select)
    
    plot_data <- metadata_with_umap
    plot_data$feature <- counts_matrix[input$gene_select, ]
    
    plot_ly(plot_data,
            x = ~umap1,
            y = ~umap2,
            z = ~feature,
            type = 'mesh3d',
            intensity = ~feature,
            colorscale = 'Viridis') %>%
      layout(title = paste(input$gene_select, "3D分布"))
  })
  
  output$feature_plot2 <- renderPlotly({
    req(input$gene_select)
    
    plot_data <- metadata_with_umap
    plot_data$feature <- counts_matrix[input$gene_select, ]
    
    plot_data %>%
      plot_ly(
        x = ~umap1,
        y = ~umap2,
        type = 'histogram2d',
        colorscale = 'Viridis'
      ) %>%
      layout(title = "细胞密度分布")
  })
  
  # 聚类统计
  output$cluster_stats <- renderPlot({
    req(input$group_var)
    
    p1 <- metadata_with_umap %>%
      count(!!sym(input$group_var)) %>%
      ggplot(aes(x = reorder(!!sym(input$group_var), n), y = n, fill = !!sym(input$group_var))) +
      geom_bar(stat = "identity") +
      theme_minimal() +
      labs(x = input$group_var, y = "细胞数", title = "细胞数分布") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    p2 <- metadata_with_umap %>%
      select(!!sym(input$group_var), starts_with("percent.")) %>%
      pivot_longer(cols = starts_with("percent.")) %>%
      ggplot(aes(x = !!sym(input$group_var), y = value, fill = !!sym(input$group_var))) +
      geom_boxplot() +
      facet_wrap(~name, scales = "free_y") +
      theme_minimal() +
      labs(x = input$group_var, y = "百分比", title = "质控指标分布") +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    cowplot::plot_grid(p1, p2, ncol = 1, rel_heights = c(1, 2))
  })
}

# 运行应用
shinyApp(ui = ui, server = server)