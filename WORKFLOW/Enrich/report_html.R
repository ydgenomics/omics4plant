# 创建交互式HTML报告
library(htmltools)
library(DT)

create_html_image_report <- function(image_dir = ".", 
                                     output_file = "image_report.html",
                                     title = "图片分析报告") {
  library(htmltools)
  library(DT)
  library(jsonlite)

  image_files <- list.files(image_dir, 
                            pattern = "\\.(png|jpg|jpeg|gif|bmp|svg|pdf)$",
                            full.names = TRUE)
  if (length(image_files) == 0) stop("未找到图片文件")

  # Copy images to a subfolder relative to the output HTML file
  img_dir <- file.path(dirname(output_file), "report_images")
  if (!dir.exists(img_dir)) dir.create(img_dir, recursive = TRUE)
  img_relative_paths <- file.path("report_images", basename(image_files))
  mapply(file.copy, image_files, file.path(img_dir, basename(image_files)), overwrite = TRUE)

  image_info <- data.frame(
    序号 = seq_along(image_files),
    文件名 = basename(image_files),
    文件类型 = tools::file_ext(image_files),
    大小_MB = round(file.info(image_files)$size / 1024^2, 3),
    修改时间 = format(file.info(image_files)$mtime, "%Y-%m-%d %H:%M"),
    stringsAsFactors = FALSE
  )

  output_table <- datatable(
    image_info,
    options = list(pageLength = 10),
    class = 'display',
    escape = FALSE
  )

  html_content <- tags$html(
    tags$head(
      tags$title(title),
      tags$style(HTML("
        body { font-family: Arial; margin: 20px; }
        .gallery { display: flex; flex-wrap: wrap; gap: 15px; }
        .gallery-item { border: 1px solid #ddd; border-radius: 5px; padding: 10px; text-align: center; width: 220px;}
        .gallery-item img { max-width: 100%; height: 120px; object-fit: contain; cursor:pointer;}
        .modal { display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background: rgba(0,0,0,0.9);}
        .modal-content { margin: auto; display: block; max-width: 90%; max-height: 90%; }
        .close { position: absolute; top: 15px; right: 35px; color: #fff; font-size: 40px; font-weight: bold; cursor: pointer; }
      "))
    ),
    tags$body(
      h1(title),
      p(paste("生成时间:", Sys.Date())),
      p(paste("图片数量:", length(image_files))),
      div(class = "gallery",
          lapply(seq_along(img_relative_paths), function(i) {
            img <- img_relative_paths[i]
            if (grepl("\\.(png|jpg|jpeg|gif)$", img, ignore.case = TRUE)) {
              div(class = "gallery-item",
                  img(src = img, onclick = sprintf("showModal('%s')", img)),
                  p(basename(img)),
                  p(paste("大小:", round(file.info(file.path(img_dir, basename(img)))$size/1024^2, 2), "MB"))
              )
            }
          })
      ),
      h3("图片信息表"),
      HTML(as.character(output_table)),
      div(id = "imageModal", class = "modal",
          span(class = "close", onclick = "closeModal()", "×"),
          img(class = "modal-content", id = "modalImage")
      ),
      tags$script(HTML(sprintf("
        function showModal(imgSrc) {
          document.getElementById('imageModal').style.display = 'block';
          document.getElementById('modalImage').src = imgSrc;
        }
        function closeModal() {
          document.getElementById('imageModal').style.display = 'none';
        }
        window.onclick = function(event) {
          var modal = document.getElementById('imageModal');
          if (event.target == modal) closeModal();
        }
      ")))
    )
  )

  save_html <- function(html_content, output_file) {
    full_html <- paste(
      '<!DOCTYPE html>',
      '<html lang="zh-CN">',
      '<head>',
      '<meta charset="UTF-8">',
      '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
      '</head>',
      as.character(html_content),
      '</html>',
      sep = "\n"
    )
    writeLines(full_html, output_file)
    message("报告已保存到: ", output_file)
  }

  save_html(html_content, output_file)
#   if (interactive()) browseURL(output_file, browser = getOption("browser", ""))
  return(output_file)
}
# 用于批量SVG图片的HTML报告生成函数
create_svg_html_report <- function(svg_dir = ".", 
                                  output_file = "svg_report.html",
                                  title = "SVG图片报告") {
  svg_files <- list.files(svg_dir, pattern = "\\.svg$", full.names = TRUE)
  if (length(svg_files) == 0) stop("未找到SVG图片文件")
  
  img_dir <- file.path(dirname(output_file), "report_svgs")
  if (!dir.exists(img_dir)) dir.create(img_dir, recursive = TRUE)
  svg_relative_paths <- file.path("report_svgs", basename(svg_files))
  mapply(file.copy, svg_files, file.path(img_dir, basename(svg_files)), overwrite = TRUE)
  
  svg_info <- data.frame(
    序号 = seq_along(svg_files),
    文件名 = basename(svg_files),
    大小_KB = round(file.info(svg_files)$size / 1024, 2),
    修改时间 = format(file.info(svg_files)$mtime, "%Y-%m-%d %H:%M"),
    stringsAsFactors = FALSE
  )
  
  output_table <- datatable(
    svg_info,
    options = list(pageLength = 10),
    class = 'display',
    escape = FALSE
  )
  
  html_content <- tags$html(
    tags$head(
      tags$title(title),
      tags$style(HTML("
        body { font-family: Arial; margin: 20px; }
        .gallery { display: flex; flex-wrap: wrap; gap: 15px; }
        .gallery-item { border: 1px solid #ddd; border-radius: 5px; padding: 10px; text-align: center; width: 320px;}
        .gallery-item object { width: 300px; height: 200px; cursor:pointer; border: none;}
        .modal { display: none; position: fixed; z-index: 1000; left: 0; top: 0; width: 100%; height: 100%; background: rgba(0,0,0,0.9);}
        .modal-content { margin: auto; display: block; max-width: 90%; max-height: 90%; }
        .close { position: absolute; top: 15px; right: 35px; color: #fff; font-size: 40px; font-weight: bold; cursor: pointer; }
      "))
    ),
    tags$body(
      h1(title),
      p(paste("生成时间:", Sys.Date())),
      p(paste("SVG图片数量:", length(svg_files))),
      div(class = "gallery",
          lapply(seq_along(svg_relative_paths), function(i) {
            svg <- svg_relative_paths[i]
            div(class = "gallery-item",
                tags$object(data = svg, type = "image/svg+xml", 
                            onclick = sprintf("showModal('%s')", svg)),
                p(basename(svg)),
                p(paste("大小:", round(file.info(file.path(img_dir, basename(svg)))$size/1024, 2), "KB"))
            )
          })
      ),
      h3("SVG图片信息表"),
      HTML(as.character(output_table)),
      div(id = "imageModal", class = "modal",
          span(class = "close", onclick = "closeModal()", "×"),
          tags$object(class = "modal-content", id = "modalImage", type = "image/svg+xml")
      ),
      tags$script(HTML("
        function showModal(svgSrc) {
          document.getElementById('imageModal').style.display = 'block';
          document.getElementById('modalImage').data = svgSrc;
        }
        function closeModal() {
          document.getElementById('imageModal').style.display = 'none';
        }
        window.onclick = function(event) {
          var modal = document.getElementById('imageModal');
          if (event.target == modal) closeModal();
        }
      "))
    )
  )
  
  save_html <- function(html_content, output_file) {
    full_html <- paste(
      '<!DOCTYPE html>',
      '<html lang="zh-CN">',
      '<head>',
      '<meta charset="UTF-8">',
      '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
      '</head>',
      as.character(html_content),
      '</html>',
      sep = "\n"
    )
    writeLines(full_html, output_file)
    message("SVG报告已保存到: ", output_file)
  }
  
  save_html(html_content, output_file)
  return(output_file)
}

# 使用示例
# 请确保 image_dir 指向包含图片文件的文件夹，而不是 CSV 文件或其他非图片文件夹
# 例如:
create_html_image_report(
  image_dir = "/data/work/yita/images",
  output_file = "my_image_report.html",
  title = "我的图片分析报告"
)