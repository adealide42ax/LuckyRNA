library(ggplot2)
library(ggpubr)

files <- list.files(pattern = "*.csv", recursive = TRUE)

lollipop <- function(file,filename,title) {
  data = read.csv(file,header = TRUE)
  if(nrow(data) <= 1) return()
  data = sapply(data, table)
  #得到计算Description频数的结果
  data_Description = data$Description
  #将结果转换为data.frame
  data_Description = as.data.frame(data_Description)
  #计算行数
  rows <- nrow(data_Description)
  #按频数排序
  data_Description = data_Description[order(data_Description$Freq,data_Description$Var1,decreasing = T),]
  
  # 创建一个名为'top'的新列
  data_Description$Top <- 0
  
  # 循环遍历数据帧的每一行
  for (i in 1:nrow(data_Description)) {
    # 计算当前行的top值
    top_value <- 10 * ceiling(i / 10)
    
    # 将当前行的top值赋给新列
    data_Description[i, "Top"] <- top_value
  }
  
  # 如果剩余的行数不足十行，则将实际的行数赋给top值
  if (nrow(data_Description) < 10) {
    data_Description$Top[nrow(data_Description)] <- nrow(data_Description)
  } else
    if (nrow(data_Description) %% 10 == 0) {
      # 如果数据帧的行数为整十的倍数，则将最后十行的 "Top" 值设为实际的行数
      data_Description$Top[(nrow(data_Description) - 9):nrow(data_Description)] <- nrow(data_Description)
    } else {
      # 如果剩余的行数超过十行，则将最后一段小于十行的部分的top值设为实际的行数
      data_Description$Top[(nrow(data_Description) - (nrow(data_Description) %% 10) + 1):nrow(data_Description)] <- nrow(data_Description)
    }
  
  data_Description <- data_Description[nrow(data_Description):1, ]
  data_Description$Description = as.character(data_Description$Var1)
  data_Description$Freq = as.numeric(data_Description$Freq)
  data_Description$Top = as.factor(data_Description$Top)
  
  freq_values <- data_Description$Freq
  types <- length(unique(freq_values))
  types <- types + 8
  colors <- colorRampPalette(c("#66ccff","#39c5bb","#FFE211","#FFA500", "#EE0000"))(types)
  lollipop_pic <- ggdotchart(data_Description, x="Description", y="Freq",
                             col = colors[data_Description$Freq],
                             sorting = "descending", add = "segments", rotate = TRUE,
                             dot.size = 7,title=title,
                             label = round(data_Description$Freq), font.label = list(color="black", size=9, vjust=0.5), ggtheme = theme_pubr())
  return(lollipop_pic)
}

for (file in files) {
  file <- file
  filenames = strsplit(file, "[/_.]")
  # 检查列表中是否存在 "GO" 元素
  if (any("GO" %in% filenames[[1]])) 
  {
    filename <- paste0(filenames[[1]][length(filenames[[1]]) - 6], "_", filenames[[1]][length(filenames[[1]]) - 4], "_", filenames[[1]][length(filenames[[1]]) - 2], "_", filenames[[1]][length(filenames[[1]]) - 1])
    title = paste0(filenames[[1]][length(filenames[[1]]) - 6],"_",filenames[[1]][length(filenames[[1]]) - 4],"_",filenames[[1]][length(filenames[[1]]) - 3],"_","pathway")
    filename = paste(dirname(file),"/",filename,".pdf",sep = "")
    if (grepl("BP", filenames)& grepl("diff", filenames)) {
      bp_diff <- lollipop(file,filename,title)
    } else if (grepl("CC", filenames)& grepl("diff", filenames)) {
      cc_diff <- lollipop(file,filename,title)
    } else if (grepl("MF", filenames)& grepl("diff", filenames)){
      mf_diff <- lollipop(file,filename,title)
      lollipop_pic_GO_diff <- ggarrange(bp_diff,cc_diff,mf_diff,ncol = 1,nrow =3,widths = c(2,2,2),heights = c(1))
      ggsave(lollipop_pic_GO_diff, filename =filename,  width = 15, height = 45)
    }
      else if (grepl("BP", filenames)& grepl("Down", filenames)) {
      bp_down <- lollipop(file,filename,title)
    } else if (grepl("CC", filenames)& grepl("Down", filenames)) {
      cc_down <- lollipop(file,filename,title)
    } else if (grepl("MF", filenames)& grepl("Down", filenames)){
      mf_down <- lollipop(file,filename,title)
      lollipop_pic_GO_down <- ggarrange(bp_down,cc_down,mf_down,ncol = 1,nrow =3,widths = c(2,2,2),heights = c(1))
      ggsave(lollipop_pic_GO_down, filename =filename,  width = 15, height = 45)
    }
      else if (grepl("BP", filenames)& grepl("Up", filenames)) {
      bp_up <- lollipop(file,filename,title)
    } else if (grepl("CC", filenames)& grepl("Up", filenames)) {
      cc_up <- lollipop(file,filename,title)
    } else if (grepl("MF", filenames)& grepl("Up", filenames)){
      mf_up <- lollipop(file,filename,title)
      lollipop_pic_GO_up <- ggarrange(bp_up,cc_up,mf_up,ncol = 1,nrow =3,widths = c(2,2,2),heights = c(1))
      ggsave(lollipop_pic_GO_up, filename =filename,  width = 15, height = 45)
    }
  }
  
  else {
    filename <- paste0(filenames[[1]][length(filenames[[1]]) - 5], "_", filenames[[1]][length(filenames[[1]]) - 3], "_", filenames[[1]][length(filenames[[1]]) - 2], "_", filenames[[1]][length(filenames[[1]]) - 1])
    title = paste(filename,"pathway",sep = "_")
    filename = paste(dirname(file),"/",filename,".pdf",sep = "")
    lollipop_pic_KEGG <-lollipop(file,filename,title)
    ggsave(lollipop_pic_KEGG, filename =filename,  width = 15, height = 10)
  }
}
