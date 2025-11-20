# Download climatic data from the region:
library(tidyverse)
library(ggrepel)

# Set plotting theme:
theme.plot <- theme_minimal() +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black", size = 6),
    axis.title = element_text(size = 6),
    panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.ticks = element_line(color = "black", size = 0.2),
    axis.ticks.length = unit(0.6, "mm")
  )

data <- read.csv("data/climatic.data.csv")

# Tity data:
data.tidy <- data %>%
  mutate(
    fecha = as.Date(fecha),
    prec = as.numeric(prec),
    tmed = as.numeric(tmed),
    tmin = as.numeric(tmin),
    tmax = as.numeric(tmax)) %>%
  separate(fecha, into = c("year", "month", "day"), sep = "-") %>%
  # split date into year, month, day
  filter(!is.na(prec))

# Calculate mean annual precipitation for each year
data.tidy %>%
  group_by(year) %>%
  summarise(annual_prec = sum(prec, na.rm = TRUE))

# Calculate mean September precipitation per years
sep.mean <- data.tidy %>%
  group_by(year, month) %>%
  summarise(annual_prec = sum(prec, na.rm = TRUE)) %>%
  filter(month == "09")

# Calculate mean September precipitation across all years
sep.mean %>%
  ungroup() %>%
  summarise(mean.sep = mean(annual_prec, na.rm = TRUE))

# Calculate total precipitation per month in 2011:
data.tidy %>%
  filter(year == "2011") %>%
  group_by(month) %>%
  summarise(monthly.prec = sum(prec, na.rm = TRUE))

# Calculate mean annual precipitation across all years:
data.tidy %>%
  group_by(year) %>%
  summarise(annual.prec = sum(prec, na.rm = TRUE)) %>%
  summarise(mean.prec = mean(annual.prec, na.rm = TRUE))

# Data only for 2011:
data_2011 <- data.tidy %>%
  filter(year == "2011")

# Calculate monthly precipitation:
avg_prec <- data.tidy %>%
  group_by(year, month) %>%
  summarise(month_prec = sum(prec, na.rm = TRUE)) %>%
  # month as numeric
  mutate(month = as.numeric(month),
         year = as.numeric(year)) 

# Prepare data for plotting:
plot_data <- data_2011 %>%
  mutate(date = make_date(year, month, day)) %>%
  arrange(date)

# Define key dates:
sampling_drought <- as.Date("2011-10-21")
sampling_recovery <- as.Date("2011-12-27")

drought_start <- as.Date("2011-06-07")
drought_end   <- as.Date("2011-10-22")

recovery_start <- as.Date("2011-10-22")
recovery_end   <- as.Date("2011-12-27")

# Create the Plot
ggplot(plot_data, aes(x = date, y = prec)) +
  annotate("rect", xmin = drought_start, xmax = drought_end, 
           ymin = 0, ymax = Inf, 
           fill = "#d24644", alpha = 0.3) + # Light red background
  annotate("rect", xmin = recovery_start, xmax = recovery_end, 
           ymin = 0, ymax = Inf, 
           fill = "#3fbc73", alpha = 0.3) + # Light red background
  
  geom_col(fill = "#3b528b", width = 1, linewidth = 0, color = NA) +
  
  # Customize the scale
  scale_x_date(
    date_breaks = "1 month", 
    date_labels = "%b",
    limits = as.Date(c("2011-01-01", "2011-12-31")),
    expand = c(0, 0)) +
  
  # Mark Sampling Dates (Vertical Dashed Lines)
  geom_vline(xintercept = sampling_drought, linetype = "dashed", color = "#A6611A", size = 0.5) +
  geom_vline(xintercept = sampling_recovery, linetype = "dashed", color = "#018571", size = 0.5) +
  
  annotate("text", x = sampling_drought, y = 20, label = "Drought\nSampling", 
           vjust = 1, hjust = 1.1, size = 2, color = "#A6611A") +
  annotate("text", x = sampling_recovery, y = 20, label = "Recovery\nSampling", 
           vjust = 1, hjust = 1.1, size = 2, color = "#018571") +
  
  geom_text_repel(
    data = subset(plot_data, prec != "14.8" & 
                    prec != "9.6"  & prec != "19" & prec != "0" &
                    date >= drought_start & date <= drought_end), 
    aes(label = paste0(prec, " mm")),
    nudge_y = 2,
    nudge_x = 1,
    box.padding = 0.5,
    point.padding = 0.2,
    segment.color = "gray50",
    segment.size = 0.1,
    size = 1.5, 
    color = "gray30") +
  
  scale_y_continuous(expand = c(0, 0), limits = c(0, 22)) +

  labs(
    x = "Month",
    y = "Daily precipitation (mm)", # Using expression for proper formatting
    title = "") +
  
  # Clean theme
  theme.plot

ggsave("exports/figures//Figure1.tiff", width = 16, height = 10, units = "cm", dpi = 1200)
ggsave("exports/figures//Figure1.pdf", width = 16, height = 10, units = "cm", dpi = 1200)
