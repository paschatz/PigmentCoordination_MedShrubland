# Download climatic data from the region:
library(tidyverse)
library(climaemet)
library(fitdistrplus)

# Find nearest station:
aemet_stations() %>%
  # find stations close to 40.000000, -3.600000
  filter(
    abs(latitud - 40.000000) < 0.5,
    abs(longitud - -3.600000) < 0.5
  )

# Extract Aranjuez code:
aranj <- "3100B"

## Get api key from AEMET
#browseURL("https://opendata.aemet.es/centrodedescargas/obtencionAPIKey")

## Use this function to register your API Key temporarly or permanently
# REceived the API in my email. 
aemet_api_key("eyJhbGciOiJIUzI1NiJ9.eyJzdWIiOiJwYXNjaGFsaXMuY2hhdHpvcG91bG9zQGdtLnVjYS5lcyIsImp0aSI6Ijg4YmVlZDM3LTdiZDUtNDYwYy1iNGEwLTg4MmQzNmE3ZjM5YyIsImlzcyI6IkFFTUVUIiwiaWF0IjoxNjg2MjE0ODk0LCJ1c2VySWQiOiI4OGJlZWQzNy03YmQ1LTQ2MGMtYjRhMC04ODJkMzZhN2YzOWMiLCJyb2xlIjoiIn0.Od2RGc1y6XFkp7dtQeGuI61G5svO7m53tfXvUrtHq-k", 
install = TRUE)

# Test data to see if are available data:
test <- aemet_daily_clim(aranj, start = "2011-01-01", end = "2011-12-31")

# Download climatic data seperately:
data_daily_0 <-
  aemet_daily_clim(aranj, start = "2005-01-01", end = "2009-12-31")

data_daily_1 <-
  aemet_daily_clim(aranj, start = "2010-01-01", end = "2016-12-31")

data_daily_2 <-
  aemet_daily_clim(aranj, start = "2017-01-01", end = "2024-10-15")

# Merge datasets
data <- bind_rows(data_daily_0, data_daily_1, data_daily_2)

#write.csv(data, "scripts//climatic_data_aranjuez.csv", row.names = FALSE)

data <- read.csv("scripts//climatic_data_aranjuez.csv")

# Tity data:
data.tidy <- data %>%
  mutate(
    prec = as.numeric(prec),
    tmed = as.numeric(tmed),
    tmin = as.numeric(tmin),
    tmax = as.numeric(tmax)) %>%
  # split date into year, month, day
  filter(!is.na(prec))

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

# Plot sum of precipitation with a continuous color scale and dashed line for 2011
plot_1 <- ggplot(avg_prec, aes(x = month, y = month_prec, group = year)) +
  # Solid lines for all years except 2011
  geom_line(aes(color = year), size = 0.8, data = subset(avg_prec, year != 2011)) +
  # Dashed line specifically for the year 2011
  geom_line(aes(color = year), size = 0.8, linetype = "dashed", data = subset(avg_prec, year == 2011)) +
  geom_point(aes(color = year), size = 1) +
  scale_x_continuous(breaks = 1:12, labels = month.abb) +
  scale_color_viridis_c() +
  labs(
    title = "",
    x = "Month",
    y = "Precipitation (mm)",
    color = "Year") +
  theme_minimal() +
  # Background blank
  theme(panel.grid.major = element_blank(), 
        panel.grid.minor = element_blank())

# Export precipitation graph:
#ggsave("exports///precipitation.png", plot_1, width = 20, height = 12, units = "cm", dpi = 300)

# Calculate the Standardized Precipitation Index (SPI) for the region:
# Step 1: Prepare the data by summarizing total monthly precipitation
precip_data <- data.tidy %>%
  group_by(year, month) %>%
  summarise(total_precip = sum(prec, na.rm = TRUE), .groups = 'drop') %>%
  mutate(year_month = as.Date(paste(year, month, "01", sep = "-")))

# Step 2: Fit the gamma distribution to the long-term precipitation data
fit <- fitdist(precip_data$total_precip, "gamma", method = "mme")

# Step 3: Calculate the CDF for each total precipitation value using the total_precip column
precip_data <- precip_data %>%
  mutate(cdf = pgamma(total_precip, shape = fit$estimate["shape"], rate = fit$estimate["rate"]))

# Step 4: Transform the CDF to a normal distribution (using the quantile function)
precip_data <- precip_data %>%
  mutate(spi = qnorm(cdf))

# Print the results
print(precip_data)

spi <- ggplot(precip_data, aes(x = year_month, y = spi)) +
  geom_line(color = "blue", size = 1) +   # Line for SPI values
  geom_hline(yintercept = 0, linetype = "dashed", color = "red") +  # Add a horizontal line at y = 0
  geom_ribbon(aes(ymin = -2, ymax = 2), fill = "lightblue", alpha = 0.3) +  # Add a shaded area for SPI range
  labs(title = "Standardized Precipitation Index (SPI)",
       x = "Date",
       y = "SPI Value") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
        axis.title.x = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text = element_text(size = 10),
        panel.grid.minor = element_blank()) +
  # increse X by one year
  scale_x_date(date_breaks = "1 year", date_labels = "%Y")

spi 
# Export SPI graph:
ggsave("exports///spi.png", spi, width = 30, height = 20, units = "cm", dpi = 300)
