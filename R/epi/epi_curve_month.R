## Code used to create figure 1a 

library(tidyverse)
library(zoo)
library(lubridate)
library(ggthemes)

## Read in deaths from 2001-2023
## Also includes Vaccincation coverage and confirmed/probable cases
dat <- read_csv("data/clean/deaths_cases_2001_2023.csv")

## Transform from wide to long data for plotting
dat_cases = dat %>%
  select(!VC) %>% 
  pivot_longer(!year, names_to = "case_status", values_to = "n")

## Change vaccination coverage to percentage
dat_vac = dat %>% 
  select(year, VC) %>% 
  na.omit() %>% 
  mutate(VC = VC*100) 

## change case status to a factor for plotting
dat_cases = dat_cases %>% 
mutate(case_status_fct = 
         factor(case_status, 
                levels = c("deaths","confirmed", "probable")))

levels(dat_cases$case_status_fct) <- 
  c("Human Death", 'Confirmed (Animal)', "Probable (Animal)")
  
## Fig 1.A
yearly_plot = 
  dat_cases %>% 
  ggplot() +
  geom_rect(aes(xmin = 2020+0.25, xmax = 2022+0.25, ymin = 0, ymax = 40),
            fill = "#E8EBED", linewidth = 0, linetype = "solid") +
  scale_fill_manual(values=c("#001219","#c1121f","#DE7002")) +
  geom_col(aes(x = year, y = n, fill = case_status_fct), alpha = 0.9) +
  geom_segment(data = dat_vac, aes(x=year-0.4, xend = year+0.4, y=VC, yend = VC), lwd = 1.5, col = "#3a86ff") +
  scale_y_continuous(
    # Features of the first axis
    name = "Number of cases",
    # Add a second axis and specify its features
    sec.axis = sec_axis(~.*1, name="Vaccination coverage (%)")
  ) + 
  geom_vline(xintercept=2019+0.25, linetype=2, linewidth = 0.8, color = "#3a86ff") +
  scale_x_continuous(breaks=seq(2001,2023,1), name = "") +
  theme_clean(base_size = 14) +
  theme(legend.position = "top",
        legend.title = element_blank(),
        legend.background = element_rect(color = NA),
        axis.title.y.left = element_text(colour="#001219"),
        axis.text.y.left = element_text(colour="#001219"),
        axis.title.y.right = element_text(colour="#3a86ff"),
        axis.text.y.right = element_text(colour="#3a86ff"),
        plot.background = element_blank()) +
  guides(size = guide_legend(nrow = 1)) 
  
yearly_plot

## Export as png
png(filename = "output/figures/fig1a.png",width = 750, height = 250)
yearly_plot
dev.off()

## Export as pdf
pdf(file = "output/figures/fig1a.pdf",width = 12, height = 4)
yearly_plot
dev.off()
